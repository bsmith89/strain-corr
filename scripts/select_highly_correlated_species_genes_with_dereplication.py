#!/usr/bin/env python3
#
# USAGE: {input.script} \
#     {input.species_depth} \
#     {wildcards.species} \
#     {params.absent_thresh} \
#     {params.present_thresh} \
#     {input.strain_frac} \
#     {params.frac_thresh} \
#     {input.gene_depth} \
#     {params.n_marker_genes} \
#     {output.species_gene} \
#     {output.species_corr} \

import sys
import pandas as pd
import xarray as xr
from scipy.spatial.distance import cdist
from lib.pandas_util import idxwhere
from lib.util import info

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]  # {input.species_depth}
    species_id = sys.argv[2]  # {wildcards.species}
    absent_thresh = float(sys.argv[3])  # {params.absent_thresh}
    present_thresh = float(sys.argv[4])  # {params.present_thresh}
    strain_frac_inpath = sys.argv[5]  # {input.strain_frac}
    frac_thresh = float(sys.argv[6])  # {params.frac_thresh}
    gene_depth_inpath = sys.argv[7]  # {input.gene_depth}
    root = float(sys.argv[8])  # {params.trnsf_root}
    n_marker_genes = int(sys.argv[9])  # {params.n_marker_genes}
    species_gene_outpath = sys.argv[10]  # {output.species_gene}
    species_corr_outpath = sys.argv[11]  # {output.species_corr}

    info("Loading input data.")
    info("Loading species depth.")
    species_depth = (
        pd.read_table(
            species_depth_inpath,
            dtype={"species_id": str},
            index_col=["sample", "species_id"],
        )
        .squeeze()
        .unstack(fill_value=0)[species_id]
        .to_xarray()
    )

    info("Loading strain composition.")
    strain_frac = pd.read_table(
        strain_frac_inpath, index_col=["sample", "strain"]
    ).community.to_xarray()

    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_inpath).fillna(0)

    info("Identifying strain-pure samples.")
    sample_to_strain = (
        (strain_frac >= frac_thresh)
        .to_series()[lambda x: x]
        .index.to_frame()[[]]
        .reset_index("strain")
        .squeeze()
        .to_xarray()
    )

    info("Aligning indexes.")
    sample_list = list(set(gene_depth.sample.values) & set(species_depth.sample.values))
    species_depth = species_depth.sel(sample=sample_list)
    gene_depth = gene_depth.sel(sample=sample_list)
    sample_to_strain = sample_to_strain.sel(
        sample=list(set(sample_to_strain.sample.values) & set(sample_list))
    )
    nospecies_samples = idxwhere(species_depth.to_series() < absent_thresh)

    info("De-replicating strain-pure samples.")
    species_depth_agg = (
        species_depth.sel(sample=sample_to_strain.sample)
        .groupby(sample_to_strain)
        .sum()
    )
    gene_depth_agg = (
        gene_depth.sel(sample=sample_to_strain.sample).groupby(sample_to_strain).sum()
    )

    if species_depth_agg.sizes["strain"] == 1:
        info(
            "WARNING: Only one strain group remaining after de-replication; this will skew the results."
        )

    info("Concatenating strain samples with no-species samples.")
    species_depth = xr.concat(
        [
            species_depth_agg.rename({"strain": "sample"}),
            species_depth.sel(sample=nospecies_samples),
        ],
        dim="sample",
    )
    gene_depth = xr.concat(
        [
            gene_depth_agg.rename({"strain": "sample"}),
            gene_depth.sel(sample=nospecies_samples),
        ],
        dim="sample",
    )

    info(f"Transforming depth data. (root-{root}).")
    trnsf = lambda x: x ** (1 / root)
    species_depth = species_depth.pipe(trnsf)
    gene_depth = gene_depth.pipe(trnsf)

    n_genes = len(gene_depth.gene_id)
    n_samples = len(gene_depth.sample)
    info(f"Calculating correlation among {n_genes} genes across {n_samples} samples.")
    corr = pd.Series(
        1
        - cdist(
            species_depth.expand_dims(dict(_=1)),
            gene_depth.transpose("gene_id", "sample"),
            metric="cosine",
        )[0],
        index=gene_depth.gene_id,
    )

    info(f"Identifying top {n_marker_genes} highly correlated genes.")
    corr_thresh = corr.sort_values(ascending=False).head(n_marker_genes).min()
    species_gene_hit = idxwhere(corr >= corr_thresh)
    nhits = len(species_gene_hit)
    info(f"Found {nhits} highly correlated genes (cosine similarity > {corr_thresh}).")

    info("Writing species gene list.")
    with open(species_gene_outpath, "w") as f:
        for g in species_gene_hit:
            print(g, file=f)
    info("Writing species correlations.")
    corr.to_csv(species_corr_outpath, sep="\t", header=False)
