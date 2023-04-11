#!/usr/bin/env python3
#
# USAGE: {input.script} \
#     {input.species_depth} \
#     {wildcards.species} \
#     {params.diss_thresh} \
#     {input.gene_depth} \
#     {params.trnsf_root} \
#     {params.n_marker_genes} \
#     {output.species_gene} \
#     {output.species_corr} \

import sys
import pandas as pd
import xarray as xr
from lib.pandas_util import idxwhere
from lib.util import info
from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.cluster import AgglomerativeClustering

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]  # {input.species_depth}
    species_id = sys.argv[2]  # {wildcards.species}
    diss_thresh = float(sys.argv[3])  # {params.diss_thresh}
    gene_depth_inpath = sys.argv[4]  # {input.gene_depth}
    root = float(sys.argv[5])  # {params.trnsf_root}
    n_marker_genes = int(sys.argv[6])  # {params.n_marker_genes}
    species_gene_outpath = sys.argv[7]  # {output.species_gene}
    species_corr_outpath = sys.argv[8]  # {output.species_corr}

    info("Loading input data.")
    info("Loading species depth.")
    species_depth = (
        pd.read_table(
            species_depth_inpath,
            dtype={"sample": str, "species_id": str, "depth": float},
            index_col=["sample", "species_id"],
        )
        .squeeze()
        .to_xarray()
        .sel(species_id=species_id)
        .fillna(0)
    )
    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_inpath).fillna(
        0
    )  # FIXME: Shouldn't be necessary.

    info("Aligning indexes.")
    sample_list = list(set(gene_depth.sample.values) & set(species_depth.sample.values))
    species_depth = species_depth.sel(sample=sample_list)
    gene_depth = gene_depth.sel(sample=sample_list)

    n_samples = len(gene_depth.sample)
    info(f"Calculating cosine dissimilarity of gene depths among {n_samples} samples.")
    cdmat = pdist(gene_depth.to_pandas(), metric="cosine")
    info(
        f"De-replicating samples using agglomorative clustering with a cosine dissimilarity threshold of {diss_thresh}."
    )
    clust = (
        pd.Series(
            AgglomerativeClustering(
                n_clusters=None,
                affinity="precomputed",
                linkage="average",
                distance_threshold=diss_thresh,
            ).fit_predict(squareform(cdmat)),
            index=gene_depth["sample"],
            name="clust",
        )
        .rename_axis(index="sample")
        .to_xarray()
    )
    species_depth = species_depth.groupby(clust).sum()
    gene_depth = gene_depth.groupby(clust).sum()
    n_clusts = len(gene_depth.clust)
    info(f"Found {n_clusts} de-replicated clusters.")

    info("Transforming input data.")
    species_depth = species_depth.pipe(lambda x: x ** (1 / root))
    gene_depth = gene_depth.pipe(lambda x: x ** (1 / root))

    n_genes = len(gene_depth.gene_id)
    info(f"Calculating correlation of {n_genes} genes.")
    corr = pd.Series(
        1
        - cdist(
            species_depth.expand_dims(dict(_=1)),
            gene_depth.transpose("gene_id", "clust"),
            metric="cosine",  # lambda x, y: sp.stats.spearmanr(x, y)[0],
        )[0],
        index=gene_depth.gene_id,
    )

    info(f"Identifying top {n_marker_genes} highly correlated genes.")
    corr_thresh = corr.sort_values(ascending=False).head(n_marker_genes).min()
    species_gene_hit = idxwhere(corr >= corr_thresh)
    nhits = len(species_gene_hit)
    info(f"Found {nhits} highly correlated genes (correlation > {corr_thresh}).")

    info("Writing species gene list.")
    with open(species_gene_outpath, "w") as f:
        for g in species_gene_hit:
            print(g, file=f)
    info("Writing species correlations.")
    corr.to_csv(species_corr_outpath, sep="\t", header=False)
