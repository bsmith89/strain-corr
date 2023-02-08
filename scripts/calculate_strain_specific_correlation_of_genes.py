#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import sys
from lib.util import info
from lib.pandas_util import idxwhere
from scipy.spatial.distance import cdist
from tqdm import tqdm
import numpy as np

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]
    species_depth_thresh_abs = float(sys.argv[2])
    species_depth_thresh_pres = float(sys.argv[3])
    strain_frac_inpath = sys.argv[4]
    strain_frac_thresh = float(sys.argv[5])
    gene_depth_inpath = sys.argv[6]
    transformation_exponent = float(sys.argv[7])
    outpath = sys.argv[8]

    info("Loading input data.")
    info("Loading species depth.")
    species_depth = (
        pd.read_table(
            species_depth_inpath, names=["sample", "depth"], index_col="sample"
        )
        .squeeze()
        .to_xarray()
    )
    info("Loading strain fractions.")
    strain_frac = (
        pd.read_table(strain_frac_inpath, index_col=["sample", "strain"])
        .squeeze()
        .to_xarray()
    )
    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_inpath).fillna(
        0
    )  # FIXME: Shouldn't be necessary.
    info("Aligning indexes.")
    shared_samples = list(
        set(species_depth.sample.values) & set(gene_depth.sample.values)
    )
    species_depth = species_depth.sel(sample=shared_samples)
    gene_depth = gene_depth.sel(sample=shared_samples)

    strain_frac = strain_frac.sel(
        sample=list(set(strain_frac.sample.values) & set(shared_samples))
    )

    info("Identifying strain-pure samples.")
    homogenous_samples = idxwhere(
        (strain_frac > strain_frac_thresh).any("strain").to_series()
        & (species_depth.to_series() > species_depth_thresh_pres)
    )
    no_species_samples = idxwhere(species_depth.to_series() < species_depth_thresh_abs)
    # strain_total_depth = (
    #     strain_frac.sel(sample=homogenous_samples)
    #     * species_depth.sel(sample=homogenous_samples)
    # ).sum("sample")
    # strain_sample_tally = (
    #     strain_frac.sel(sample=homogenous_samples) > strain_frac_thresh
    # ).sum("sample")
    # strain_list = idxwhere(
    #     ((strain_total_depth > 1.0) & (strain_sample_tally >= 2)).to_series()
    # )
    # strain_order = (
    #     strain_sample_tally.to_series()
    #     .loc[strain_list]
    #     .sort_values(ascending=False)
    #     .index.to_list()
    # )
    # nstrains = len(strain_order)
    # info(f"Found {nstrains} strains with total depth > 1.0 and pure in >= 2 samples.")

    info("Transforming depths.")
    trnsfm = lambda x: x**transformation_exponent
    gene_depth = trnsfm(gene_depth)
    species_depth = trnsfm(species_depth)
    info("Iterating strains.")
    corr = {}
    for strain in tqdm(strain_frac.strain.to_series().to_list()):
        strain_pure_samples = idxwhere(
            strain_frac.sel(strain=strain).to_series() > strain_frac_thresh
        )
        if strain_pure_samples:
            focal_samples = strain_pure_samples + no_species_samples
            y = gene_depth.sel(sample=focal_samples)
            x = species_depth.sel(sample=focal_samples)
            corr[strain] = pd.Series(
                1
                - cdist(
                    x.expand_dims(dict(_=1)),
                    y.transpose("gene_id", "sample"),
                    metric="cosine",
                )[0],
                index=gene_depth.gene_id,
            )
        else:
            corr[strain] = pd.Series(np.nan, index=gene_depth.gene_id)

    info("Compiling correlations table.")
    corr = (
        pd.DataFrame(corr)
        .rename_axis(index="gene_id", columns="strain")
        .stack()
        .rename("correlation")
    )
    info("Writing output.")
    corr.to_csv(outpath, sep="\t")
