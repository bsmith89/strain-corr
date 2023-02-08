#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import sys
from lib.util import info
from lib.pandas_util import idxwhere
from tqdm import tqdm
import numpy as np

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]
    strain_frac_inpath = sys.argv[2]
    strain_frac_thresh = float(sys.argv[3])
    gene_depth_inpath = sys.argv[4]
    outpath = sys.argv[5]

    info("Loading data.")
    gene_depth = xr.load_dataarray(gene_depth_inpath)
    strain_frac = (
        pd.read_table(strain_frac_inpath, index_col=["sample", "strain"])
        .squeeze()
        .to_xarray()
    )
    species_depth = (
        pd.read_table(
            species_depth_inpath, names=["sample", "depth"], index_col="sample"
        )
        .squeeze()
        .to_xarray()
    )
    info("Aligning datasets.")
    shared_samples = list(
        set(species_depth.sample.values)
        & set(strain_frac.sample.values)
        & set(gene_depth.sample.values)
    )
    species_depth = species_depth.sel(sample=shared_samples)
    strain_frac = strain_frac.sel(sample=shared_samples)
    gene_depth = gene_depth.sel(sample=shared_samples)

    info("Identifying strain-pure samples.")
    homogenous_samples = idxwhere(
        (strain_frac > strain_frac_thresh).any("strain").to_series()
    )
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
    #
    info("Iterating strains.")
    depth_ratio = {}
    for strain in tqdm(strain_frac.strain.to_series().to_list()):
        strain_pure_samples = idxwhere(
            (strain_frac.sel(strain=strain) > strain_frac_thresh).to_series()
        )
        if strain_pure_samples:
            depth_ratio[strain] = (
                gene_depth.sel(sample=strain_pure_samples).sum("sample")
                / species_depth.sel(sample=strain_pure_samples).sum()
            ).to_series()
        else:
            depth_ratio[strain] = pd.Series(np.nan, index=gene_depth.gene_id)

    info("Compiling depth ratios table.")
    depth_ratio = (
        pd.DataFrame(depth_ratio)
        .rename_axis(columns="strain", index="gene_id")
        .stack()
        .rename("depth_ratio")
    )
    info("Writing output.")
    depth_ratio.to_csv(outpath, sep="\t")
