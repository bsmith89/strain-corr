#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import sys
from lib.util import info
from lib.pandas_util import idxwhere
from scipy.spatial.distance import cdist
from tqdm import tqdm


if __name__ == "__main__":
    all_species_depth_inpath = sys.argv[1]
    species_id = sys.argv[2]
    species_depth_thresh_abs = float(sys.argv[3])
    species_depth_thresh_pres = float(sys.argv[4])
    strain_frac_inpath = sys.argv[5]
    strain_frac_thresh = float(sys.argv[6])
    gene_depth_inpath = sys.argv[7]
    transformation_root = float(sys.argv[8])
    outpath = sys.argv[9]

    info("Loading input data.")
    info("Loading species depth.")
    all_species_depth = (
        pd.read_table(
            all_species_depth_inpath,
            dtype={"sample": str, "species_id": str, "depth": float},
            index_col=["sample", "species_id"],
        )
        .squeeze()
        .to_xarray()
        .fillna(0)
    )
    info("Loading strain fractions.")
    strain_frac = (
        pd.read_table(strain_frac_inpath, index_col=["sample", "strain"])
        .squeeze()
        .to_xarray()
    )
    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_inpath)
    info("Aligning indexes.")
    shared_samples = list(
        set(all_species_depth.sample.values)
        & set(strain_frac.sample.values)
        & set(gene_depth.sample.values)
    )
    all_species_depth = all_species_depth.sel(sample=shared_samples)
    strain_frac = strain_frac.sel(sample=shared_samples)
    gene_depth = gene_depth.sel(sample=shared_samples)

    info("Identifying strain-pure samples.")
    homogenous_samples = idxwhere(
        (strain_frac > strain_frac_thresh).any("strain").to_series()
        & (
            all_species_depth.sel(species_id=species_id) > species_depth_thresh_pres
        ).to_series()
    )
    no_species_samples = idxwhere(
        all_species_depth.sel(species_id=species_id).to_series()
        < species_depth_thresh_abs
    )
    # strain_total_depth = (
    #     strain_frac.sel(sample=homogenous_samples)
    #     * all_species_depth.sel(species_id=species_id, sample=homogenous_samples)
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
    info("Transforming depths.")
    trnsfm = lambda x: x ** (1 / transformation_root)
    gene_depth = trnsfm(gene_depth)
    all_species_depth = trnsfm(all_species_depth)
    info("Iterating strains.")
    corr = {}
    for strain in tqdm(strain_frac.strain.to_series().to_list()):
        strain_pure_samples = idxwhere(
            strain_frac.sel(strain=strain).to_series() > strain_frac_thresh
        )
        focal_samples = strain_pure_samples + no_species_samples
        y = gene_depth.sel(sample=focal_samples)
        x = all_species_depth.sel(sample=focal_samples)
        corr[strain] = (
            pd.DataFrame(
                1
                - cdist(
                    x.T,
                    y,
                    metric="cosine",
                ),
                columns=gene_depth.gene_id,
                index=all_species_depth.species_id,
            )
            .fillna(0)
            .stack()
        )
    info("Compiling correlations table.")
    corr = (
        pd.DataFrame(corr)
        .rename_axis(index=["species_id", "gene_id"], columns="strain")
        .stack()
        .rename("correlation")
    )
    corr.to_xarray().to_dataset(name="correlation").to_netcdf(
        outpath, encoding=dict(correlation=dict(zlib=True, complevel=5))
    )
