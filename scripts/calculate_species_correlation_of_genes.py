#!/usr/bin/env python3

import sys

import pandas as pd
import xarray as xr
from scipy.spatial.distance import cdist
from lib.pandas_util import align_indexes, idxwhere
from lib.util import info
from scipy.stats import trim_mean

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]
    species_id = sys.argv[2]
    gene_depth_inpath = sys.argv[3]
    # gene_meta_inpath = sys.argv[4]
    # aggregate_genes_by = sys.argv[5]
    transformation_root = float(sys.argv[4])
    corr_thresh = float(sys.argv[5])
    trim_frac = float(sys.argv[6])
    corr_outpath = sys.argv[7]
    sample_depth_outpath = sys.argv[8]
    gene_depth_outpath = sys.argv[9]

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
    gene_depth = xr.load_dataarray(gene_depth_inpath)

    info("Transforming input data.")
    info("Aligning indexes.")
    sample_list = list(
        set(gene_depth.sample.values) & set(species_depth.sample.values)
    )
    species_depth = species_depth.sel(sample=sample_list)
    gene_depth = gene_depth.sel(sample=sample_list)

    info("Calculating correlation.")
    trnsfm = lambda x: x ** (1 / transformation_root)
    corr = pd.Series(
        1
        - cdist(
            trnsfm(species_depth.expand_dims(dict(_=1))),
            trnsfm(gene_depth),
            metric="cosine",
        )[0],
        index=gene_depth.gene_id,
    )

    info(f"Identifying highly correlated genes (> {corr_thresh}).")
    gene_hits = idxwhere(corr > corr_thresh)
    nhits = len(gene_hits)
    info(f"Found {nhits} highly correlated genes.")
    info("Calculating mean depth of samples.")
    sample_mean_depth = xr.apply_ufunc(
        trim_mean,
        gene_depth.sel(gene_id=gene_hits),
        kwargs=dict(proportiontocut=trim_frac),
        input_core_dims=[["gene_id", "sample"]],
        output_core_dims=[["sample"]],
    )
    info("Calculating mean depth of genes.")
    gene_mean_depth = gene_depth.sum("sample") / sample_mean_depth.sum()

    info("Writing correlation.")
    corr.to_csv(corr_outpath, sep="\t", header=False)
    info("Writing sample mean depth.")
    sample_mean_depth.to_series().to_csv(sample_depth_outpath, sep="\t", header=False)
    info("Writing gene mean depth.")
    gene_mean_depth.to_series().to_csv(gene_depth_outpath, sep="\t", header=False)
