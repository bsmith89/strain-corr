#!/usr/bin/env python3

import sys

import pandas as pd
import xarray as xr
from scipy.spatial.distance import cdist
from lib.pandas_util import idxwhere
from lib.util import info
from scipy.stats import trim_mean

if __name__ == "__main__":
    species_depth_inpath = sys.argv[1]
    species_id = sys.argv[2]
    gene_depth_inpath = sys.argv[3]
    # gene_meta_inpath = sys.argv[4]
    # aggregate_genes_by = sys.argv[5]
    transformation_exponent = float(sys.argv[4])
    n_marker_genes = int(sys.argv[5])
    # corr_thresh = float(sys.argv[5])
    trim_frac = float(sys.argv[6])
    corr_outpath = sys.argv[7]
    sample_depth_outpath = sys.argv[8]
    gene_depth_outpath = sys.argv[9]
    threshold_outpath = sys.argv[10]

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
    sample_list = list(set(gene_depth.sample.values) & set(species_depth.sample.values))
    species_depth = species_depth.sel(sample=sample_list)
    gene_depth = gene_depth.sel(sample=sample_list)

    n_genes = len(gene_depth.gene_id)
    n_samples = len(gene_depth.sample)
    info(f"Calculating correlation among {n_genes} genes across {n_samples} samples.")
    trnsfm = lambda x: x ** transformation_exponent
    corr = pd.Series(
        1
        - cdist(
            trnsfm(species_depth.expand_dims(dict(_=1))),
            trnsfm(gene_depth),
            metric="cosine",
        )[0],
        index=gene_depth.gene_id,
    )

    info(f"Identifying top {n_marker_genes} highly correlated genes.")
    corr_thresh = corr.sort_values(ascending=False).head(n_marker_genes + 1).min()
    gene_hits = idxwhere(corr > corr_thresh)
    nhits = len(gene_hits)
    info(f"Found {nhits} highly correlated genes (cosine similarity > {corr_thresh}).")
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
    info("Writing species correlation threshold.")
    with open(threshold_outpath, "w") as f:
        print(corr_thresh, file=f)
