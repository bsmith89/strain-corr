#!/usr/bin/env python3

import pandas as pd
import sys
from lib.pandas_util import align_indexes
import warnings


if __name__ == "__main__":
    strain_corr_path = sys.argv[1]
    strain_depth_path = sys.argv[2]
    thresholds_path = sys.argv[3]
    outpath = sys.argv[4]

    strain_corr = pd.read_table(
        strain_corr_path, index_col=["gene_id", "strain"]
    ).correlation.unstack("strain")
    strain_depth = pd.read_table(
        strain_depth_path, index_col=["gene_id", "strain"]
    ).depth.unstack("strain")
    thresholds = pd.read_table(thresholds_path, index_col="strain")

    # Short-circuit on empty data:
    emptyness = {
        k: v.empty
        for k, v in {
            "strain_corr": strain_corr,
            "strain_depth": strain_depth,
            "thresholds": thresholds,
        }.items()
    }
    if any(emptyness.values()):
        warnings.warn(
            f"Writing empty output because one or more input tables were empty\n{emptyness}"
        )
        with open(outpath, "w") as f:
            print("gene_id", file=f)
        exit(0)

    # Align data:
    strain_corr, strain_depth = align_indexes(
        strain_corr, strain_depth, axis="index", how="outer", fill_value=0
    )
    thresholds_t, strain_corr, strain_depth = align_indexes(
        thresholds.T, strain_corr, strain_depth, axis="columns", how="inner"
    )
    thresholds = thresholds_t.T
    hits = (strain_corr >= thresholds.correlation) & (
        strain_depth >= thresholds.depth_low
    )

    hits[lambda x: x.sum(1) > 0].astype(int).to_csv(outpath, sep="\t")
