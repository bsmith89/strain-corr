#!/usr/bin/env python3

import pandas as pd
import sys
from lib.pandas_util import align_indexes


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
