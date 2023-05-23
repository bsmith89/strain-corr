#!/usr/bin/env python3
# USAGE: {input.script} {input.seed} {input.n_samples} {input.samples} {output}


import sys
import pandas as pd


if __name__ == "__main__":
    n_samples = int(sys.argv[1])
    depth_inpath = sys.argv[2]
    samples_inpath = sys.argv[3]
    outpath = sys.argv[4]

    depth = pd.read_table(
        depth_inpath, names=["sample", "depth"], index_col="sample"
    ).depth
    samples = pd.read_table(samples_inpath)
    result = (
        samples.join(depth, on="sample")
        .sort_values("depth", ascending=False)
        .groupby("strain")
        .head(n_samples)
        .reset_index()
    )

    result[["strain", "sample"]].to_csv(outpath, sep="\t", index=False)
