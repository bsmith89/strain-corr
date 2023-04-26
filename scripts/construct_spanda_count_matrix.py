#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    outpath = sys.argv[1]
    sample_args = sys.argv[2:]

    data = {}
    for arg in sample_args:
        sample_name, sample_path = arg.split("=")
        data[sample_name] = pd.read_table(
            sample_path, names=["gene_id", "tally"], index_col="gene_id"
        ).tally
    data = pd.DataFrame(data).fillna(0).rename_axis(index="Gene")

    data.to_csv(outpath)
