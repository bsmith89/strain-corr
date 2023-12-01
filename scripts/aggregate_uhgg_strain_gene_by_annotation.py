#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    uhgg_inpath = sys.argv[1]
    agg_inpath = sys.argv[2]
    agg_column = sys.argv[3]
    outpath = sys.argv[4]

    gene_table = pd.read_table(uhgg_inpath, index_col="gene_id")
    agg = pd.read_table(agg_inpath, index_col="gene_id").squeeze()
    result = (
        gene_table.astype(bool).join(agg).groupby(agg_column).any().astype(int)
    ).rename_axis(index="gene_id")
    result.to_csv(outpath, sep="\t")
