#!/usr/bin/env python3

import pandas as pd
from lib.pandas_util import read_table_lz4
import sys


EMPTY_DATA = pd.Series([], name='tally').rename_axis("gene_id")


if __name__ == "__main__":
    outpath = sys.argv[1]
    sample_args = sys.argv[2:]

    data = {}
    for arg in sample_args:
        sample_name, sample_path = arg.split("=")
        try:
            data[sample_name] = read_table_lz4(sample_path, index_col="gene_id").tally
        except pd.errors.EmptyDataError:
            data[sample_name] = EMPTY_DATA
    data = pd.DataFrame(data).fillna(0).astype(int).rename_axis(index="Gene")

    data.to_csv(outpath)
