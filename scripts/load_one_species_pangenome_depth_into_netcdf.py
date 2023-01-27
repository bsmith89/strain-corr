#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    outpath = sys.argv[1]
    data = pd.read_table(sys.stdin, index_col=["sample", "gene_id"]).squeeze()
    data.to_xarray().to_netcdf(outpath)
