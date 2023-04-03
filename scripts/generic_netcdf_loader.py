#!/usr/bin/env python3

import sys
import pandas as pd


if __name__ == "__main__":
    index_cols = sys.argv[1].split(",")
    outpath = sys.argv[2]
    data = pd.read_table(sys.stdin, index_col=index_cols)
    data = data.squeeze()
    data = data.to_xarray()
    data.to_netcdf(outpath)
