#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    data = pd.read_table(sys.argv[1])
    data = data.rename(columns={"Unnamed: 0": "gene_id"})
    data = data.set_index("gene_id").rename_axis(columns="sample")
    data = data.T.stack().rename("depth").to_xarray()
    data.to_netcdf(sys.argv[2])
