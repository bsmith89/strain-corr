#!/usr/bin/env python3
# {input.script} {input.norm} {input.depth} {output}

import xarray as xr
import pandas as pd
import sys


if __name__ == "__main__":
    norm_path = sys.argv[1]
    depth_path = sys.argv[2]
    out_path = sys.argv[3]

    norm = (
        pd.read_table(norm_path, names=["sample", "norm"], index_col="sample")
        .squeeze()
        .to_xarray()
    )
    depth = xr.load_dataarray(depth_path)
    (depth / norm).to_netcdf(out_path)
