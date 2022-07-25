#!/usr/bin/env python3

import sys
import pandas as pd
from lib.util import info
import xarray as xr

if __name__ == "__main__":
    depth_inpath = sys.argv[1]
    meta_inpath = sys.argv[2]
    aggregate_genes_by = sys.argv[3]
    outpath = sys.argv[4]

    info("Loading input data.")
    info("Loading gene depth.")
    depth = xr.load_dataarray(depth_inpath)
    info("Loading gene metadata.")
    meta = (
        pd.read_table(meta_inpath)
        .set_index("centroid_99", drop=False)
        .rename_axis(index="gene_id")
        .loc[depth.gene_id]
    )
    info("Aggregating gene depth.")
    depth = depth.groupby(meta[aggregate_genes_by].to_xarray()).sum()
    info("Writing output.")
    depth.to_dataset(name="depth").to_netcdf(
        outpath, encoding=dict(depth=dict(zlib=True, complevel=5))
    )
