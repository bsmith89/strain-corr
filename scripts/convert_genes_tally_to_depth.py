#!/usr/bin/env python3

import sys
import pandas as pd
from lib.util import info

if __name__ == "__main__":
    tallypath = sys.argv[1]
    metapath = sys.argv[2]
    assumed_read_length = int(sys.argv[3])
    outpath = sys.argv[4]

    info("Loading input data.")
    info("Loading gene tallies.")
    tally = (
        pd.read_table(
            tallypath,
            index_col="gene_id",
        )
        .astype(int)
        .rename_axis(columns="sample", index="gene_id")
        .stack()
        .to_xarray()
    )
    info("Loading gene metadata.")
    meta = (
        pd.read_table(metapath)
        .set_index("centroid_99", drop=False)
        .rename_axis(index="gene_id")
        .loc[tally.gene_id]
    )
    gene_length = meta.centroid_99_length.loc[tally.gene_id].to_xarray()
    info("Calculating gene depth.")
    depth = tally * assumed_read_length / gene_length
    info("Writing output.")
    depth.to_dataset(name="depth").to_netcdf(
        outpath, encoding=dict(depth=dict(zlib=True, complevel=5))
    )
