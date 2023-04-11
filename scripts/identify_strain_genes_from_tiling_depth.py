#!/usr/bin/env python3

import xarray as xr
import sys
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    depth_inpath = sys.argv[1]
    strain_id = sys.argv[2]
    depth_thresh = float(sys.argv[3])
    outpath = sys.argv[4]

    depth = xr.load_dataarray(depth_inpath).sel(sample=sys.argv[2]).to_series()
    gene_hits = idxwhere(depth > depth_thresh)
    with open(outpath, "w") as f:
        for i, g in enumerate(gene_hits):
            print(i, g, sep="\t", file=f)
