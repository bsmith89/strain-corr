#!/usr/bin/env python3

import xarray as xr
import sys
from lib.pandas_util import idxwhere
import pandas as pd


if __name__ == "__main__":
    depth_inpath = sys.argv[1]
    strain_id = sys.argv[2]
    depth_thresh = float(sys.argv[3])
    outpath = sys.argv[4]

    depth = xr.load_dataarray(depth_inpath).sel(sample=strain_id).to_series()
    gene_hits = idxwhere(depth > depth_thresh)
    output = pd.DataFrame(1, index=gene_hits, columns=[strain_id]).rename_axis(index="gene_id")
    output.to_csv(outpath, sep='\t')
