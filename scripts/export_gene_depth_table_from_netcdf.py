#!/usr/bin/env python3

import sys
import xarray as xr

if __name__ == "__main__":
    inpath = sys.argv[1]
    outpath = sys.argv[2]

    data = xr.load_dataarray(inpath)
    data = data.to_pandas().T.rename_axis(index="gene")
    data.to_csv(outpath, sep='\t')
