#!/usr/bin/env python3
# USAGE:
#        {input.script} {params.trim} {input.mgen} {output}

import xarray as xr
import scipy as sp
import scipy.stats
import sys


if __name__ == "__main__":
    trim_frac = float(sys.argv[1])
    inpath = sys.argv[2]
    outpath = sys.argv[3]
    species_data = xr.load_dataarray(inpath)
    depth = (
        species_data.sum("allele")
        .to_pandas()
        .apply(lambda x: sp.stats.trim_mean(x, trim_frac), axis=1)
    ).rename_axis("sample")
    depth.to_csv(outpath, sep="\t", header=False)
