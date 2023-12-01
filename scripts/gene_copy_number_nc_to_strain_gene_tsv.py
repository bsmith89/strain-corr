#!/usr/bin/env python3
# "{input.script} {input.copies} {output}"

import xarray as xr
import sys

if __name__ == "__main__":
    copy_number_inpath = sys.argv[1]
    outpath = sys.argv[2]

    copy_number = xr.load_dataarray(copy_number_inpath)
    gene_presence = (copy_number >= 1).astype(int).to_pandas().T
    gene_presence.to_csv(outpath, sep='\t')
