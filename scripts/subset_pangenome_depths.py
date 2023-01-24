#!/usr/bin/env python3

import xarray as xr
import sys


if __name__ == "__main__":
    genes_path = sys.argv[1]
    depth_path = sys.argv[2]
    out_path = sys.argv[3]

    with open(genes_path) as f:
        select_genes = [line.strip() for line in f]

    depth = xr.open_dataarray(depth_path)
    select_genes = set(select_genes) & set(depth.gene_id.values)

    depth = depth.sel(gene_id=list(select_genes))
    depth.to_netcdf(out_path)
