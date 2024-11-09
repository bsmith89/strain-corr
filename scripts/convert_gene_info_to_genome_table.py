#!/usr/bin/env python3
"""
Take a gene_info.txt.lz4 file and convert it into a genome-by-centroid matrix.

Assumes that input file is lz4 compressed.
Assumes that input file has a header, and is indexed by a column named gene_id.
Writes output as an xarray DataArray to NetCDF.
"""

import pandas as pd
import subprocess
import sys

if __name__ == "__main__":
    gene_path = sys.argv[1]
    column = sys.argv[2]
    out_path = sys.argv[3]

    # TODO: Drop the lz4 decompression, yeah?
    with subprocess.Popen(["lz4cat", gene_path], stdout=subprocess.PIPE) as f:
        reference_gene = pd.read_table(f.stdout, usecols=["gene_id", column])
    reference_gene = (
        reference_gene.assign(
            genome_id=lambda x: x.gene_id.str.rsplit("_", n=1).str[0]
        )[["genome_id", column]]
        .value_counts()
        .rename_axis(["genome_id", "gene_id"])
        .to_xarray()
        .fillna(0)
        .astype(int)
    )
    reference_gene.to_dataset(name="copy_number").to_netcdf(
        out_path, encoding=dict(copy_number=dict(zlib=True, complevel=5))
    )
