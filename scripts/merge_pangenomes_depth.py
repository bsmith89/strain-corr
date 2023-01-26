#!/usr/bin/env python3

# import pandas as pd
import sys
import subprocess
import numpy as np
from itertools import product
import xarray as xr
from tqdm import tqdm
from lib.util import info


def read_table(file):
    out = {}
    _ = next(file)  # Skip header line.
    for line in file:
        gene_id, depth = line.strip().split(b"\t")
        out[gene_id.decode("utf-8")] = float(depth)
    return out


def compile_dataarray_from_dict_of_dicts(data):
    info("Merging indexes.")
    row_idx = list(data.keys())
    col_idx = set()
    for r in row_idx:
        col_idx |= set(data[r].keys())
    col_idx = list(col_idx)

    info("Entering data into matrix.")
    out = np.zeros((len(row_idx), len(col_idx)))
    for (i, r), (j, c) in tqdm(
        product(enumerate(row_idx), enumerate(col_idx)),
        total=len(row_idx) * len(col_idx),
        ncols=0,
    ):
        row = data[r]
        if c in row:
            out[i, j] = row[c]

    info("Wrapping data structure.")
    out = xr.DataArray(out, coords=[row_idx, col_idx], dims=["sample", "gene_id"])
    return out


if __name__ == "__main__":
    outpath = sys.argv[1]
    num_inputs = len(sys.argv[2:])

    data = {}
    info(f"Loading data from {num_inputs} files.")
    for arg in tqdm(sys.argv[2:], ncols=0):
        sample, path = arg.split("=")
        with subprocess.Popen(["bzip2", "-dc", path], stdout=subprocess.PIPE) as proc:
            data[sample] = read_table(proc.stdout)
    info("Compiling inputs into a uniform matrix.")
    data = compile_dataarray_from_dict_of_dicts(data)
    info(f"Writing output file: {outpath}")
    data.to_netcdf(outpath)
