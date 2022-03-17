#!/usr/bin/env python3
"""Merge r1 and r2 reads from previously concatenated libraries
for a single species and then write the result to compressed NetCDF
file.

Usage:

```bash
merge_gtpro_to_netcdf.py <R1> <R2> <OUTFILE>
```

Then use xarray from within a python script to read this data quickly:

```python
xarray.open_dataarray('<OUTFILE>')
```

"""

import pandas as pd
import sys
import numpy as np
from lib.util import info
import warnings
import math
from functools import reduce

DEFAULT_COMPRESSION_LEVEL = 6
INT_TYPES = ["uint8", "uint16", "uint32", "uint64"]


def _load_read(path):
    return (
        pd.read_table(
            path,
            names=[
                "library_id",
                "species_id",
                "position",
                "_3",
                "_4",
                "_5",
                "_6",
                "ref",
                "alt",
            ],
            index_col=["library_id", "species_id", "position"],
        )[["ref", "alt"]]
        .rename_axis(columns="allele")
        .stack()
        .astype(int)
        .squeeze()
    )


def _write_netcdf(path, data, compression_level=DEFAULT_COMPRESSION_LEVEL):
    (
        data.to_dataset(name="tally").to_netcdf(
            path,
            encoding=dict(tally=dict(zlib=True, complevel=compression_level)),
        )
    )


def _concat_reads_to_xr(r1, r2):
    r1 = r1.to_frame("r1").rename_axis(columns="read").stack()
    r2 = r2.to_frame("r2").rename_axis(columns="read").stack()
    return pd.concat([r1, r2]).to_xarray().fillna(0).astype(int)


def is_empty(x):
    mul = lambda x, y: x * y
    return reduce(mul, x.shape, 1) == 0


if __name__ == "__main__":
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        module="numpy.lib.arraysetops",
        lineno=580,
    )

    r1path = sys.argv[1]
    r2path = sys.argv[2]
    outpath = sys.argv[3]

    info("Loading r1")
    r1 = _load_read(r1path)
    info("Loading r2")
    r2 = _load_read(r2path)

    info("Concatenating reads")
    data = _concat_reads_to_xr(r1, r2)
    info("Summing r1 and r2")
    data = data.sum("read")
    if is_empty(data):
        info("Data is empty. Will write an empty dataarray.")
        max_value = 0
    else:
        max_value = data.max().values
    info(f"Checking smallest int type that can store maximum value: {max_value}")
    for dtype in INT_TYPES:
        if max_value < np.iinfo(dtype).max:
            break
    else:
        raise ValueError(f"Largest count is too big for all dtypes checked (final: {dtype})")
    info(f"Casting to {dtype}")
    data = data.astype(dtype)

    info("Writing NetCDF")
    _write_netcdf(outpath, data)
    info("Finished")
