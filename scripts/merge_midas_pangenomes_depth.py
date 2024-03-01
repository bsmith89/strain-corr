#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info
from lib.pandas_util import read_table_lz4
from tqdm import tqdm


if __name__ == "__main__":
    outpath = sys.argv[1]
    sample_paths = dict((arg.split("=", maxsplit=1) for arg in sys.argv[2:]))

    results = []
    for i, (sample, path) in tqdm(enumerate(sample_paths.items())):
        results.append(
            read_table_lz4(path)[["mean_depth", "cluster_75_id"]]
            .assign(sample=sample)
            .set_index(["sample", "cluster_75_id"])
            .mean_depth
        )
    results = pd.concat(results)

    info("Writing output.")
    results = results.to_xarray()
    breakpoint()
    results.to_netcdf(outpath)
