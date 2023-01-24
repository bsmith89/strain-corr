#!/usr/bin/env python3

import pandas as pd
import sys
import subprocess


if __name__ == "__main__":
    outpath = sys.argv[1]

    data = []
    for arg in sys.argv[2:]:
        sample, path = arg.split("=")
        with subprocess.Popen(["bzip2", "-dc", path], stdout=subprocess.PIPE) as proc:
            data.append(
                pd.read_table(proc.stdout).assign(
                    sample=sample,
                )[["gene_id", "sample", "depth"]]
            )
    data = (
        pd.concat(data).set_index(["gene_id", "sample"]).squeeze().to_xarray().fillna(0)
    )
    data.to_netcdf(outpath)
