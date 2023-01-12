#!/usr/bin/env python3

import pandas as pd
import sys
import subprocess

if __name__ == "__main__":
    outpath = sys.argv[1]

    data = []
    for arg in sys.argv[2:]:
        sample, path = arg.split("=")
        with subprocess.Popen(["lz4", "-dc", path], stdout=subprocess.PIPE) as proc:
            data.append(
                pd.read_table(proc.stdout).assign(
                    sample=sample, depth=lambda x: x.mean_coverage * x.fraction_covered
                )[["sample", "gene_id", "depth"]]
            )
    data = pd.concat(data).set_index(["sample", "gene_id"]).squeeze().to_xarray().fillna(0)
    data.to_netcdf(outpath)
