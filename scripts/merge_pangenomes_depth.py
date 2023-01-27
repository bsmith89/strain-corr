#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info

if __name__ == "__main__":
    inpath = sys.argv[1]
    outpath = sys.argv[2]

    info("Loading data.")
    data = pd.read_table(
        inpath,
        names=["sample", "gene_id", "depth"],
        index_col=["sample", "gene_id"],
    )
    info("Reshaping data.")
    data = data.squeeze().to_xarray().fillna(0)
    info(f"Writing output file: {outpath}")
    data.to_netcdf(outpath)
