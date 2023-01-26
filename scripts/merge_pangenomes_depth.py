#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info
import sqlite3

if __name__ == "__main__":
    temp_db = sys.argv[1]
    outpath = sys.argv[2]
    num_inputs = len(sys.argv[3:])

    info("Loading data from db.")
    con = sqlite3.connect(temp_db)
    data = (
        pd.read_sql(
            "SELECT sample, gene_id, depth FROM main ORDER BY sample, gene_id;",
            con=con,
            index_col=["sample", "gene_id"],
        )
        .squeeze()
        .to_xarray()
        .fillna(0)
    )
    info(f"Writing output file: {outpath}")
    data.to_netcdf(outpath)
