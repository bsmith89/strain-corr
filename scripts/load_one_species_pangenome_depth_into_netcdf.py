#!/usr/bin/env python3

import pandas as pd
import sys
import sqlite3


if __name__ == "__main__":
    dbpath = sys.argv[1]
    species = sys.argv[2]
    outpath = sys.argv[3]

    con = sqlite3.connect(f"file:{dbpath}?mode=ro", uri=True)
    con.executescript(
        """
PRAGMA journal_mode=OFF;
PRAGMA synchronous=OFF;
PRAGMA locking_mode=EXCLUSIVE;
PRAGMA temp_store=MEMORY;
PRAGMA cache_size=1000000;
"""
    )

    data = pd.read_sql(
        (
            f"SELECT sample, gene_id, depth "
            f"FROM gene JOIN sample_x_gene USING (gene_id) "
            f"WHERE species = '{species}';"
        ),
        con=con,
        index_col=["sample", "gene_id"],
    ).squeeze()
    data.to_xarray().fillna(0).to_netcdf(outpath)
