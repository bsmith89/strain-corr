#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    gene_info_path = sys.argv[1]
    agg_by = sys.argv[2]
    hits_path = sys.argv[3]
    outpath = sys.argv[4]

    agg_col = pd.read_table(gene_info_path, index_col="gene_id")[agg_by]
    hits = pd.read_table(hits_path, names=["orf_id", "gene_id"])
    hits = hits.assign(agg=lambda x: x.gene_id.map(agg_col))
    hits = hits.drop(columns=["gene_id"]).drop_duplicates()
    hits = hits.rename(columns={"agg": "gene_id"})
    hits = hits[["gene_id", "orf_id"]]
    hits.to_csv(outpath, sep="\t", index=False)
