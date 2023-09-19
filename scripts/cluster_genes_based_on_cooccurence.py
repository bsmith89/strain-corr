#!/usr/bin/env python3

from fastcluster import linkage
from scipy.cluster.hierarchy import fcluster
import pandas as pd
from lib.pandas_util import idxwhere
import sys

if __name__ == "__main__":
    gene_path = sys.argv[1]
    filt_path = sys.argv[2]
    thresh = float(sys.argv[3])
    outpath = sys.argv[4]

    filt_list = idxwhere(
        pd.read_table(filt_path, index_col="genome_id")
        .rename(str)
        .passes_filter.astype(bool)
    )
    data = pd.read_table(gene_path, index_col="gene_id").astype(bool).loc[:, filt_list]

    drop_nohit_genes_list = idxwhere(data.sum(1) == 0)
    drop_ubiq_genes_list = idxwhere((~data).sum(1) == 0)
    drop_single_hit_genes_list = idxwhere(data.sum(1) == 1)
    drop_only_one_missing_genes_list = idxwhere((~data).sum(1) == 1)
    data_filt = data.drop(
        drop_nohit_genes_list
        + drop_ubiq_genes_list
        + drop_single_hit_genes_list
        + drop_only_one_missing_genes_list
    )

    gene_linkage = linkage(data_filt, metric="correlation", method="average")
    cluster = pd.Series(
        fcluster(gene_linkage, thresh, criterion="distance"), index=data_filt.index
    )
    cluster = pd.concat(
        [
            cluster,
            pd.Series(-1, index=drop_ubiq_genes_list),
            pd.Series(-2, index=drop_only_one_missing_genes_list),
            pd.Series(-3, index=drop_single_hit_genes_list),
            pd.Series(-4, index=drop_nohit_genes_list),
        ]
    )

    cluster.to_csv(outpath, sep="\t", header=False)
