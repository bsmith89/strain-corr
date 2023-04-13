#!/usr/bin/env python3

import pandas as pd
import sys
import scipy as sp
import scipy.stats
from lib.pandas_util import idxwhere
from tqdm import tqdm
import numpy as np


def calculate_asymetric_statistics(match_matrix):
    # TODO: Inline these calculations and drop approach to "asymetric" matching.
    tp_orf = idxwhere(match_matrix.any(1))
    fn_orf = idxwhere(~match_matrix.any(1))
    tp_gene = idxwhere(match_matrix.any(0))
    fp_gene = idxwhere(~match_matrix.any(0))

    n_tp_gene = len(tp_gene)
    n_fp_gene = len(fp_gene)
    n_fn_orf = len(fn_orf)
    n_tp_orf = len(tp_orf)

    if (n_tp_gene + n_fp_gene) == 0:
        precision = 0  # np.nan
    else:
        precision = n_tp_gene / (n_tp_gene + n_fp_gene)
    if (n_tp_orf + n_fn_orf) == 0:
        recall = 0  # np.nan
    else:
        recall = n_tp_orf / (n_tp_orf + n_fn_orf)

    if np.isnan(precision) or np.isnan(recall):
        f1 = 0  # np.nan
    else:
        f1 = sp.stats.hmean([precision, recall])

    return precision, recall, f1


if __name__ == "__main__":
    gene_matching_path = sys.argv[1]
    corr_path = sys.argv[2]
    depth_path = sys.argv[3]
    thresh_path = sys.argv[4]
    outpath = sys.argv[5]

    # Load data
    # TODO: Change "matching gene" input to just a list.
    matching_gene = (
        pd.read_table(
            gene_matching_path, names=["orf", "gene"], index_col=["orf", "gene"]
        )
        .assign(hit=True)
        .unstack("orf", fill_value=False)
    )
    corr = pd.read_table(corr_path, index_col=["gene_id", "strain"]).squeeze().unstack()
    depth = (
        pd.read_table(depth_path, index_col=["gene_id", "strain"]).squeeze().unstack()
    )
    thresh = pd.read_table(thresh_path, index_col=["strain"])

    # Align data
    gene_list = list(set(corr.index) | set(depth.index) | set(matching_gene.index))
    depth = depth.reindex(index=gene_list, fill_value=0)
    corr = corr.reindex(index=gene_list, fill_value=0)
    matching_gene = matching_gene.reindex(gene_list, fill_value=False).T

    # Align strains
    strain_list = list(set(corr.columns) & set(depth.columns))
    depth = depth[strain_list]
    corr = corr[strain_list]

    results = {}
    for strain in tqdm(corr.columns):
        # Classify genes
        depth_and_corr_hit = idxwhere(
            (corr[strain] >= thresh.correlation[strain])
            & (depth[strain] >= thresh.depth_low[strain])
        )

        precision, recall, f1 = calculate_asymetric_statistics(
            matching_gene.loc[:, depth_and_corr_hit]
        )
        results[strain] = [
            precision,
            recall,
            f1,
        ]

    # Compile results
    results = (
        pd.DataFrame(
            results,
            index=[
                "precision",
                "recall",
                "f1",
            ],
        )
        .T.rename_axis(index="strain")
        .sort_values("f1", ascending=False)
    )

    # Write output
    results.to_csv(outpath, sep="\t")
