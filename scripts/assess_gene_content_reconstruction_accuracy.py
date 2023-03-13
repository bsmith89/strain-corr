#!/usr/bin/env python3

import pandas as pd
import sys
import scipy as sp
import scipy.stats
from lib.pandas_util import idxwhere
from tqdm import tqdm


def calculate_asymetric_statistics(match_matrix):
    tp_orf = idxwhere(match_matrix.any(1))
    fn_orf = idxwhere(~match_matrix.any(1))
    tp_gene = idxwhere(match_matrix.any(0))
    fp_gene = idxwhere(~match_matrix.any(0))

    n_tp_gene = len(tp_gene)
    n_fp_gene = len(fp_gene)
    n_fn_orf = len(fn_orf)
    n_tp_orf = len(tp_orf)

    precision = n_tp_gene / (n_tp_gene + n_fp_gene)
    recall = n_tp_orf / (n_tp_orf + n_fn_orf)

    return precision, recall


if __name__ == "__main__":
    gene_matching_path = sys.argv[1]
    corr_path = sys.argv[2]
    depth_path = sys.argv[3]
    thresh_path = sys.argv[4]
    # corr_q_thresh = float(sys.argv[4])
    # depth_q_thresh = float(sys.argv[5])
    outpath = sys.argv[5]

    # Load data
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

    # Align strains FIXME: Why does this need to happen?
    # assert set(corr_q_all.columns) == set(depth_q_all.columns)
    strain_list = list(set(corr.columns) & set(depth.columns))
    depth = depth[strain_list]
    corr = corr[strain_list]

    scaled_matching_gene = (
        matching_gene.T * (matching_gene * matching_gene.sum()).sum(1)
    ).T
    gene_1to1 = idxwhere(scaled_matching_gene.sum() <= 1)
    orf_1to1 = idxwhere(scaled_matching_gene.sum(1) == 1)

    results = {}
    for strain in tqdm(corr.columns):
        # Classify genes
        depth_hit = idxwhere(depth[strain] >= thresh.depth_low[strain])
        depth_and_corr_hit = idxwhere(
            (corr[strain] >= thresh.correlation[strain])
            & (depth[strain] >= thresh.depth_low[strain])
        )

        precision, recall = calculate_asymetric_statistics(
            matching_gene.loc[:, depth_and_corr_hit]
        )
        precision_depth_only, recall_depth_only = calculate_asymetric_statistics(
            matching_gene.loc[:, depth_hit]
        )
        precision_1to1, recall_1to1 = calculate_asymetric_statistics(
            matching_gene.loc[
                orf_1to1,
                list(set(gene_1to1) & set(depth_and_corr_hit)),
            ]
        )
        (
            precision_depth_only_1to1,
            recall_depth_only_1to1,
        ) = calculate_asymetric_statistics(
            matching_gene.loc[
                orf_1to1,
                list(set(gene_1to1) & set(depth_hit)),
            ]
        )

        # Calculate F1 scores
        f1 = sp.stats.hmean([precision, recall])
        f1_depth_only = sp.stats.hmean([precision_depth_only, recall_depth_only])
        f1_1to1 = sp.stats.hmean([precision_1to1, recall_1to1])
        f1_depth_only_1to1 = sp.stats.hmean(
            [precision_depth_only_1to1, recall_depth_only_1to1]
        )

        results[strain] = [
            precision,
            recall,
            f1,
            precision_depth_only,
            recall_depth_only,
            f1_depth_only,
            precision_1to1,
            recall_1to1,
            f1_1to1,
            precision_depth_only_1to1,
            recall_depth_only_1to1,
            f1_depth_only_1to1,
        ]

    # Compile results
    results = (
        pd.DataFrame(
            results,
            index=[
                "precision",
                "recall",
                "f1",
                "precision_depth_only",
                "recall_depth_only",
                "f1_depth_only",
                "precision_1to1",
                "recall_1to1",
                "f1_1to1",
                "precision_depth_only_1to1",
                "recall_depth_only_1to1",
                "f1_depth_only_1to1",
            ],
        )
        .T.rename_axis(index="strain")
        .sort_values("f1", ascending=False)
    )

    # Write output
    results.to_csv(outpath, sep="\t")
