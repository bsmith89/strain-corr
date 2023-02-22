#!/usr/bin/env python3

import pandas as pd
import sys
import scipy as sp
import scipy.stats
from lib.pandas_util import idxwhere
from tqdm import tqdm


if __name__ == "__main__":
    gene_matching_path = sys.argv[1]
    corr_q_path = sys.argv[2]
    depth_q_path = sys.argv[3]
    corr_q_thresh = float(sys.argv[4])
    depth_q_thresh = float(sys.argv[5])
    outpath = sys.argv[6]

    # Load data
    matching_gene = (
        pd.read_table(gene_matching_path, index_col=["orf", "gene"])
        .assign(hit=True)
        .unstack("orf", fill_value=False)
    )
    corr_q_all = pd.read_table(corr_q_path, index_col=["gene_id"])
    depth_q_all = pd.read_table(depth_q_path, index_col=["gene_id"])

    # Align data
    gene_list = list(
        set(corr_q_all.index) | set(depth_q_all.index) | set(matching_gene.index)
    )
    depth_q_all = depth_q_all.reindex(gene_list, fill_value=0)
    corr_q_all = corr_q_all.reindex(gene_list, fill_value=0)
    matching_gene = matching_gene.reindex(gene_list, fill_value=False).T

    # Which genes are 1-to-1 (non-ambiguous mapping)?
    gene_multi_hit = matching_gene.sum(0) > 1
    gene_1to1 = ~gene_multi_hit
    orf_1to1 = (matching_gene.sum(1) == 1) & ~matching_gene.loc[:, gene_multi_hit].any(
        1
    )
    gene_1to1 = set(idxwhere(gene_1to1))
    orf_1to1 = set(idxwhere(orf_1to1))

    assert set(corr_q_all.columns) == set(depth_q_all.columns)
    precision_results = []
    recall_results = []
    f1_results = []
    precision_1to1_results = []
    recall_1to1_results = []
    f1_1to1_results = []
    for strain in tqdm(corr_q_all.columns):
        corr_q = corr_q_all[strain]
        depth_q = depth_q_all[strain]

        # Classify genes and orfs
        depth_and_corr_hit = (corr_q >= corr_q_thresh) & (depth_q >= depth_q_thresh)
        tp_gene = set(
            idxwhere(matching_gene.loc[:, depth_and_corr_hit].any(axis=0))
        )  # Genes hit by both.
        fp_gene = set(
            idxwhere(~(matching_gene.loc[:, depth_and_corr_hit].any()))
        )  # Genes were hit by SPGC but never by BLAST?
        fn_orf = set(
            idxwhere(~(matching_gene.loc[:, depth_and_corr_hit].any(axis=1)))
        )  # How many ORFs were hit by BLAST but no matching SPGC hits?
        tp_orf = set(
            idxwhere((depth_and_corr_hit & matching_gene).any(axis=1))
        )  # How many ORFs were hit by BLAST and by SPGC?

        # Calculate accuracy
        n_tp_gene = len(tp_gene)
        n_fp_gene = len(fp_gene)
        n_fn_orf = len(fn_orf)
        n_tp_orf = len(tp_orf)

        n_tp_1to1 = len(tp_gene & gene_1to1)
        n_fp_1to1 = len(fp_gene & gene_1to1)
        n_fn_1to1 = len(fn_orf & orf_1to1)
        n_tn_1to1 = len(set(depth_and_corr_hit) & gene_1to1) - (
            n_tp_1to1 + n_fp_1to1 + n_fn_1to1
        )

        precision = n_tp_gene / (n_tp_gene + n_fp_gene)
        recall = n_tp_orf / (n_tp_orf + n_fn_orf)
        precision_1to1 = n_tp_1to1 / (n_tp_1to1 + n_fp_1to1)
        recall_1to1 = n_tp_1to1 / (n_tp_1to1 + n_fn_1to1)
        f1 = sp.stats.hmean([precision, recall])
        f1_1to1 = sp.stats.hmean([precision_1to1, recall_1to1])

        precision_results.append(precision)
        recall_results.append(recall)
        f1_results.append(f1)
        precision_1to1_results.append(precision_1to1)
        recall_1to1_results.append(recall_1to1)
        f1_1to1_results.append(f1_1to1)

    # Compile results
    results = pd.DataFrame(
        dict(
            precision=precision_results,
            recall=recall_results,
            f1=f1_results,
            precision_1to1=precision_1to1_results,
            recall_1to1=recall_1to1_results,
            f1_1to1=f1_1to1_results,
        ),
        index=corr_q_all.columns,
    ).sort_values('f1_1to1', ascending=False)

    # Write output
    results.to_csv(outpath, sep='\t')
