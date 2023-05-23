#!/usr/bin/env python3
# {input.script} {input.truth} {input.infer} {output}

import pandas as pd
import sys
import scipy as sp
import scipy.stats
from lib.pandas_util import idxwhere
import numpy as np
from tqdm import tqdm


def accuracy_stats(truth, predicted):
    truth = set(truth)
    predicted = set(predicted)
    false_negatives = truth - predicted
    false_positives = predicted - truth
    true_positives = truth & predicted

    n_t = len(truth)
    n_p = len(predicted)
    n_fn = len(false_negatives)
    n_fp = len(false_positives)
    n_tp = len(true_positives)

    if n_p > 0:
        precision = n_tp / n_p
    else:
        precision = np.nan

    if n_t > 0:
        recall = n_tp / n_t
    else:
        recall = np.nan

    if np.isnan([precision, recall]).any():
        f1 = np.nan
    else:
        f1 = sp.stats.hmean([precision, recall])

    return precision, recall, f1


if __name__ == "__main__":
    truth_path = sys.argv[1]
    infer_path = sys.argv[2]
    outpath = sys.argv[3]

    truth_genes = pd.read_table(truth_path, index_col=0)
    assert truth_genes.shape[1] == 1
    truth_genes = truth_genes.iloc[:, 0]
    infer_genes = pd.read_table(infer_path, index_col=0).rename_axis(
        columns="strain"
    )
    assert truth_genes.index.name == infer_genes.index.name
    truth = set(idxwhere(truth_genes.astype(bool)))
    out = dict()
    for strain in tqdm(infer_genes.columns):
        predicted = set(idxwhere(infer_genes[strain] == 1))
        precision, recall, f1 = accuracy_stats(truth, predicted)
        out[strain] = dict(
            precision=precision,
            recall=recall,
            f1=f1,
        )

    out = pd.DataFrame(out).T.rename_axis(index="strain")
    out.to_csv(outpath, sep="\t")
