#!/usr/bin/env python
# "{input.script} {input.pickle} {input.spgc_filt} {input.ref_filt} {params.thresh} {output}"

import pandas as pd
import sys
from scipy.spatial.distance import pdist
import pickle
from lib.pandas_util import idxwhere
import numpy as np

if __name__ == "__main__":
    spgc_gene_inpath = sys.argv[1]
    ref_gene_inpath = sys.argv[2]
    spgc_prev_inpath = sys.argv[3]
    ref_prev_inpath = sys.argv[4]
    maximum_prevalence_ratio = float(sys.argv[5])
    outpath = sys.argv[6]

    spgc_prev = pd.read_table(spgc_prev_inpath, names=['gene_id', 'prevalence'], index_col="gene_id").prevalence
    ref_prev = pd.read_table(ref_prev_inpath, names=['gene_id', 'prevalence'], index_col="gene_id").prevalence

    prev = (
        pd.DataFrame(dict(spgc=spgc_prev, ref=ref_prev))
        .fillna(0)
        .assign(
            log_spgc_to_ref_ratio=lambda x: np.log(x.spgc / x.ref),
            passes_filter=lambda x: np.abs(x.log_spgc_to_ref_ratio)
            < np.abs(np.log(maximum_prevalence_ratio)),
        )
    )

    gene_list = idxwhere(prev.passes_filter)
    data = pd.concat(
        [
            pd.read_table(spgc_gene_inpath, index_col="gene_id").reindex(gene_list),
            pd.read_table(ref_gene_inpath, index_col="gene_id").reindex(gene_list),
        ],
        axis=1,
    ).fillna(0).T.astype(bool)

    cdmat = pdist(data, metric="jaccard")
    with open(outpath, "wb") as f:
        pickle.dump(dict(cdmat=cdmat, labels=data.index.values), file=f)
