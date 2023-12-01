#!/usr/bin/env python
# "{input.script} {input.pickle} {input.spgc_filt} {input.ref_filt} {params.thresh} {output}"

import pandas as pd
import sys
from scipy.spatial.distance import pdist, squareform
import pickle
from lib.pandas_util import idxwhere
import numpy as np
from lib.dissimilarity import dump_dmat_as_pickle


def prevalence(x, pseudo=0):
    "Calculate the prevalence in a binary vector, potentially with a pseudocounts added."
    pseudo = 1
    return (x.sum() + pseudo) / (len(x) + 2 * pseudo)


if __name__ == "__main__":
    spgc_gene_inpath = sys.argv[1]
    ref_gene_inpath = sys.argv[2]
    spgc_filt_inpath = sys.argv[3]
    ref_filt_inpath = sys.argv[4]
    maximum_prevalence_ratio = float(sys.argv[5])
    metric = sys.argv[6]
    outpath = sys.argv[7]

    spgc_gene = pd.read_table(spgc_gene_inpath, index_col="gene_id")
    ref_gene = pd.read_table(ref_gene_inpath, index_col="gene_id")
    spgc_filt = pd.read_table(spgc_filt_inpath, index_col="genome_id").rename(str)
    ref_filt = pd.read_table(ref_filt_inpath, index_col="genome_id")

    gene_list = list(set(ref_gene.index) | set(spgc_gene.index))

    spgc_list = idxwhere(spgc_filt.passes_filter)
    ref_list = idxwhere(ref_filt.passes_filter)

    prevalence_pc = pd.DataFrame(
        dict(
            ref=ref_gene.reindex(index=gene_list, columns=ref_list, fill_value=0).apply(
                prevalence, pseudo=1, axis=1
            ),
            spgc=spgc_gene.reindex(index=gene_list, columns=spgc_list, fill_value=0).apply(
                prevalence, pseudo=1, axis=1
            ),
        )
    ).assign(log_ratio=lambda x: np.log(x.spgc) - np.log(x.ref))
    gene_list = idxwhere(np.abs(prevalence_pc.log_ratio) < np.abs(np.log(maximum_prevalence_ratio)))

    assert len(gene_list) > 0

    data = (
        pd.concat(
            [spgc_gene.reindex(gene_list), ref_gene.reindex(gene_list)],
            axis=1,
        )
        .fillna(0)
        .T.astype(bool)
    )

    cdmat = pdist(data, metric=metric)
    pdmat = pd.DataFrame(squareform(cdmat), index=data.index, columns=data.index)
    dump_dmat_as_pickle(pdmat, outpath)
