#!/usr/bin/env python3
# "{input.script} {input.spgc_gene} {input.spgc_filt} {input.ref_gene} {input.ref_filt} {input.pdist} {output}"

import pandas as pd
from lib.dissimilarity import load_dmat_as_pickle
import numpy as np
from lib.pandas_util import idxwhere
import sys


def morans_i(x, w, return_zscore=False, centered=False):
    """TODO

    w should be a weight matrix, e.g.

        w = np.exp(-(strain_geno_pdist.loc[_genome_list, _genome_list].values))
        w[np.diag_indices_from(w)] = 0

    """
    # Moran's I.
    x_centered = x - x.mean()
    denomenator = (x_centered ** 2).sum()
    if denomenator == 0:
        return np.nan

    n = len(x)
    w_sum = w.sum()

    numerator = np.einsum("i,j,ij->", x_centered, x_centered, w)
    normalizer = n / w_sum
    observed = numerator * normalizer / denomenator

    # Expected value of Moran's I under the null hypothesis:
    expected = -1 / (len(x) - 1)
    if centered:
        observed_out = observed - expected
    else:
        observed_out = observed

    if return_zscore:
        # Expected variance under the null
        s1 = 0.5 * ((w + w.T) ** 2).sum()
        s2 = ((w.sum(0) + w.sum(1)) ** 2).sum()
        s3 = ((1 / n) * (x_centered ** 4).sum()) / (
            (1 / n) * x_centered ** 2
        ).sum() ** 2
        s4 = (n ** 2 - 3 * n + 3) * s1 - n * s2 + 3 * w_sum ** 2
        s5 = (n ** 2 - n) * s1 - 2 * n * s2 + 6 * w_sum ** 2

        i_var = (n * s4 - s3 * s5) / (
            (n - 1) * (n - 2) * (n - 3) * w_sum ** 2
        ) - expected ** 2

        zscore = (observed - expected) / np.sqrt(i_var)
        return observed_out, zscore
    else:
        return observed_out


def geno_pdist_to_weights(pdist, rate=2):
    pdist = np.array(pdist)
    w = np.exp(-rate * pdist)
    w[np.diag_indices_from(w)] = 0
    return w


if __name__ == "__main__":
    spgc_gene_inpath = sys.argv[1]
    spgc_filt_inpath = sys.argv[2]
    ref_gene_inpath = sys.argv[3]
    ref_filt_inpath = sys.argv[4]
    pdist_inpath = sys.argv[5]
    outpath = sys.argv[6]

    dmat = load_dmat_as_pickle(pdist_inpath)

    spgc_meta = pd.read_table(spgc_filt_inpath, index_col="genome_id").rename(str)
    spgc_gene = pd.read_table(spgc_gene_inpath, index_col="gene_id").rename_axis(
        columns="genome_id"
    )
    ref_meta = pd.read_table(ref_filt_inpath, index_col="genome_id")
    ref_gene = pd.read_table(ref_gene_inpath, index_col="gene_id").rename_axis(
        columns="genome_id"
    )

    # SPGC
    _genome_list = idxwhere(spgc_meta.passes_filter)
    _meta = spgc_meta.loc[_genome_list]
    _gene = spgc_gene.loc[:, _genome_list]
    _dmat = dmat.loc[_genome_list, _genome_list]
    _weights = geno_pdist_to_weights(_dmat, rate=1)
    results = {}
    for gene_id in _gene.index:
        results[gene_id] = morans_i(
            _gene.loc[gene_id].values,
            _weights,
            return_zscore=False,
            centered=True,
        )
    spgc_results = pd.Series(results)

    # Ref
    _genome_list = idxwhere(ref_meta.passes_filter)
    _meta = ref_meta.loc[_genome_list]
    _gene = ref_gene.loc[:, _genome_list]
    _dmat = dmat.loc[_genome_list, _genome_list]
    _weights = geno_pdist_to_weights(_dmat, rate=1)
    results = {}
    for gene_id in _gene.index:
        results[gene_id] = morans_i(
            _gene.loc[gene_id].values,
            _weights,
            return_zscore=False,
            centered=True,
        )
    ref_results = pd.Series(results)

    results = pd.DataFrame(dict(spgc=spgc_results, ref=ref_results)).rename_axis(index="gene_id")
    results.to_csv(outpath, sep="\t")
