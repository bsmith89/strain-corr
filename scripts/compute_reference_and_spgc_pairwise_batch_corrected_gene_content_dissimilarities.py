#!/usr/bin/env python3

import pandas as pd
import sys
from lib.pandas_util import idxwhere, align_indexes
from lib.dissimilarity import dump_dmat_as_pickle, dmatrix


if __name__ == "__main__":
    spgc_gene_inpath = sys.argv[1]
    ref_gene_inpath = sys.argv[2]
    spgc_filt_inpath = sys.argv[3]
    ref_filt_inpath = sys.argv[4]
    metric = sys.argv[5]
    outpath = sys.argv[6]

    # Load and align gene content and metadata for SPGC and Ref
    spgc_gene = pd.read_table(spgc_gene_inpath, index_col="gene_id").T
    spgc_filt = pd.read_table(spgc_filt_inpath, index_col="genome_id").rename(str)
    spgc_gene, spgc_filt = align_indexes(spgc_gene, spgc_filt)

    ref_gene = pd.read_table(ref_gene_inpath, index_col="gene_id").T
    ref_filt = pd.read_table(ref_filt_inpath, index_col="genome_id")
    ref_gene, ref_filt = align_indexes(ref_gene, ref_filt)

    # Align genes across the columns
    spgc_gene, ref_gene = align_indexes(
        spgc_gene, ref_gene, axis="columns", how="outer", fill_value=0
    )

    spgc_strain_filt_list = idxwhere(spgc_filt.passes_filter)
    ref_strain_filt_list = idxwhere(ref_filt.passes_filter)

    correction = spgc_gene.loc[spgc_strain_filt_list].mean() - ref_gene.loc[ref_strain_filt_list].mean()

    data = pd.concat(
        [spgc_gene - correction, ref_gene],
    )

    pdmat = dmatrix(data, metric=metric)
    dump_dmat_as_pickle(pdmat, outpath)
