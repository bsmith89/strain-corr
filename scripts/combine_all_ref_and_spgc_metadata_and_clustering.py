#!/usr/bin/env python3
# "{input.script} {input.spgc} {input.ref} {input.geno_pdist} {input.gene_pdist} {input.clust} {outpath}"

import sys
from lib.pandas_util import idxwhere
from lib.dissimilarity import load_dmat_as_pickle
import pandas as pd
import numpy as np


if __name__ == "__main__":
    spgc_meta_inpath = sys.argv[1]
    ref_meta_inpath = sys.argv[2]
    geno_pdmat_inpath = sys.argv[3]
    gene_raw_pdmat_inpath = sys.argv[4]
    gene_filt_pdmat_inpath = sys.argv[5]
    clust_inpath = sys.argv[6]
    outpath = sys.argv[7]

    # Load data
    spgc_meta = pd.read_table(spgc_meta_inpath, index_col="genome_id").rename(str)
    ref_meta = pd.read_table(ref_meta_inpath, index_col="genome_id").rename(
        columns={"passes_num_positions": "passes_geno_positions"}
    )  # NOTE (2023-09-19): This is a stopgap and should be fixed upstream.
    geno_pdmat = load_dmat_as_pickle(geno_pdmat_inpath)
    gene_raw_pdmat = load_dmat_as_pickle(gene_raw_pdmat_inpath)
    gene_filt_pdmat = load_dmat_as_pickle(gene_filt_pdmat_inpath)
    clust = pd.read_table(
        clust_inpath, names=["genome_id", "clust"], index_col="genome_id"
    ).clust

    ref_list = idxwhere(ref_meta.passes_filter)

    # Find closest reference genotype distance for each genome (excluding self).
    geno_pdmat_mask_self = geno_pdmat + np.eye(*geno_pdmat.shape)

    if len(ref_list) >= 2:  # Otherwise the sole reference doesn't have anything to be compared to.
        ref_nearest_neighbor = geno_pdmat_mask_self.loc[ref_list].idxmin()
        # Query the closest genotype and gene content distance.
        ref_nn_geno_diss = (
            geno_pdmat.stack()
            .loc[ref_nearest_neighbor.reset_index().apply(tuple, axis=1)]
            .reset_index()
            .set_index("level_0")[0]
        )

        # Raw gene content
        ref_nn_gene_raw_diss = (
            gene_raw_pdmat.stack()
            .loc[ref_nearest_neighbor.reset_index().apply(tuple, axis=1)]
            .reset_index()
            .set_index("level_0")[0]
        )
        # Not the nearest neighbor, but just the closest gene distance.
        gene_raw_pdmat_mask_self = gene_raw_pdmat + np.eye(*gene_raw_pdmat.shape)
        minimum_ref_gene_raw_diss = gene_raw_pdmat_mask_self.loc[ref_list].min()

        # Filtered gene content
        ref_nn_gene_filt_diss = (
            gene_filt_pdmat.stack()
            .loc[ref_nearest_neighbor.reset_index().apply(tuple, axis=1)]
            .reset_index()
            .set_index("level_0")[0]
        )
        # Not the nearest neighbor, but just the closest gene distance.
        gene_filt_pdmat_mask_self = gene_filt_pdmat + np.eye(*gene_filt_pdmat.shape)
        minimum_ref_gene_filt_diss = gene_filt_pdmat_mask_self.loc[ref_list].min()
    else:
        ref_nearest_neighbor = np.nan
        ref_nn_geno_diss = np.nan

        ref_nn_gene_raw_diss = np.nan
        minimum_ref_gene_raw_diss = np.nan

        ref_nn_gene_filt_diss = np.nan
        minimum_ref_gene_filt_diss = np.nan

    # Compile metadata
    data = pd.concat([spgc_meta, ref_meta]).assign(
        ref_nn_genome_id=ref_nearest_neighbor,
        min_ref_geno_diss=ref_nn_geno_diss,

        ref_nn_gene_raw_diss=ref_nn_gene_raw_diss,
        min_ref_gene_raw_diss=minimum_ref_gene_raw_diss,

        ref_nn_gene_filt_diss=ref_nn_gene_filt_diss,
        min_ref_gene_filt_diss=minimum_ref_gene_filt_diss,

        clust=clust,
    )

    # Write output
    data.to_csv(outpath, sep="\t")
