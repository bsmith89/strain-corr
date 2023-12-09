#!/usr/bin/env python3
# "{input.script} {input.meta} {input.geno_pdist} {input.gene_pdist} {outpath}"

import sys
from lib.pandas_util import idxwhere
from lib.dissimilarity import load_dmat_as_pickle
import pandas as pd
import numpy as np


if __name__ == "__main__":
    meta_inpath = sys.argv[1]
    geno_pdmat_inpath = sys.argv[2]
    uhgg_pdmat_inpath = sys.argv[3]
    eggnog_pdmat_inpath = sys.argv[4]
    outpath = sys.argv[5]

    # Load data
    meta = pd.read_table(meta_inpath, index_col="genome_id")
    geno_pdmat = load_dmat_as_pickle(geno_pdmat_inpath)
    uhgg_pdmat = load_dmat_as_pickle(uhgg_pdmat_inpath)
    eggnog_pdmat = load_dmat_as_pickle(eggnog_pdmat_inpath)

    ref_list = idxwhere(meta.genome_type.isin(['MAG', 'Isolate']) & meta.passes_filter)

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
        ref_nn_uhgg_diss = (
            uhgg_pdmat.stack()
            .loc[ref_nearest_neighbor.reset_index().apply(tuple, axis=1)]
            .reset_index()
            .set_index("level_0")[0]
        )
        # Not the nearest neighbor, but just the closest gene distance.
        uhgg_pdmat_mask_self = uhgg_pdmat + np.eye(*uhgg_pdmat.shape)
        minimum_ref_uhgg_diss = uhgg_pdmat_mask_self.loc[ref_list].min()

        # Filtered gene content
        ref_nn_eggnog_diss = (
            eggnog_pdmat.stack()
            .loc[ref_nearest_neighbor.reset_index().apply(tuple, axis=1)]
            .reset_index()
            .set_index("level_0")[0]
        )
        # Not the nearest neighbor, but just the closest gene distance.
        eggnog_pdmat_mask_self = eggnog_pdmat + np.eye(*eggnog_pdmat.shape)
        minimum_ref_eggnog_diss = eggnog_pdmat_mask_self.loc[ref_list].min()
    else:
        ref_nearest_neighbor = np.nan
        ref_nn_geno_diss = np.nan

        ref_nn_uhgg_diss = np.nan
        minimum_ref_uhgg_diss = np.nan

        ref_nn_eggnog_diss = np.nan
        minimum_ref_eggnog_diss = np.nan

    # Compile metadata
    data = meta[[]].assign(
        ref_nn_genome_id=ref_nearest_neighbor,
        min_ref_geno_diss=ref_nn_geno_diss,

        ref_nn_uhgg_diss=ref_nn_uhgg_diss,
        min_ref_uhgg_diss=minimum_ref_uhgg_diss,

        ref_nn_eggnog_diss=ref_nn_eggnog_diss,
        min_ref_eggnog_diss=minimum_ref_eggnog_diss,
    )

    # Write output
    data.to_csv(outpath, sep="\t")
