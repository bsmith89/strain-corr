#!/usr/bin/env python3

import pandas as pd
import sfacts as sf
import sys
import scipy as sp
from lib.dissimilarity import load_dmat_as_pickle


if __name__ == "__main__":
    dmat_inpath = sys.argv[1]
    spgc_agg_mgtp_inpath = sys.argv[2]
    outpath = sys.argv[3]

    dmat = load_dmat_as_pickle(dmat_inpath)

    spgc_mgtp = sf.Metagenotype.load(spgc_agg_mgtp_inpath)
    spgc_list = list(spgc_mgtp.sample.astype(str).values)
    ref_list = list(set(dmat.index) - set(spgc_list))

    result = dmat.loc[ref_list, spgc_list].stack().rename_axis(['genome_id', 'strain']).rename("genotype_dissimilarity")
    result.to_csv(outpath, sep="\t")
