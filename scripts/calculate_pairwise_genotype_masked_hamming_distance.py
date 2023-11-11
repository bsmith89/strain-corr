#!/usr/bin/env python3

import sfacts as sf
import sys
from scipy.spatial.distance import pdist, squareform
import numpy as np
from lib.dissimilarity import dump_dmat_as_pickle
import pandas as pd


def native_masked_hamming_distance_pdist(X, pseudo=0):
    isnan_X = np.isnan(X)
    displaced_cityblock = pdist(np.nan_to_num(X, nan=0.5), metric="cityblock")
    mismatched_nan_count = pdist(isnan_X, metric="cityblock")
    cityblock = displaced_cityblock - (mismatched_nan_count / 2)
    anynan_count = (
        np.nan_to_num(pdist(isnan_X, metric="kulczynski1"), nan=0)
        * mismatched_nan_count
        + mismatched_nan_count
    )
    unmasked_positions = X.shape[1] - anynan_count
    masked_hamming_distance = (cityblock + pseudo) / (unmasked_positions + pseudo)
    return masked_hamming_distance


if __name__ == "__main__":
    mgtp_inpath = sys.argv[1]
    ambiguity_threshold = float(
        sys.argv[2]
    )  # Maximum minor allele frequency to replace with NaN
    pseudo = float(sys.argv[3])
    outpath = sys.argv[4]

    est_geno = (
        sf.Metagenotype.load(mgtp_inpath)
        .to_estimated_genotype()
        .discretized(max_ambiguity=ambiguity_threshold)
        .rename_coords(strain=str)
    )

    geno_cdmat = native_masked_hamming_distance_pdist(est_geno.values, pseudo=pseudo)
    labels = est_geno.strain.values
    dmat = pd.DataFrame(squareform(geno_cdmat), index=labels, columns=labels)
    dump_dmat_as_pickle(dmat, path=outpath)
