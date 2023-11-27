#!/usr/bin/env python3

import xarray as xr
import sfacts as sf
import sys
from scipy.spatial.distance import pdist, squareform
import numpy as np
from lib.util import info
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
    spgc_agg_mgtp_inpath = sys.argv[1]
    ref_geno_inpath = sys.argv[2]
    ambiguity_threshold = float(
        sys.argv[3]
    )  # Maximum minor allele frequency to replace with NaN
    pseudo = float(sys.argv[4])
    outpath = sys.argv[5]

    spgc_est_geno = (
        sf.Metagenotype.load(spgc_agg_mgtp_inpath)
        .to_estimated_genotype()
        .discretized(max_ambiguity=ambiguity_threshold)
        .rename_coords(strain=str)
    )
    ref_geno = (
        sf.Metagenotype.load(ref_geno_inpath)
        .to_estimated_genotype()
        .discretized(max_ambiguity=ambiguity_threshold)
        # .rename_coords(strain=lambda x: "UHGG" + x[len("GUT_GENOME") :])
    )
    # ref_geno.data["strain"] = (
    #     ref_geno.data["strain"]
    #     .to_series()
    #     .map(lambda x: "UHGG" + x[len("GUT_GENOME") :])
    # )

    # Debug num positions in each dataset.
    spgc_positions = set(spgc_est_geno.position.values)
    ref_positions = set(ref_geno.position.values)
    info(
        "Postions:",
        "spgc-private",
        "shared",
        "ref-private",
    )
    info(
        len(spgc_positions - ref_positions),
        len(spgc_positions & ref_positions),
        len(ref_positions - spgc_positions),
    )
    assert len(spgc_positions - ref_positions) == 0

    geno = sf.Genotype.concat(
        dict(
            ref=ref_geno.sel(position=spgc_est_geno.position),
            spgc=spgc_est_geno,
        ),
        dim="strain",
        rename=False,
    )
    geno_cdmat = native_masked_hamming_distance_pdist(geno.values, pseudo=pseudo)
    labels = geno.strain.values
    dmat = pd.DataFrame(squareform(geno_cdmat), index=labels, columns=labels)
    dump_dmat_as_pickle(dmat, path=outpath)
