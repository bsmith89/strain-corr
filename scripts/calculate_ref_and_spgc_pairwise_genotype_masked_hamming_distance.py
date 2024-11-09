#!/usr/bin/env python3

import sfacts as sf
import sys
from scipy.spatial.distance import squareform
from lib.util import info
from lib.dissimilarity import dump_dmat_as_pickle
from lib.thisproject.genotype_dissimilarity import native_masked_hamming_distance_pdist
import pandas as pd
from warnings import warn


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
    if len(spgc_positions - ref_positions) == 0:
        warn(
            "Are you sure you passed a reference genotype? Not all positions are found in ref."
        )

    geno = sf.Genotype.concat(
        dict(
            ref=ref_geno,  # NOTE: We could trim this geno down to just existing positions: .sel(position=spgc_est_geno.position),
            spgc=spgc_est_geno,
        ),
        dim="strain",
        rename=False,
    )
    geno_cdmat = native_masked_hamming_distance_pdist(geno.values, pseudo=pseudo)
    labels = geno.strain.values
    dmat = pd.DataFrame(squareform(geno_cdmat), index=labels, columns=labels)
    dump_dmat_as_pickle(dmat, path=outpath)
