#!/usr/bin/env python3

import xarray as xr
import sfacts as sf
import sys
from scipy.spatial.distance import pdist
import numpy as np
import pickle
from lib.util import info


if __name__ == "__main__":
    ref_geno_inpath = sys.argv[1]
    ambiguity_threshold = float(sys.argv[2])
    outpath = sys.argv[3]

    ref_geno = (
        sf.Metagenotype.load(ref_geno_inpath)
        .to_estimated_genotype()
        .discretized(max_ambiguity=ambiguity_threshold)
        # .rename_coords(strain=lambda x: "UHGG" + x[len("GUT_GENOME") :])
    )

    # Count num positions.
    ref_positions = (~ref_geno.data.to_pandas().isna()).sum(1)

    # Write output
    ref_positions.to_frame("num_geno_positions").rename_axis("genome_id").to_csv(
        outpath, sep="\t"
    )
