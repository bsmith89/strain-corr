#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    species_depth_path = sys.argv[1]
    strain_frac_path = sys.argv[2]
    frac_thresh = float(sys.argv[3])
    absent_thresh = float(sys.argv[4])
    present_thresh = float(sys.argv[5])
    output_nospecies_path = sys.argv[6]
    output_strain_path = sys.argv[7]

    species_depth = pd.read_table(
        species_depth_path, index_col=["sample", "species_id"]
    )

    assert len(species_depth.species_id.unique()) == 1
    species_depth = species_depth.squeeze().unstack("species_id").squeeze()
    strain_frac = (
        pd.read_table(strain_frac_path, index_col=["sample", "strain"])
        .squeeze()
        .unstack("strain")
    )

    # species_depth <
