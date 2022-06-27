#!/usr/bin/env python3

import pandas as pd
import sys
from tqdm import tqdm
import numpy as np


DUMMY_STRAIN_LABEL = '-1'
MINIMUM_FRAC = 1e-3

if __name__ == "__main__":
    species_path = sys.argv[1]
    outpath = sys.argv[2]
    species_path_map = {
        k: v for k, v in [arg.strip().split("=") for arg in sys.argv[3:]]
    }

    # Load species level coverage.
    species_depth = (
        pd.read_table(
            species_path,
            dtype=dict(species_id=str, sample=str, depth=float),
            index_col=["sample", "species_id"],
        )
        .squeeze()
        .unstack("species_id", fill_value=0)
    )
    species_depth.columns = species_depth.columns.astype(str)
    assert species_depth.index.is_unique
    assert species_depth.columns.is_unique

    sample_list = species_depth.index
    species_list = species_depth.columns

    # Compile intra-species relative abundance multiplied by the species
    # total coverage.
    strain_depth = {}
    for species in tqdm(species_list):
        if species in species_path_map:
            strain_frac = (
                pd.read_table(
                    species_path_map[species],
                    dtype=dict(strain=str, sample=str, community=float),
                    index_col=["sample", "strain"],
                )
                .squeeze().unstack("strain", fill_value=0)
            )
            assert np.allclose(strain_frac.sum(1), 1.0, atol=1e-4)
            # Drop existing dummy strain fraction (we'll add it back later).
            strain_frac = (
                strain_frac
                .reindex(sample_list, fill_value=0)
                .drop(columns=[DUMMY_STRAIN_LABEL], errors='ignore')
            )
            # Drop any strain fractions less than MINIMUM_FRAC.
            strain_frac = strain_frac.where(strain_frac > MINIMUM_FRAC, 0)
            # Add back the dummy strain fraction as everything needed to add
            # up to 1.0.
            strain_frac[DUMMY_STRAIN_LABEL] = 1 - strain_frac.sum("columns")
            assert np.allclose(strain_frac.sum(1), 1.0, atol=1e-4)
        else:
            # When no strains were provided, mark everything as the dummy strain.
            strain_frac = pd.DataFrame({DUMMY_STRAIN_LABEL: 1}, index=sample_list)
        # Rename strains by prepending the species_id.
        strain_frac = strain_frac.rename(columns=lambda s: species + '-' + s)
        # Multiply by the species_depth to get the strain depth.
        strain_depth[species] = (strain_frac.T * species_depth[species]).T

    strain_depth = pd.concat(strain_depth.values(), axis=1)
    assert np.allclose(strain_depth.sum("columns"), species_depth.sum("columns"), atol=1e-4)
    # Write stacked table (dropping 0s to keep the file-size small).
    strain_depth.where(strain_depth > 0, np.nan).stack().dropna().rename("depth").to_csv(outpath, sep='\t')
