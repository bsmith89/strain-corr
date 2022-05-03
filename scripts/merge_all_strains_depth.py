#!/usr/bin/env python3

import pandas as pd
import sys
from tqdm import tqdm


if __name__ == "__main__":
    species_path = sys.argv[1]
    outpath = sys.argv[2]
    species_strains_path_map = {
        k: v for k, v in [arg.strip().split("=") for arg in sys.argv[3:]]
    }

    # Load species level coverage.
    motu_cv = pd.read_table(
        species_path, index_col=["sample", "species_id"], squeeze=True,
    ).unstack("species_id", fill_value=0)
    motu_cv.columns = motu_cv.columns.astype(str)
    assert motu_cv.index.is_unique
    assert motu_cv.columns.is_unique

    sample_list = motu_cv.index

    # Compile intra-species relative abundance multiplied by the species
    # total coverage.
    sotu_cv = {}
    for motu in tqdm(species_strains_path_map):
        path = species_strains_path_map[motu]
        frac = (
            pd.read_table(
                path, index_col=["sample", "strain"], squeeze=True,
            )
            .unstack("strain", fill_value=0)
            .reindex(sample_list, fill_value=0)
        )
        if motu not in motu_cv.columns:
            motu_cv[motu] = 0
        sotu_cv[motu] = (
            frac.assign(other=1 - frac.sum(1))
            .rename(columns=lambda strain_i: motu + "-" + str(strain_i))
            .rename_axis(columns="sotu")
            .multiply(motu_cv[motu], axis=0)
            .stack()
        )
    # Add entries for all species that didn't have strain splits.
    motu_cv_missing = motu_cv.rename(columns=lambda s: s + "-other").drop(
        columns=[k + "-other" for k in sotu_cv.keys()]
    )
    sotu_cv = (
        pd.concat(sotu_cv.values())
        .unstack("sotu", fill_value=0)
        .join(motu_cv_missing)
    )
    # FIXME: Some of the '-other' strains are < 0 because of rounding error.
    sotu_cv[(sotu_cv < 0)] = 0

    assert sotu_cv.index.is_unique
    assert sotu_cv.columns.is_unique
    assert (sotu_cv >= 0).all().all()

    sotu_cv.stack()[lambda x: x > 0].to_csv(outpath, sep="\t", header=True)
