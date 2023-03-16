#!/usr/bin/env python3

import pandas as pd
import sys
from lib.pandas_util import idxwhere, align_indexes


if __name__ == "__main__":
    species_depth_path = sys.argv[1]
    strain_frac_path = sys.argv[2]
    frac_thresh = float(sys.argv[3])
    absent_thresh = float(sys.argv[4])
    present_thresh = float(sys.argv[5])
    output_nospecies_path = sys.argv[6]
    output_strain_path = sys.argv[7]

    assert (
        frac_thresh >= 0.5
    ), "frac_thresh < 0.5 might prevent complete partitioning of samples."

    species_depth = pd.read_table(
        species_depth_path,
        names=["sample", "depth"],
        index_col="sample",
    ).squeeze()
    strain_frac = (
        pd.read_table(strain_frac_path, index_col=["sample", "strain"])
        .squeeze()
        .unstack("strain")
    )
    with open(output_nospecies_path, "w") as f:
        for sample in idxwhere(species_depth < absent_thresh):
            print(sample, file=f)

    species_depth, strain_frac = align_indexes(
        species_depth, strain_frac, axis="index", how="inner"
    )
    strain_samples = (strain_frac > frac_thresh).apply(
        lambda x: x & (species_depth > present_thresh)
    )

    strain_sample_output = strain_samples.T.stack()[lambda x: x].index
    strain_sample_output.to_frame().to_csv(output_strain_path, sep="\t", index=False)
