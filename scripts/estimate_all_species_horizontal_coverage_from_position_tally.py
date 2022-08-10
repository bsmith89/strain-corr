#!/usr/bin/env python3

import sys
import pandas as pd


def info(*msg):
    print(*msg, file=sys.stderr)


if __name__ == "__main__":
    snp_dict_path = sys.argv[1]
    tally_path = sys.argv[2]
    outpath = sys.argv[3]

    info("Loading reference table.")
    ref_species_snp_count = (
        pd.read_table(
            snp_dict_path,
            names=[
                "species_id",
                "species_position_id",
                "_2",
                "_3",
                "_4",
                "_5",
            ],
        )
        .groupby(["species_id"])
        .species_position_id.count()
        .rename("tally")
    )

    info("Loading data.")
    species_snp_count = (
        pd.read_table(
            tally_path,
            names=["sample_id", "species_id", "tally"],
            index_col=["sample_id", "species_id"],
        )
        .squeeze()
        .unstack("species_id")
        .fillna(0)
        .astype(int)
    )

    info("Calculating horizontal coverage.")
    horizontal_coverage = (
        (species_snp_count / ref_species_snp_count)
        .stack()
        .rename("horizontal_coverage")
    )

    info("Writing output.")
    horizontal_coverage.to_csv(outpath, header=False, sep="\t")
