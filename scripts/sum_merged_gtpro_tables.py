#!/usr/bin/env python3

import pandas as pd
import sys
from tqdm import tqdm


if __name__ == "__main__":
    *inpaths, outpath = sys.argv[1:]
    mgen_list = []
    for path in tqdm(inpaths, bar_format="{l_bar}{r_bar}"):
        mgen_list.append(
            pd.read_table(
                path,
                names=[
                    "sample_id",
                    "species_id",
                    "global_pos",
                    "contig",
                    "local_pos",
                    "ref_allele",
                    "alt_allele",
                    "ref_count",
                    "alt_count",
                ],
            )
        )
    mgen = (
        pd.concat(mgen_list)
        .groupby(
            [
                "sample_id",
                "species_id",
                "global_pos",
                "contig",
                "local_pos",
                "ref_allele",
                "alt_allele",
            ]
        )
        .sum()
    )
    mgen.to_csv(outpath, sep="\t", header=False)
