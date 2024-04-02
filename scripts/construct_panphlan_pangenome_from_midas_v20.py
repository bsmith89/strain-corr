#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    gene_info_path = sys.argv[1]
    centroid = sys.argv[2]
    outpath = sys.argv[3]

    centroid_col = f"centroid_{centroid}"

    gene_info = (pd.read_table(gene_info_path).assign(
        genome_id=lambda x: x.gene_id.str.rsplit("_", n=1).str[0]
        )
        .assign(
            dummy_contig="dummy_contig",  # Probably doesn't matter, right?
            dummy_left=0,  # And right will be the length of the gene?
        )
    )

    gene_info[
        [
            centroid_col,
            "gene_id",
            "genome_id",
            "dummy_contig",
            "dummy_left",
            "gene_length",
        ]
    ].to_csv(outpath, sep="\t", header=False, index=False)
