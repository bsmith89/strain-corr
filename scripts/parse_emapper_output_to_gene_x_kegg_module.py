#!/usr/bin/env python3

import pandas as pd
import sys
import numpy as np

EGGNOG_COLUMNS = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
    "contig_id",
    "start",
    "end",
    "strand",
    "gene_type",
    "contig_length",
]

INPUT_ANNOT_NAME = "KEGG_Module"
OUTPUT_ANNOT_NAME = "kegg_module"

if __name__ == "__main__":
    # Load emapper table
    emapper = (
        pd.read_table(
            sys.argv[1],
            comment="#",
            names=EGGNOG_COLUMNS,
            index_col="query",
        )
        .rename_axis(index="centroid_99")
        .replace({"-": np.nan})
    )

    # Parse multiple entries in annot
    gene_x_annot = (
        emapper.assign(annot_list=lambda x: x[INPUT_ANNOT_NAME].str.split(","))[
            ["annot_list"]
        ]
        .explode("annot_list")
        .dropna()
    )

    # Output
    gene_x_annot['annot_list'].rename(OUTPUT_ANNOT_NAME).to_csv(sys.argv[2], sep="\t")
