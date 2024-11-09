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


if __name__ == "__main__":
    eggnog = (
        pd.read_table(sys.argv[1], comment="#", names=EGGNOG_COLUMNS, index_col="query")
        .rename_axis(index="centroid_99")
        .replace({"-": np.nan})
        .assign(eggNOG_OGs=lambda x: x.eggNOG_OGs.str.split(","))
        .explode("eggNOG_OGs")
        .assign(annot_lvl=lambda x: x.eggNOG_OGs.str.split("@").str[1])[
            lambda x: x.annot_lvl == x.max_annot_lvl
        ]
    )
    gene_x_eggnog = eggnog.eggNOG_OGs.dropna().rename("top_eggnog")
    gene_x_eggnog.to_csv(sys.argv[2], sep="\t")
