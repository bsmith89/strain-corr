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
    )
    gene_x_ko = eggnog.KEGG_ko.dropna().str.split(",").explode().rename("ko")
    gene_x_ko.to_csv(sys.argv[2], sep="\t")
