#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    annot_inpath = sys.argv[1]
    outpath = sys.argv[2]

    # NOTE: The index column is labeled "centroid_99" but it is in reality "gene_id".
    annot_table = (
        pd.read_table(annot_inpath)
        .rename(columns={"centroid_99": "gene_id"})[["gene_id", "accession_no"]]
        .drop_duplicates()
    )
    annot_table.to_csv(outpath, sep="\t", index=False)
