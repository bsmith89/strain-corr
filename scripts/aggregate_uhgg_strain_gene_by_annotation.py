#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    uhgg_inpath = sys.argv[1]
    agg_inpath = sys.argv[2]
    agg_column = sys.argv[3]
    outpath = sys.argv[4]

    gene_table = pd.read_table(uhgg_inpath, index_col="gene_id")

    # In order to be generic across centroid levels, the agg/annotation
    # data does not specify the name of the index column.
    # It's probably, however, centroid_NN (75 usually) because this
    # is how we aggregated the voting.
    # It's NOT gene_id.
    # On the other hand, the gene_table does use that
    # name for its index, regardless of what centroid_NN was being used by SPGC.
    agg = pd.read_table(agg_inpath, index_col=0).squeeze()
    result = (
        gene_table.astype(bool).join(agg).groupby(agg_column).any().astype(int)
    ).rename_axis(index="gene_id") # Here "gene_id" is the generic name for annotation.
    result.to_csv(outpath, sep="\t")
