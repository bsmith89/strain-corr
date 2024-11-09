#!/usr/bin/env python3

# USAGE: {input.script} {input.prevalence} {params.threshold} {output}

import pandas as pd
import sys
from lib.util import info
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    prevalence_inpath = sys.argv[1]
    threshold = float(sys.argv[2])
    outpath = sys.argv[3]

    info("Loading reference table.")
    prevalence = pd.read_table(
        prevalence_inpath, names=["gene_id", "prevalence"], index_col="gene_id"
    ).prevalence

    high_prevalence_genes = idxwhere((prevalence > threshold))
    # Write output.
    with open(outpath, "w") as f:
        for gene in high_prevalence_genes:
            print(gene, file=f)
