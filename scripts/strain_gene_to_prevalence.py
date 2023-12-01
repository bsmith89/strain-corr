#!/usr/bin/env python3
# "{input.script} {input.gene} {output}"

import sys
import pandas as pd
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    gene_presence_inpath = sys.argv[1]
    genome_filt_inpath = sys.argv[2]
    pseudo = float(sys.argv[3])
    outpath = sys.argv[4]

    genome_list = idxwhere(
        pd.read_table(genome_filt_inpath, index_col="genome_id")
        .rename(str)
        .passes_filter.astype(bool)
    )

    presence = pd.read_table(gene_presence_inpath, index_col="gene_id")
    prevalence = (presence[genome_list].sum(1) + pseudo) / (len(genome_list) + pseudo)

    prevalence.to_csv(outpath, sep="\t", header=False)
