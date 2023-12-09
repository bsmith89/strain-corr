#!/usr/bin/env python3
# "{input.script} {input.spgc} {input.ref} {input.clust} {outpath}"

import sys
from lib.pandas_util import idxwhere
import pandas as pd


if __name__ == "__main__":
    spgc_meta_inpath = sys.argv[1]
    ref_meta_inpath = sys.argv[2]
    clust_inpath = sys.argv[3]
    outpath = sys.argv[4]

    # Load data
    spgc_meta = pd.read_table(spgc_meta_inpath, index_col="genome_id").rename(str)
    ref_meta = pd.read_table(ref_meta_inpath, index_col="genome_id").rename(
        columns={"passes_num_positions": "passes_geno_positions"}
    )  # NOTE (2023-09-19): This is a stopgap and should be fixed upstream.
    clust = pd.read_table(
        clust_inpath, names=["genome_id", "clust"], index_col="genome_id"
    ).clust

    ref_list = idxwhere(ref_meta.passes_filter)

    # Compile metadata
    data = pd.concat([spgc_meta, ref_meta]).assign(clust=clust)

    # Write output
    data.to_csv(outpath, sep="\t")
