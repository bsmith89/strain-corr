#!/usr/bin/env python3
import sys
import pandas as pd

if __name__ == "__main__":
    inpath = sys.argv[1]
    agg_name = sys.argv[2]
    strain_name = sys.argv[3]
    outpath = sys.argv[4]

    strain_gene_x_agg = pd.read_table(inpath)
    result = (
        (strain_gene_x_agg.groupby(agg_name).apply(len) > 0)
        .astype(int)
        .to_frame(name=strain_name)
    ).rename_axis(index="gene_id")
    result.to_csv(outpath, sep="\t")
