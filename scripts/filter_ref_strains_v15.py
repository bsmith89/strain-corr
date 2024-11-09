#!/usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    meta_inpath = sys.argv[1]
    min_completeness = float(sys.argv[2])
    max_contamination = float(sys.argv[3])
    min_positions = int(sys.argv[4])
    outpath = sys.argv[5]

    genome_meta = pd.read_table(meta_inpath, index_col="genome_id").assign(
        passes_completeness=lambda x: (x.completeness / 100) >= min_completeness,
        passes_contamination=lambda x: (x.contamination / 100) <= max_contamination,
        passes_num_positions=lambda x: x.num_geno_positions >= min_positions,
        passes_filter=lambda x: x[
            [
                "passes_completeness",
                "passes_contamination",
                "passes_num_positions",
            ]
        ].all(1),
    )

    (
        genome_meta[
            [
                "genome_type",
                "completeness",
                "contamination",
                "num_genes",
                "num_geno_positions",
                "passes_completeness",
                "passes_contamination",
                "passes_num_positions",
                "passes_filter",
            ]
        ].to_csv(outpath, sep="\t")
    )
