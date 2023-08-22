#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    results = []
    for arg in sys.argv[1:]:
        genome_id, inpath = arg.split("=")
        results.append(
                pd.read_table(inpath).assign(genome_id=genome_id)
                )
    pd.concat(results).to_csv(sys.stdout, sep="\t", index=False)

