#!/usr/bin/env python3
# USAGE: {input.script} {input.seed} {input.n_samples} {input.samples} {output}


import sys
import pandas as pd
import numpy as np


if __name__ == "__main__":
    seed = int(sys.argv[1])
    n_samples = int(sys.argv[2])
    inpath = sys.argv[3]
    outpath = sys.argv[4]

    np.random.seed(seed)

    indata = pd.read_table(inpath)
    result = indata.sample(frac=1).groupby("strain").head(n_samples)

    result.to_csv(outpath, sep="\t", index=False)
