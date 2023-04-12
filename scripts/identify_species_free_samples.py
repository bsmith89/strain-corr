#!/usr/bin/env python3
# USAGE: {input.script} {params.thresh} {input.species_depth} {output}

import sys
import pandas as pd
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    thresh = float(sys.argv[1])
    inpath = sys.argv[2]
    outpath = sys.argv[3]
    depth = pd.read_table(inpath, names=["sample", "depth"], index_col=["sample"]).depth
    sample_list = idxwhere(depth < thresh)
    with open(outpath, "w") as f:
        for sample in sample_list:
            print(sample, file=f)
