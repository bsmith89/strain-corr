#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info
from lib.pandas_util import read_table_lz4
from tqdm import tqdm


if __name__ == "__main__":
    outpath = sys.argv[1]
    sample_paths = dict((arg.split("=", maxsplit=1) for arg in sys.argv[2:]))

    data = {}
    for i, (sample, path) in tqdm(list(enumerate(sample_paths.items()))):
        d = read_table_lz4(path, index_col='cluster_75_id')
        data[sample] = d.total_depth
    data = (
        pd.DataFrame(data)
        .rename_axis(columns="sample", index="gene_id")
        .fillna(0)
        .astype(int)
    )

    info("Writing output.")
    data.to_csv(outpath)
