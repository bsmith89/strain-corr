#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    comm_path = sys.argv[1]
    mgen_subset_path = sys.argv[2]
    hits_path = sys.argv[3]
    comm = pd.read_table(comm_path, dtype={'sample': str, 'strain': str, 'community': float}, index_col=['sample', 'strain']).community.unstack()
    mgen_list = pd.read_table(mgen_subset_path, index_col='mgen_id').index.to_list()
    hits = pd.read_table(hits_path, index_col="gene_id").rename(columns=str)

    comm = comm.loc[list(set(mgen_list) & set(comm.index))]
    dominant_strain = comm.mean().sort_values(ascending=False).index[0]
    output = hits[[dominant_strain]]
    output.to_csv(sys.stdout, sep='\t')
