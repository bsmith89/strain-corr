#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info
from lib.pandas_util import read_table_lz4_filter, idxwhere
import numpy as np
from tqdm import tqdm


if __name__ == "__main__":
    gene_info_path = sys.argv[1]
    outpath = sys.argv[2]
    sample_paths = dict((arg.split("=", maxsplit=1) for arg in sys.argv[3:]))

    info("Loading gene info.")
    gene_info = pd.read_table(gene_info_path, index_col="gene_id")

    info("Pre-calculating gene lists.")
    species_gene_set = set(gene_info.index)
    species_gene_list = list(species_gene_set)
    num_output_genes = len(species_gene_list)

    data = np.empty((len(sample_paths), num_output_genes))
    info(f"Loading samples into matrix with final dims: {data.shape}.")
    for i, (sample, path) in tqdm(list(enumerate(sample_paths.items()))):
        _filt = lambda i, s: (i == 0) or (
            s.split(maxsplit=1, sep="\t")[0] in species_gene_set
        )
        d = read_table_lz4_filter(path, filt=_filt, index_col="gene_id").tally
        d = d.reindex(species_gene_list, fill_value=0)
        data[i] = d.values
    data = (
        pd.DataFrame(data, index=list(sample_paths), columns=species_gene_list)
        .rename_axis(index="sample", columns="gene_id")
        .astype(int)
    )
    data = data[idxwhere(data.sum() > 0)].T

    info("Writing output.")
    data.to_csv(outpath)
