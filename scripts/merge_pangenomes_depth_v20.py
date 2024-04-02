#!/usr/bin/env python3

# import pandas as pd
import sys
import pandas as pd
from lib.util import info
from lib.pandas_util import read_table_lz4_filter
import numpy as np
from tqdm import tqdm


if __name__ == "__main__":
    genes_info_path = sys.argv[1]
    groupby_col = sys.argv[2]
    outpath = sys.argv[3]
    sample_paths = dict((arg.split("=", maxsplit=1) for arg in sys.argv[4:]))

    info("Loading gene info.")
    gene_info = pd.read_table(genes_info_path, index_col="gene_id")
    gene_groupby = gene_info[groupby_col]
    info("Loading gene lengths.")
    gene_len = gene_info.gene_length

    info("Pre-calculating gene lists.")
    species_gene_set = set(gene_len.index)
    output_gene_list = sorted(gene_groupby.unique())
    num_output_genes = len(output_gene_list)

    data = np.empty((len(sample_paths), num_output_genes))
    info(f"Loading samples into matrix with final dims: {data.shape}.")
    for i, (sample, path) in tqdm(enumerate(sample_paths.items())):
        _filt = lambda i, s: (i == 0) or (
            s.split(maxsplit=1, sep="\t")[0] in species_gene_set
        )
        d = read_table_lz4_filter(path, filt=_filt, index_col="gene_id").tally
        sample_gene_list = list(set(d.index) & species_gene_set)
        d = d.loc[sample_gene_list] / gene_len.loc[sample_gene_list]
        d = d.groupby(gene_groupby.loc[sample_gene_list]).sum()
        d = d.reindex(output_gene_list, fill_value=0)
        data[i] = d.values
    data = pd.DataFrame(data, index=list(sample_paths), columns=output_gene_list).rename_axis(index="sample", columns="gene_id")

    info("Writing output.")
    data.stack().to_xarray().to_netcdf(outpath)
