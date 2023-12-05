#!/usr/bin/env python3

import sys
import xarray as xr
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    data = xr.open_dataset(sys.argv[1])
    species_genes = idxwhere(data['is_core_gene'].to_pandas())
    with open(sys.argv[2], 'w') as f:
        for gene in species_genes:
            print(gene, file=f)
