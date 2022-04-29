#!/usr/bin/env python3

import xarray as xr
import seaborn as sns
import numpy as np
import scipy as sp
import scipy.stats
from tqdm import tqdm
from lib.util import is_empty
import pandas as pd
import sys
import math


if __name__ == '__main__':
    trim_frac = float(sys.argv[1])
    cvrg = {}
    pbar = tqdm(sys.argv[2:])
    for arg in pbar:
        species_id, path = arg.split('=')
        pbar.set_postfix({"species_id": species_id})
        species_id = int(species_id)
        species_data = xr.load_dataarray(path)
        if is_empty(species_data):
            continue
        cvrg[species_id] = (
            species_data
            .sel(species_id=species_id)
            .sum('allele')
            .to_pandas()
            .apply(lambda x: sp.stats.trim_mean(x, trim_frac), axis=1)
        )
    cvrg = pd.DataFrame(cvrg).rename_axis(index='mgen_id', columns='species_id')
    cvrg.stack()[lambda x: x > 0].rename('coverage').to_csv(sys.stdout, sep='\t')
