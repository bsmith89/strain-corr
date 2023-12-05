#!/usr/bin/env python3

import sys
import xarray as xr


if __name__ == "__main__":
    data = xr.open_dataset(sys.argv[1])
    data = data['species_depth'].to_pandas()
    data.to_csv(sys.argv[2], sep='\t', header=False)
