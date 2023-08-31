#!/usr/bin/env python3

import sys
import pandas as pd
import sfacts as sf


if __name__ == "__main__":
    sample_to_strain_inpath = sys.argv[1]
    mgtp_inpath = sys.argv[2]
    outpath = sys.argv[3]

    mgtp = sf.data.Metagenotype.load(mgtp_inpath)
    sample_to_strain = pd.read_table(
        sample_to_strain_inpath, index_col="sample"
    ).squeeze()

    observed_mgtp = sf.Metagenotype(
        mgtp.data.sel(sample=sample_to_strain.index)
        .groupby(sample_to_strain.to_xarray())
        .sum()
        .rename(strain="sample")
    )

    observed_mgtp.dump(outpath)
