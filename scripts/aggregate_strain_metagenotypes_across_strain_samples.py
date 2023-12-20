#!/usr/bin/env python3

import sys
import pandas as pd
import sfacts as sf

EMPTY_METAGENOTYPE = sf.data.Metagenotype(
    pd.DataFrame([], columns=["sample", "position", "allele", "tally"])
    .set_index(["sample", "position", "allele"])
    .tally.to_xarray()
)


if __name__ == "__main__":
    sample_to_strain_inpath = sys.argv[1]
    mgtp_inpath = sys.argv[2]
    outpath = sys.argv[3]

    mgtp = sf.data.Metagenotype.load(mgtp_inpath)
    sample_to_strain = pd.read_table(
        sample_to_strain_inpath, names=["sample", "strain"], index_col="sample"
    ).strain

    if sample_to_strain.empty:
        # NOTE: (2023-12-13) This edge case is necessary
        # only because some empty xr.DataArray objects
        # don't like to dump to NetCDF.
        observed_mgtp = EMPTY_METAGENOTYPE
    else:
        observed_mgtp = sf.Metagenotype(
            mgtp.data.sel(sample=sample_to_strain.index)
            .groupby(sample_to_strain.to_xarray())
            .sum()
            .rename(strain="sample")
        )

    observed_mgtp.dump(outpath)
