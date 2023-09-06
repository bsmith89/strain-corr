#!/usr/bin/env python3

import xarray as xr
import sfacts as sf
import sys
from scipy.spatial.distance import pdist
import numpy as np
import pickle

def genotype_dissimilarity_transformed_values(x, y):
    dist = np.abs((x - y) / 2)
    weight = np.abs(x * y)
    wmean_dist = (weight * dist).sum() / weight.sum()
    # NOTE: Why not finish up by powering it by (1 / q)?
    # I don't do this part because it loses the city-block distance
    # interpretation when x and y are both discrete (i.e. one of {0, 1}).

    # While the basic function is undefined where weight.sum() == 0
    # (and this is only true when one of x or y is always exactly 0.5 at every
    # index),
    # the limit approaches the same value from both directions.
    # We therefore redefine the dissimilarity as a piecewise function,
    # but one that is nonetheless everywhere smooth and defined.
    if np.isnan(wmean_dist):
        return dist.mean()
    return wmean_dist

if __name__ == "__main__":
    spgc_agg_mgtp_inpath = sys.argv[1]
    ref_geno_inpath = sys.argv[2]
    outpath = sys.argv[3]

    spgc_est_mgtp = sf.Metagenotype.load(spgc_agg_mgtp_inpath)
    spgc_est_mgtp.data["sample"] = spgc_est_mgtp.data["sample"].to_series().apply(str)
    spgc_est_geno = spgc_est_mgtp.to_estimated_genotype(pseudo=0)
    spgc_est_mgtp.sizes

    # "Re-calculated GT-Pro genotype"
    ref_geno = sf.data.Genotype(
        xr.load_dataarray(ref_geno_inpath).fillna(0.5).rename({"genome_id": "strain"}).T
    )

    geno = sf.Genotype.concat(
        dict(
            ref=ref_geno.sel(position=spgc_est_geno.position),
            spgc=spgc_est_geno,
        ),
        dim="strain",
        rename=False,
    ).mlift("fillna", 0.5)
    geno_trsfm = sf.math.genotype_binary_to_sign(geno.values).astype(np.float32)
    geno_cdmat = pdist(geno_trsfm, metric=genotype_dissimilarity_transformed_values)
    labels = geno.strain.values
    with open(outpath, 'wb') as f:
        pickle.dump(dict(cdmat=geno_cdmat, labels=labels), file=f)
