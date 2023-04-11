#!/usr/bin/env python3
# USAGE: {input.script} {input.species_gene} {input.gene_depth} {params.trim_frac} {output.species_depth}

import sys
import xarray as xr
from lib.util import info
from scipy.stats import trim_mean

if __name__ == "__main__":
    species_gene_inpath = sys.argv[1]
    gene_depth_inpath = sys.argv[2]
    trim_frac = float(sys.argv[3])
    species_depth_outpath = sys.argv[4]

    info("Loading species core gene list.")
    with open(species_gene_inpath) as f:
        species_gene_hit = []
        for line in f:
            species_gene_hit.append(line.strip())
    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_inpath)
    species_gene_depth = gene_depth.reindex(gene_id=species_gene_hit, fill_value=0)
    info("Calculating mean depth of samples.")
    species_mean_depth = xr.apply_ufunc(
        trim_mean,
        species_gene_depth,
        kwargs=dict(proportiontocut=trim_frac),
        input_core_dims=[["gene_id", "sample"]],
        output_core_dims=[["sample"]],
    )
    info("Writing sample mean depth.")
    species_mean_depth.to_series().to_csv(species_depth_outpath, sep="\t", header=False)
