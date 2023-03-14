#!/usr/bin/env python3

import pandas as pd
import xarray as xr
import sys
from lib.util import info
from scipy.spatial.distance import cdist
from tqdm import tqdm
import numpy as np
from collections import defaultdict

if __name__ == "__main__":
    nospecies_path = sys.argv[1]
    partition_path = sys.argv[2]
    species_depth_path = sys.argv[3]
    gene_depth_path = sys.argv[4]
    corr_outpath = sys.argv[5]
    depth_outpath = sys.argv[6]

    info("Loading input data.")
    # Collect list of species-free samples.
    with open(nospecies_path) as f:
        no_species_samples = [line.strip() for line in f]

    # Build partition dictionary of {strain: [sample_list]}.
    strain_partition = defaultdict(lambda: [])
    with open(partition_path) as f:
        for line in f:
            strain, sample = line.strip().split("\t")
            strain_partition[strain].append(sample)

    info("Loading species depth.")
    species_depth = (
        pd.read_table(species_depth_path, names=["sample", "depth"], index_col="sample")
        .squeeze()
        .to_xarray()
    )
    info("Loading gene depth.")
    gene_depth = xr.load_dataarray(gene_depth_path).fillna(
        0
    )  # FIXME: Shouldn't be necessary.
    info("Aligning indexes.")
    shared_samples = set(species_depth.sample.values) & set(gene_depth.sample.values)
    # # NOTE: Dropping these lines because I align focal samples
    # # for each strain individually.
    # species_depth = species_depth.sel(sample=shared_samples)
    # gene_depth = gene_depth.sel(sample=shared_samples)
    assert not (
        set(no_species_samples) - shared_samples
    ), "Some no-species samples are missing from shared samples list."

    info("Iterating strains.")
    corr = {}
    depth = {}
    for strain, strain_pure_samples in tqdm(strain_partition.items()):
        strain_pure_samples = list(set(strain_pure_samples) & shared_samples)
        if strain_pure_samples:
            focal_samples = strain_pure_samples + no_species_samples
            corr[strain] = pd.Series(
                1
                - cdist(
                    species_depth.sel(sample=focal_samples).expand_dims(dict(_=1)),
                    gene_depth.sel(sample=focal_samples).transpose("gene_id", "sample"),
                    metric="cosine",
                )[0],
                index=gene_depth.gene_id,
            ).fillna(0)
            depth[strain] = (
                gene_depth.sel(sample=strain_pure_samples).sum("sample")
                / species_depth.sel(sample=strain_pure_samples).sum("sample")
            ).to_series()
        else:  # TODO: Is this a no-op? Won't all strains have a non-empty list?
            corr[strain] = pd.Series(np.nan, index=gene_depth.gene_id)
            depth[strain] = pd.Series(np.nan, index=gene_depth.gene_id)

    info("Compiling correlations table.")
    corr = (
        pd.DataFrame(corr)
        .rename_axis(index="gene_id", columns="strain")
        .stack()
        .rename("correlation")
    )
    info("Compiling depth table.")
    depth = (
        pd.DataFrame(depth)
        .rename_axis(index="gene_id", columns="strain")
        .stack()
        .rename("depth")
    )
    info("Writing outputs.")
    corr.to_csv(corr_outpath, sep="\t")
    depth.to_csv(depth_outpath, sep="\t")
