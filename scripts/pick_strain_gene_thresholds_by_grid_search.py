#!/usr/bin/env python3
"""

USAGE: {input.script} \
         {input.species_gene} \
         {input.strain_corr} \
         {input.strain_depth} \
         {params.resolution_corr} \
         {params.resolution_depth} \
         {params.alpha} \
         {output}
"""

import pandas as pd
import numpy as np
import sys
from itertools import product
from lib.util import info
from tqdm import tqdm

if __name__ == "__main__":
    species_gene_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_depth_path = sys.argv[3]
    corr_resolution = int(sys.argv[4])
    depth_resolution = int(sys.argv[5])
    alpha = float(sys.argv[6])
    outpath = sys.argv[7]

    # Select highly correlated species marker genes.
    info("Loading data.")
    with open(species_gene_path) as f:
        species_gene_hit = [line.strip() for line in f]

    strain_corr = (
        pd.read_table(strain_corr_path, index_col=["gene_id", "strain"])
        .squeeze()
        .unstack(fill_value=0)
    )
    strain_depth = (
        pd.read_table(strain_depth_path, index_col=["gene_id", "strain"])
        .squeeze()
        .unstack(fill_value=0)
    )

    info("Aligning data.")
    gene_list = list(
        set(species_gene_hit) | set(strain_corr.index) | set(strain_depth.index)
    )
    strain_list = list(set(strain_corr.columns) | set(strain_depth.columns))

    strain_corr = strain_corr.reindex(
        index=gene_list, columns=strain_list, fill_value=0
    )
    strain_depth = strain_depth.reindex(
        index=gene_list, columns=strain_list, fill_value=0
    )

    ref_strain_corr = strain_corr.reindex(species_gene_hit)
    ref_strain_depth = strain_depth.reindex(species_gene_hit)

    info(f"Filling {corr_resolution}-by-{depth_resolution} grid.")
    corr_grid = np.linspace(0, 1, num=corr_resolution)
    depth_grid = np.linspace(0, 1, num=depth_resolution)
    both_excluded = np.empty((corr_resolution, depth_resolution, len(strain_list)))
    refs_included = np.empty((corr_resolution, depth_resolution, len(strain_list)))
    for (i, corr_thresh), (j, depth_thresh) in tqdm(
        product(enumerate(depth_grid), enumerate(depth_grid))
    ):
        both_excluded[i, j, :] = (
            (strain_corr < corr_thresh) & (strain_depth < depth_thresh)
        ).sum()
        refs_included[i, j, :] = (
            (ref_strain_corr >= corr_thresh) & (ref_strain_depth > depth_thresh)
        ).sum()

    info("Identifying argmax for each strain.")
    score = np.log(both_excluded + 1) + alpha * np.log(refs_included + 1)
    correlation_threshold = []
    depth_threshold = []
    for i, strain in enumerate(strain_list):
        correlation_argmax, depth_argmax = (
            pd.DataFrame(score[:, :, i], index=corr_grid, columns=depth_grid)
            .stack()
            .idxmax()
        )
        correlation_threshold.append(correlation_argmax)
        depth_threshold.append(depth_argmax)

    out = pd.DataFrame(
        dict(
            correlation=correlation_threshold,
            depth_low=depth_threshold,
        ),
        index=strain_list,
    ).rename_axis(index="strain")
    out.to_csv(outpath, sep="\t")
