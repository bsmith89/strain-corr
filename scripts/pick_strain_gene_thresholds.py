#!/usr/bin/env python3

import pandas as pd

# from lib.pandas_util import idxwhere
import sys

if __name__ == "__main__":
    species_gene_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_corr_quantile = float(sys.argv[3])
    min_corr = float(sys.argv[4])
    max_corr = float(sys.argv[5])
    strain_depth_path = sys.argv[6]
    strain_depth_quantile = float(sys.argv[7])
    min_depth = float(sys.argv[8])
    max_depth = float(sys.argv[9])
    outpath = sys.argv[10]

    # Select highly correlated species marker genes.
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

    # gene_list = strain_corr.index.to_list()
    # TODO: Confirm that all entries in gene_list are also found in strain_corr and strain_depth.
    strain_corr = strain_corr.reindex(species_gene_hit, fill_value=0)
    strain_depth = strain_depth.reindex(species_gene_hit, fill_value=0)

    strain_gene_mask = (strain_corr > min_corr) & (strain_depth > min_depth)

    # Calculate the strain correlation threshold for each strain at which strain_corr_quantile
    # of the species genes (defined as those passing the species_corr_threshold)
    # are also assigned to the strain.
    strain_corr_threshold = strain_corr[strain_gene_mask].quantile(
        strain_corr_quantile
    )
    (_, strain_depth_threshold_low), (_, strain_depth_threshold_high) = (
        strain_depth[strain_gene_mask]
        .quantile([strain_depth_quantile, 1 - strain_depth_quantile])
        .iterrows()
    )

    # Clip values above or below max/min thresholds.
    strain_corr_threshold = strain_corr_threshold.clip(min_corr, max_corr)
    strain_depth_threshold_low = strain_depth_threshold_low.clip(min_depth, max_depth)

    out = pd.DataFrame(
        dict(
            correlation=strain_corr_threshold,
            depth_low=strain_depth_threshold_low,
            depth_high=strain_depth_threshold_high,
        ),
    )
    out.to_csv(outpath, sep="\t")
