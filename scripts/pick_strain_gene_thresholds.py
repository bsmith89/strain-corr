#!/usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    species_gene_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_corr_quantile = float(sys.argv[3])
    strain_depth_path = sys.argv[4]
    strain_depth_quantile = float(sys.argv[5])
    outpath = sys.argv[6]

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

    gene_list = strain_corr.index.to_list()
    # TODO: Confirm that all entries in gene_list are also found in strain_corr and strain_depth.
    strain_corr = strain_corr.reindex(gene_list, fill_value=0)
    strain_depth = strain_depth.reindex(gene_list, fill_value=0)

    # Calculate the strain correlation threshold for each strain at which strain_corr_quantile
    # of the species genes (defined as those passing the species_corr_threshold)
    # are also assigned to the strain.
    strain_corr_threshold = strain_corr.loc[species_gene_hit].quantile(
        strain_corr_quantile
    )
    (_, strain_depth_threshold_low), (_, strain_depth_threshold_high) = (
        strain_depth.loc[species_gene_hit]
        .quantile([strain_depth_quantile, 1 - strain_depth_quantile])
        .iterrows()
    )
    out = pd.DataFrame(
        dict(
            correlation=strain_corr_threshold,
            depth_low=strain_depth_threshold_low,
            depth_high=strain_depth_threshold_high,
        ),
    )
    out.to_csv(outpath, sep="\t")
