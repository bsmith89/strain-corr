#!/usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    species_gene_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_corr_quantile_strict = float(sys.argv[3])
    strain_corr_quantile_moderate = float(sys.argv[4])
    strain_corr_quantile_lenient = float(sys.argv[5])
    strain_depth_path = sys.argv[6]
    strain_depth_quantile = float(sys.argv[7])
    outpath = sys.argv[8]

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
    (
        (_, strain_corr_threshold_strict),
        (_, strain_corr_threshold_moderate),
        (_, strain_corr_threshold_lenient),
    ) = (
        strain_corr.loc[species_gene_hit]
        .quantile(
            [
                strain_corr_quantile_strict,
                strain_corr_quantile_moderate,
                strain_corr_quantile_lenient,
            ]
        )
        .iterrows()
    )
    (_, strain_depth_threshold_low), (_, strain_depth_threshold_high) = (
        strain_depth.loc[species_gene_hit]
        .quantile([strain_depth_quantile, 1 - strain_depth_quantile])
        .iterrows()
    )
    out = pd.DataFrame(
        dict(
            correlation_strict=strain_corr_threshold_strict,
            correlation_moderate=strain_corr_threshold_moderate,
            correlation_lenient=strain_corr_threshold_lenient,
            depth_low=strain_depth_threshold_low,
            depth_high=strain_depth_threshold_high,
        ),
    )
    out.to_csv(outpath, sep="\t")
