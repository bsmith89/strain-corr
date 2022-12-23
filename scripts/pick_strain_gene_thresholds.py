#!/usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    species_corr_path = sys.argv[1]
    n_marker_genes = int(sys.argv[2])
    strain_corr_path = sys.argv[3]
    strain_corr_quantile_strict = float(sys.argv[4])
    strain_corr_quantile_moderate = float(sys.argv[5])
    strain_corr_quantile_lenient = float(sys.argv[6])
    strain_depth_path = sys.argv[7]
    strain_depth_quantile = float(sys.argv[8])
    outpath = sys.argv[9]

    species_corr = pd.read_table(
        species_corr_path, names=["gene_id", "correlation"], index_col="gene_id"
    ).squeeze()
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
    # strain_corr = (
    #     xr.load_dataarray(strain_corr_path)
    #     .sel(species_id=species_id)
    #     .to_series()
    #     .unstack("strain")
    # )
    # Align genes and strains.
    gene_list = species_corr.index.to_list()
    strain_corr = strain_corr.reindex(gene_list, fill_value=0)
    strain_depth = strain_depth.reindex(gene_list, fill_value=0)

    # Select highly correlated species marker genes.
    species_gene_hit = species_corr.sort_values(ascending=False).head(n_marker_genes).index

    # Calculate the strain correlation threshold for each strain at which strain_corr_quantile
    # of the species genes (defined as those passing the species_corr_threshold)
    # are also assigned to the strain.
    strain_corr_threshold_strict = strain_corr.loc[species_gene_hit].quantile(
        strain_corr_quantile_strict
    )
    strain_corr_threshold_moderate = strain_corr.loc[species_gene_hit].quantile(
        strain_corr_quantile_moderate
    )
    strain_corr_threshold_lenient = strain_corr.loc[species_gene_hit].quantile(
        strain_corr_quantile_lenient
    )
    strain_depth_threshold_low = strain_depth.loc[species_gene_hit][
        strain_corr > strain_corr_threshold_moderate
    ].quantile(strain_depth_quantile)
    strain_depth_threshold_high = strain_depth.loc[species_gene_hit][
        strain_corr > strain_corr_threshold_moderate
    ].quantile(1 - strain_depth_quantile)

    out = pd.DataFrame(
        dict(
            correlation_strict=strain_corr_threshold_strict,
            correlation_moderate=strain_corr_threshold_moderate,
            correlation_lenient=strain_corr_threshold_lenient,
            depth_low=strain_depth_threshold_low,
            depth_high=strain_depth_threshold_high,
        )
    )
    out.to_csv(outpath, sep="\t")
