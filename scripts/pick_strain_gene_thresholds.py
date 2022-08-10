#!/usr/bin/env python3

import pandas as pd
from lib.pandas_util import idxwhere
import sys

if __name__ == "__main__":
    species_corr_path = sys.argv[1]
    species_corr_threshold = float(sys.argv[2])
    strain_corr_path = sys.argv[3]
    strain_corr_quantile = float(sys.argv[4])
    strain_depth_path = sys.argv[5]
    strain_depth_quantile = float(sys.argv[6])
    outpath = sys.argv[7]

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
    species_gene_hit = idxwhere(species_corr > species_corr_threshold)

    # Calculate the strain correlation threshold for each strain at which strain_corr_quantile
    # of the species genes (defined as those passing the species_corr_threshold)
    # are also assigned to the strain.
    strain_corr_threshold = strain_corr.loc[species_gene_hit].quantile(
        strain_corr_quantile
    )
    strain_depth_threshold_low = strain_depth.loc[species_gene_hit].quantile(
        strain_depth_quantile
    )
    strain_depth_threshold_high = strain_depth.loc[species_gene_hit].quantile(
        1 - strain_depth_quantile
    )

    out = pd.DataFrame(
        dict(
            correlation=strain_corr_threshold,
            depth_low=strain_depth_threshold_low,
            depth_high=strain_depth_threshold_high,
        )
    )
    out.to_csv(outpath, sep="\t")
