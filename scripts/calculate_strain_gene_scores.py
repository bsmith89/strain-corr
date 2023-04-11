#!/usr/bin/env python3

import pandas as pd
import sys
from statsmodels.distributions.empirical_distribution import ECDF


if __name__ == "__main__":
    species_gene_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_depth_path = sys.argv[3]
    corr_outpath = sys.argv[4]
    depth_outpath = sys.argv[5]

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
    strain_corr = strain_corr.reindex(gene_list, fill_value=0)
    strain_depth = strain_depth.reindex(gene_list, fill_value=0)
    strain_corr_quantile = strain_corr.apply(lambda x: ECDF(x.reindex(species_gene_hit, fill_value=0))(x))
    strain_depth_quantile = strain_depth.apply(lambda x: ECDF(x.reindex(species_gene_hit, fill_value=0))(x))

    strain_corr_quantile.to_csv(corr_outpath, sep="\t")
    strain_depth_quantile.to_csv(depth_outpath, sep="\t")
