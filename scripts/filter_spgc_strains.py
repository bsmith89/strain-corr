#!/usr/bin/env python3
# "{input.script} {input.meta} {input.spgc} {input.sample_to_strain} {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {params.max_log_gene_depth_ratio_std} {params.min_geno_positions} {output}"

import pandas as pd
import sys
import numpy as np
import scipy as sp
import scipy.stats
import xarray as xr

if __name__ == "__main__":
    meta_inpath = sys.argv[1]
    # spgc_inpath = sys.argv[2]
    # sample_to_strain_inpath = sys.argv[3]
    min_species_genes_frac = float(sys.argv[2])
    min_total_depth = float(sys.argv[3])
    gene_count_outlier_alpha = float(sys.argv[4])
    max_log_selected_gene_depth_ratio_std = float(sys.argv[5])
    min_geno_positions = int(sys.argv[6])
    outpath = sys.argv[7]

    # sample_to_strain = pd.read_table(
    #     sample_to_strain_inpath, names=["sample", "strain"], index_col="sample"
    # ).strain.astype(str)
    # # with open(sample_list_inpath) as f:
    # #     sample_list = [s.strip() for s in f]
    # sample_list = list(
    #     sample_to_strain.index
    # )  # FIXME: (2023-12-04) This is superfluous and an artifact of previously explicitly splitting out portions of a larger group.
    # strain_x_sample_list_count = (
    #     sample_to_strain.to_frame("genome_id")
    #     .assign(in_sample_list=lambda x: x.index.isin(sample_list))
    #     .value_counts()
    #     .unstack("in_sample_list", fill_value=0)
    #     .reindex(columns=[False, True], fill_value=0)
    #     .rename(
    #         columns={
    #             True: "num_samples_in_sample_list",
    #             False: "num_samples_not_in_sample_list",
    #         }
    #     )
    # )

    strain_meta = pd.read_table(meta_inpath, index_col="genome_id").rename(str)

    # Fit a heavy-tailed model (t-dist with df=2) to the data
    observed_gene_counts = strain_meta[
        lambda x: x.species_gene_frac > min_species_genes_frac
    ].num_gene.values

    if len(observed_gene_counts) > 0:
        # Use the location and scale parameters to define the underlying copy
        # number distribution.
        _df, _loc, _scale = sp.stats.t.fit(observed_gene_counts, fix_df=2)
        _dist = sp.stats.norm(_loc, _scale)
        gene_count_cdf = strain_meta.num_gene.map(_dist.cdf)
        gene_count_score = np.minimum(gene_count_cdf, 1 - gene_count_cdf)
        # Filter genomes with too few or too many using a threshold set by
        # alpha.
        thresh_min_num_genes, thresh_max_num_genes = _dist.ppf(
            [gene_count_outlier_alpha, 1 - gene_count_outlier_alpha]
        )
    else:
        thresh_min_num_genes = np.nan
        thresh_max_num_genes = np.nan
        gene_count_score = np.nan

    strain_filter = strain_meta.assign(
        genome_type="SPGC",
        gene_count_score=lambda x: gene_count_score,
        passes_total_depth=lambda x: (x.sum_strain_depth >= min_total_depth),
        passes_species_gene_frac=lambda x: (
            x.species_gene_frac >= min_species_genes_frac
        ),
        passes_gene_count=lambda x: (x.gene_count_score >= gene_count_outlier_alpha),
        passes_log_selected_gene_depth_ratio_std=lambda x: (
            x.log_selected_gene_depth_ratio_std <= max_log_selected_gene_depth_ratio_std
        ),
        passes_geno_positions=lambda x: (x.num_geno_positions >= min_geno_positions),
        # passes_in_sample_list=lambda x: (x.num_samples_in_sample_list >= 1),
        passes_filter=lambda x: x[
            [
                "passes_total_depth",
                "passes_species_gene_frac",
                "passes_gene_count",
                "passes_geno_positions",
                # "passes_in_sample_list",
                "passes_log_selected_gene_depth_ratio_std",
            ]
        ].all(1),
    )

    # Write output.
    (
        strain_filter
        # .astype(
        #     dict(
        #         passes_total_depth=int,
        #         passes_species_gene_frac=int,
        #         passes_gene_count=int,
        #         passes_geno_positions=int,
        #         passes_filter=int,
        #     )
        # )
        .to_csv(outpath, sep="\t")
    )
