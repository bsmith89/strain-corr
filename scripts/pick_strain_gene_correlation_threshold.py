#!/usr/bin/env python3

import pandas as pd
from lib.pandas_util import align_indexes, idxwhere
import sys

if __name__ == "__main__":
    species_corr_path = sys.argv[1]
    strain_corr_path = sys.argv[2]
    strain_depth_path = sys.argv[3]
    gene_meta_path = sys.argv[4]
    aggregate_genes_by = sys.argv[5]
    species_corr_threshold = float(sys.argv[6])
    strain_corr_quantile = float(sys.argv[7])
    trim_frac = float(sys.argv[8])
    outpath = sys.argv[9]

    species_corr_agg = pd.read_table(
        species_corr_path, names=["gene_id", "correlation"], index_col="gene_id"
    ).squeeze()
    strain_corr = (
        pd.read_table(strain_corr_path, index_col=["gene_id", "strain"])
        .squeeze()
        .unstack(fill_value=0)
    )
    gene_meta = (
        pd.read_table(gene_meta_path)
        .set_index("centroid_99", drop=False)
        .rename_axis(index="gene_id")
    )
    strain_depth = (
        pd.read_table(strain_depth_path, index_col=["gene_id", "strain"])
        .squeeze()
        .unstack()
    )

    # Align genes and strains.
    strain_corr, strain_depth = align_indexes(
        *align_indexes(strain_corr, strain_depth), axis="columns"
    )

    # Aggregate data.
    strain_corr_agg = strain_corr.groupby(gene_meta[aggregate_genes_by]).max()
    strain_depth_agg = strain_depth.groupby(gene_meta[aggregate_genes_by]).sum()
    # Align aggs.
    species_corr_agg, strain_corr_agg, strain_depth_agg = align_indexes(
        species_corr_agg, strain_corr_agg, strain_depth_agg
    )

    # Calculate the strain correlation threshold for each strain at which strain_corr_quantile
    # of the species genes (defined as those passing the species_corr_threshold)
    # are also assigned to the strain.
    species_agg_hit = species_corr_agg > species_corr_threshold
    number_species_agg_hit = species_agg_hit.sum()
    strain_selection_threshold = strain_corr_agg.reindex(
        idxwhere(species_agg_hit)
    ).quantile(1 - strain_corr_quantile)

    # Stats on aggregated gene hits (aggs where best hit is over threshold)
    strain_agg_hit = strain_corr_agg.gt(strain_selection_threshold)
    number_strain_agg_hit = strain_agg_hit.sum()
    total_depth_ratio_strain_agg_hit = (strain_depth_agg * strain_agg_hit).sum()

    # Stats on aggregated gene hits (aggs where best hit is over threshold) that were also species hits
    strain_species_agg_hit = (strain_agg_hit.T & species_agg_hit).T
    number_strain_species_agg_hit = strain_species_agg_hit.sum()
    total_depth_ratio_species_agg_hit = (strain_depth_agg.T * species_agg_hit).T.sum()
    total_depth_ratio_strain_species_agg_hit = (
        strain_depth_agg * strain_species_agg_hit
    ).sum()
    mean_depth_ratio_species_agg_hit = (strain_depth_agg[species_agg_hit]).mean()
    frac_total_depth_ratio_species_agg_hit = (
        total_depth_ratio_strain_species_agg_hit / total_depth_ratio_species_agg_hit
    )

    # Stats on gene hits (not aggregated)
    strain_gene_hit = strain_corr.gt(strain_selection_threshold)
    agg_depth_ratio_strain_gene_hit = (
        strain_depth[strain_gene_hit]
        .groupby(gene_meta[aggregate_genes_by])
        .sum()
        .reindex(species_corr_agg.index)
        .fillna(0)
    )
    agg_tally_strain_gene_hit = (
        strain_depth[strain_gene_hit]
        .notna()
        .groupby(gene_meta[aggregate_genes_by])
        .sum()
        .reindex(species_corr_agg.index)
        .fillna(0)
    )
    number_strain_gene_hit = agg_tally_strain_gene_hit.sum()
    total_depth_ratio_strain_gene_hit = agg_depth_ratio_strain_gene_hit.sum()

    # Stats on gene hits (not aggregated) that were also species hits
    number_strain_species_gene_hit = agg_tally_strain_gene_hit[species_agg_hit].sum()
    mean_agg_depth_ratio_strain_species_gene_hit = agg_depth_ratio_strain_gene_hit[
        species_agg_hit
    ].mean()

    out = pd.DataFrame(
        dict(
            strain_selection_threshold=strain_selection_threshold,
            number_strain_agg_hit=number_strain_agg_hit,
            total_depth_ratio_strain_agg_hit=total_depth_ratio_strain_agg_hit,
            number_species_agg=number_species_agg_hit,
            number_strain_species_agg_hit=number_strain_species_agg_hit,
            mean_depth_ratio_species_agg_hit=mean_depth_ratio_species_agg_hit,
            number_strain_gene_hit=number_strain_gene_hit,
            total_depth_ratio_strain_gene_hit=total_depth_ratio_strain_gene_hit,
            number_strain_species_gene_hit=number_strain_species_gene_hit,
            mean_agg_depth_ratio_strain_species_gene_hit=mean_agg_depth_ratio_strain_species_gene_hit,
        )
    )
    out.to_csv(outpath, sep="\t")
