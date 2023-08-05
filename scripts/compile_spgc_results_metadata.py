#!/usr/bin/env python3

import pandas as pd
import sfacts as sf
from lib.pandas_util import read_list
import sys
import warnings


if __name__ == "__main__":
    species_gene_list_path = sys.argv[1]
    species_depth_path = sys.argv[2]
    strain_fit_path = sys.argv[3]
    sample_to_strain_path = sys.argv[4]
    strain_gene_path = sys.argv[5]
    outpath = sys.argv[6]

    # Load data
    species_gene_list = read_list(species_gene_list_path)
    species_depth = pd.read_table(
        species_depth_path, names=["sample", "depth"], index_col="sample"
    ).depth
    strain_fit = sf.World.load(strain_fit_path)
    sample_to_strain = pd.read_table(sample_to_strain_path, index_col="sample").strain
    strain_gene = pd.read_table(strain_gene_path, index_col="gene_id").rename_axis(
        columns="strain"
    )

    # Short-circuit on empty sample_to_strain table:
    # Short-circuit on empty data:
    emptyness = {
        k: v.empty
        for k, v in {
            "sample_to_strain": sample_to_strain,
            "species_depth": species_depth,
            "strain_gene": strain_gene,
        }.items()
    }
    if any(emptyness.values()):
        warnings.warn(
            f"Writing empty output because one or more input tables were empty\n{emptyness}"
        )
        with open(outpath, "w") as f:
            print(
                "strain",
                "num_sample",
                "max_depth",
                "sum_depth",
                "species_gene_frac",
                "num_genes",
                "strain_metagenotype_entropy",
                sep="\t",
                file=f,
            )
        exit(0)

    # Strain genotype entropy
    strain_total_mgtp = sf.Metagenotype(
        strain_fit.metagenotype.data.sel(sample=sample_to_strain.index)
        .groupby(sample_to_strain.to_xarray())
        .sum()
        .rename(strain="sample")
    )
    strain_metagenotype_entropy = strain_total_mgtp.entropy()

    # Species genes
    species_gene_frac = strain_gene.reindex(species_gene_list, fill_value=0).mean(0)

    # Total number of genes
    num_genes = strain_gene.sum(0)

    # Compile data
    strain_sample_meta = (
        species_depth.groupby(sample_to_strain)
        .apply(
            lambda x: pd.Series(
                dict(num_sample=len(x), max_depth=x.max(), sum_depth=x.sum())
            )
        )
        .unstack()
        .rename(lambda x: str(int(x)))  # Strain IDs get annoyingly retyped as floats.
        .assign(
            num_sample=lambda x: x.num_sample.astype(
                int
            ),  # Strain counts get annoyingly retyped as floats.
            species_gene_frac=species_gene_frac,
            num_genes=num_genes,
            strain_metagenotype_entropy=strain_metagenotype_entropy,
        )
    )

    # Write output
    strain_sample_meta.to_csv(outpath, sep="\t")