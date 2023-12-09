#!/usr/bin/env python3

import pandas as pd
import sys

if __name__ == "__main__":
    meta_inpath = sys.argv[1]
    genome_to_species_inpath = sys.argv[2]
    pos_inpath = sys.argv[3]
    genes_inpath = sys.argv[4]
    species = sys.argv[5]
    min_completeness = float(sys.argv[6])
    max_contamination = float(sys.argv[7])
    min_positions = int(sys.argv[8])
    outpath = sys.argv[9]

    npositions = pd.read_table(pos_inpath, index_col="genome_id")
    genes = pd.read_table(genes_inpath, index_col="gene_id").rename_axis(
        columns="genome_id"
    )
    assert genes.isin([0, 1]).values.all()
    genome_to_species = pd.read_table(
        genome_to_species_inpath, index_col="genome", dtype=str,
    ).species
    genome_meta = (
        pd.read_table(
            meta_inpath,
            names=[
                "Genome",
                "Genome_type",
                "Length",
                "N_contigs",
                "N50",
                "GC_content",
                "Completeness",
                "Contamination",
                "rRNA_5S",
                "rRNA_16S",
                "rRNA_23S",
                "tRNAs",
                "Genome_accession",
                "Species_rep",
                "Lineage",
                "Sample_accession",
                "Study_accession",
                "Country",
                "Continent",
                "FTP_download",
                "_20",
                "_21",
            ],
            index_col=["Genome_accession"],
        )
        .assign(
            species_id=genome_to_species,
        )[lambda x: x.species_id == species]
        .rename_axis(index="genome_id")
        # .rename(index=lambda x: "UHGG" + x[len("GUT_GENOME") :])
        .rename(
            columns={
                "Genome_type": "genome_type",
                "Completeness": "completeness",
                "Contamination": "contamination",
            }
        )
        .assign(
            passes_completeness=lambda x: (x.completeness / 100) >= min_completeness,
            passes_contamination=lambda x: (x.contamination / 100) <= max_contamination,
            num_genes=genes.sum(),
            num_geno_positions=npositions,
            passes_num_positions=lambda x: x.num_geno_positions >= min_positions,
            passes_filter=lambda x: x[
                [
                    "passes_completeness",
                    "passes_contamination",
                    "passes_num_positions",
                ]
            ].all(1),
        )
    )

    (
        genome_meta[
            [
                "genome_type",
                "completeness",
                "contamination",
                "num_genes",
                "num_geno_positions",
                "passes_completeness",
                "passes_contamination",
                "passes_num_positions",
                "passes_filter",
            ]
        ]
        # .astype(
        #     dict(passes_completeness=int, passes_contamination=int, passes_filter=int)
        # )
        .to_csv(outpath, sep="\t")
    )
