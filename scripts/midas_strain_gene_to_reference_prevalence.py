#!/usr/bin/env python3
# USAGE: "{input.script} \
# {input.meta} \
# {input.genome_to_species} \
# {input.genes} \
# {wildcards.species} \
# {params.complete_thresh} \
# {params.contam_thresh} \
# {params.pseudo} \
# {output}"

import pandas as pd
import sys
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    meta_inpath = sys.argv[1]
    genome_to_species_inpath = sys.argv[2]
    genes_matrix_inpath = sys.argv[3]
    species_name = sys.argv[4]
    complete_thresh = float(sys.argv[5])
    contam_thresh = float(sys.argv[6])
    pseudo = float(sys.argv[7])
    outpath = sys.argv[8]

    genome_to_species = pd.read_table(
        genome_to_species_inpath,
        index_col="genome",
        dtype=str,
    ).species
    genes_matrix = pd.read_table(genes_matrix_inpath, index_col="gene_id")

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
        )[lambda x: x.species_id == species_name]
        .rename_axis(index="genome_id")
        .assign(
            completeness=lambda x: x.Completeness.astype(float),
            contamination=lambda x: x.Contamination.astype(float),
        )
        .assign(
            passes_completeness=lambda x: (x.completeness / 100) >= complete_thresh,
            passes_contamination=lambda x: (x.contamination / 100) <= contam_thresh,
            passes_filter=lambda x: x[
                [
                    "passes_completeness",
                    "passes_contamination",
                ]
            ].all(1),
        )
    )

    genome_list = idxwhere(genome_meta.passes_filter.astype(bool))

    gene_prevalence = (genes_matrix[genome_list].sum(1) + pseudo) / (
        len(genome_list) + pseudo
    )

    gene_prevalence.to_csv(outpath, sep="\t", header=False)
