#!/usr/bin/env python3

import sys
import pandas as pd

BLASTP_HEADER_NAMES = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


if __name__ == "__main__":
    orf_x_orf_path = sys.argv[1]
    orf_x_pangenome_path = sys.argv[2]
    gene_clust_path = sys.argv[3]
    aggregate_genes_by = sys.argv[4]
    bitscore_ratio_thresh = float(sys.argv[5])
    outpath = sys.argv[6]

    gene_cluster = (
        pd.read_table(gene_clust_path)
        .set_index("centroid_99", drop=False)
        .rename_axis(index="gene_id")
    )

    orf_x_orf = pd.read_table(orf_x_orf_path, names=BLASTP_HEADER_NAMES)
    # Identify the best possible bitscore for each ORF (we can assume this is self-alignment.)
    max_bitscore = orf_x_orf.groupby(["qseqid"]).bitscore.max()
    orf_x_pangenome = (
        pd.read_table(orf_x_pangenome_path, names=BLASTP_HEADER_NAMES)
        # Compare each alignment's bitscore to the maximum possible.
        .assign(bitscore_ratio=lambda x: x.bitscore / x.qseqid.map(max_bitscore))
        .assign(centroid=lambda x: x.sseqid.map(gene_cluster[aggregate_genes_by]))
        # Score each pairing of ORFs by pangome-centroids by their best bitscore ratios.
        .groupby(["qseqid", "centroid"])
        .bitscore_ratio.max()
        .rename_axis(["orf", "gene"])
    )
    # Write pairs that pass the bitscore ratio threshold.
    orf_x_pangenome[lambda x: x > bitscore_ratio_thresh].to_frame().drop(
        columns=["bitscore_ratio"]
    ).to_csv(outpath, sep="\t")
