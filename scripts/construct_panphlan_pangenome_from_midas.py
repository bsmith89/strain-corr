#!/usr/bin/env python3

import pandas as pd
import sys


if __name__ == "__main__":
    gene_info_path = sys.argv[1]
    cluster_info_path = sys.argv[2]
    centroid = sys.argv[3]
    fillna_gene_length = sys.argv[4]
    outpath = sys.argv[5]

    centroid_col = f"centroid_{centroid}"

    _gene_info = pd.read_table(gene_info_path).assign(
        genome_id=lambda x: x.gene_id.str.rsplit("_", n=1).str[0]
    )
    cluster_info = pd.read_table(cluster_info_path)
    gene_info = (
        _gene_info.join(
            cluster_info.set_index("centroid_99").centroid_99_length, on="centroid_99"
        )
        .fillna({"centroid_99_length": fillna_gene_length})
        .astype({"centroid_99_length": int})
        .assign(
            dummy_contig="dummy_contig",  # Probably doesn't matter, right?
            dummy_left=0,  # And right will be the length of the gene?
        )
    )

    gene_info[
        [
            centroid_col,
            "gene_id",
            "genome_id",
            "dummy_contig",
            "dummy_left",
            "centroid_99_length",
        ]
    ].to_csv(outpath, sep="\t", header=False, index=False)
