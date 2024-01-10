#!/usr/bin/env python3
# "{input.script} {input.annot} {input.clust} {params.agg} {output}"

import sys
import pandas as pd


if __name__ == "__main__":
    annot_inpath = sys.argv[1]
    gene_info_inpath = sys.argv[2]
    agg_column = sys.argv[3]
    outpath = sys.argv[4]

    annot_table = pd.read_table(annot_inpath)
    assert annot_table.columns[0] == "centroid_99"
    annot_name = annot_table.columns[1]

    all_gene_agg = pd.read_table(gene_info_inpath)[
        ["gene_id", "centroid_99", agg_column]
    ]
    assert all_gene_agg.gene_id.is_unique

    # Count the total number of genes per aggregate level (e.g. centroid_75).
    out_of = all_gene_agg[agg_column].value_counts().rename_axis(agg_column)

    # Pare-down the table before merging.
    gene_agg = all_gene_agg[lambda x: x.centroid_99.isin(annot_table.centroid_99)]
    votes_for = pd.merge(annot_table, gene_agg, on="centroid_99")[
        [agg_column, annot_name]
    ].value_counts()

    result = (
        (votes_for / out_of)[lambda x: x >= 0.5]
        .rename_axis(index=(agg_column, annot_name))
        .to_frame("fraction")
    )
    result.reset_index()[[agg_column, annot_name]].to_csv(
        outpath, sep="\t", index=False
    )
