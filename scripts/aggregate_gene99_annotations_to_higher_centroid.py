#!/usr/bin/env python3
# "{input.script} {input.annot} {input.clust} {params.agg} {output}"

import sys
import pandas as pd


def vote_on_annotations(df, annot_col):
    num_genes = len(df.gene_id.unique())
    counts = df[annot_col].value_counts()
    return counts / num_genes


if __name__ == "__main__":
    annot_inpath = sys.argv[1]
    gene_info_inpath = sys.argv[2]
    agg_column = sys.argv[3]
    outpath = sys.argv[4]

    annot_table = pd.read_table(annot_inpath)
    assert annot_table.columns[0] == "centroid_99"
    annot_name = annot_table.columns[1]
    gene_agg = pd.read_table(gene_info_inpath)[["gene_id", "centroid_99", agg_column]]
    all_votes = pd.merge(annot_table, gene_agg, on="centroid_99")
    agg_annot_votes = all_votes.groupby(agg_column).apply(
        vote_on_annotations, annot_name
    )

    result = (
        agg_annot_votes[lambda x: x >= 0.5]
        .rename_axis(index=(agg_column, annot_name))
        .to_frame("fraction")
        .reset_index()[[agg_column, annot_name]]
    )
    result.to_csv(outpath, sep="\t", index=False)
