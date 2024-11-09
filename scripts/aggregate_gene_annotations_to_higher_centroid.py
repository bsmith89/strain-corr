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

    # Collect and check the gene annotations
    # NOTE: index_col != "gene_id" because it may not be unique;
    # this is a many-to-many mapping.
    annot_table = pd.read_table(annot_inpath)
    assert annot_table.columns[0] == "gene_id"
    annot_name = annot_table.columns[1]

    gene_agg = pd.read_table(gene_info_inpath, index_col="gene_id")[[agg_column]]

    agg_genes = set(gene_agg.index)
    annot_genes = set(annot_table.gene_id)
    shared_genes = list(agg_genes & annot_genes)
    annot_only_genes = annot_genes - agg_genes
    num_annot_only_genes = len(annot_only_genes)
    if num_annot_only_genes > 0:
        print(
            f"WARNING: Found {num_annot_only_genes} genes in annotations that were missing in the gene_info.",
            file=sys.stderr,
        )

    count_votes = (
        annot_table[["gene_id", annot_name]]
        .join(gene_agg.loc[shared_genes], on="gene_id")
        .dropna(subset=[agg_column])[[agg_column, annot_name]]
        .value_counts()
    )
    out_of = (
        gene_agg[agg_column].value_counts().rename("tally").rename_axis("centroid_75")
    )
    vote_share = count_votes / out_of

    result = (
        vote_share[lambda x: x >= 0.5]
        .to_frame("fraction")
        .reset_index()[[agg_column, annot_name]]
    )
    result.to_csv(outpath, sep="\t", index=False)
