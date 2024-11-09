#!/usr/bin/env python3
"{input.script} {input.emapper} {input.gene_info} {input.cog_category} {output.meta} {output.cog_matrix}"

import pandas as pd
import sys
import numpy as np

EGGNOG_COLUMN_NAMES = (
    "query seed_ortholog evalue score eggNOG_OGs max_annot_lvl "
    "COG_category Description Preferred_name GOs EC KEGG_ko "
    "KEGG_inpathway KEGG_Module KEGG_Reaction KEGG_rclass BRITE KEGG_TC CAZy "
    "BiGG_Reaction PFAMs"
).split(" ")

if __name__ == "__main__":
    emapper_inpath = sys.argv[1]
    gene_info_inpath = sys.argv[2]
    cog_category_inpath = sys.argv[3]

    gene_meta_outpath = sys.argv[4]
    gene_x_cog_outpath = sys.argv[5]

    gene_info = pd.read_table(gene_info_inpath, index_col="gene_id", encoding="latin1")

    gene_nlength = gene_info.groupby("centroid_75").gene_length.mean()

    gene_annotation = (
        pd.read_table(
            emapper_inpath,
            index_col="#query",
        )
        .rename_axis(index="gene_id")
        .replace({"-": np.nan})
    )
    gene_annotation = gene_nlength.to_frame("nlength").join(gene_annotation)

    # COG Category assigned by emapper.
    gene_x_cog_category1 = (
        gene_annotation.COG_category.fillna("-")
        .apply(list)
        .explode()[lambda x: x != "-"]
    )
    gene_x_cog_category1

    # COG Category from alternative source.
    cog_x_category = pd.read_table(
        cog_category_inpath,
        names=["cog", "cog_category", "description", "short_name", "_4", "_5", "_6"],
        index_col="cog",
        encoding="latin1",
    ).cog_category

    gene_x_cog = (
        gene_annotation.eggNOG_OGs.fillna("")
        .str.split(",")
        .explode()[lambda x: x.str.startswith("COG")]
        .str.split("@")
        .str[0]
    )
    gene_x_cog_category2 = gene_x_cog.map(cog_x_category).dropna().apply(list).explode()

    gene_x_cog_category = (
        pd.concat(
            [
                gene_x_cog_category1,
                gene_x_cog_category2,
            ]
        )
        .reset_index()
        .drop_duplicates()
    )

    gene_x_cog_category.columns = ["gene_id", "cog_category"]
    gene_x_cog_category = gene_x_cog_category.set_index("gene_id").cog_category[
        lambda x: ~x.isin(["S"])
    ]
    no_cog_category_gene_list = list(
        set(gene_annotation.index) - set(gene_x_cog_category.index)
    )
    gene_x_cog_category = pd.concat(
        [
            gene_x_cog_category,
            pd.Series("no_category", index=no_cog_category_gene_list),
        ]
    ).rename_axis("gene_id")

    gene_annotation = gene_annotation.assign(
        COG_category=gene_x_cog_category.groupby("gene_id").apply(lambda x: x.sum())
    )

    # Write outputs
    gene_annotation.to_csv(gene_meta_outpath, sep="\t")
    gene_x_cog_category.to_frame("cog_category").to_csv(gene_x_cog_outpath, sep="\t")
