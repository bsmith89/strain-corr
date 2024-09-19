# note: this rule takes the entire eggnog output and assigns raw annotations
# to the features. in a later step, i'll aggregate these annotations
# to the higher-level-centroid by voting or something.
# note: the */eggnog.tsv file is annotatins on c99s, not unique features.
# todo: rename all these scripts to "..._emapper_output_to_gene99_x...".
rule parse_midasdb_emapper_annotations_to_gene99_x_unit_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.emapper.gene99_x_{unit}.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_{unit}.py",
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
    shell:
        "{input.script} {input.emapper} {output}"


# Take annotations at the c99 level and combine them into centroidNN
# annotations using voting at the unique feature (gene_id) level.
rule aggregate_gene99_annotations_to_higher_centroid_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.emapper.gene{centroid}_x_{unit}.tsv",
    input:
        script="scripts/aggregate_gene99_annotations_to_higher_centroid.py",
        annot="data/species/sp-{species}/midasdb_v20.emapper.gene99_x_{unit}.tsv",
        clust="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    wildcard_constraints:
        centroid="95|90|85|80|75",  # Note no 99; this would potential result in a circular dependency.
    params:
        agg=lambda w: {
            95: "centroid_95",
            90: "centroid_90",
            85: "centroid_85",
            80: "centroid_80",
            75: "centroid_75",
        }[int(w.centroid)],
    resources:
        walltime_hr=10,
    shell:
        "{input.script} {input.annot} {input.clust} {params.agg} {output}"
