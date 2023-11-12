# EEN Metagenomes


rule process_raw_zotus_table:
    output:
        "data/group/een/a.proc.zotu_counts.tsv",
    input:
        "raw/een-mgen/2023-09-25_deborah.haecker@tum-create.edu.sg/zOTUs-Table_byron_final2.tab",
    shell:
        "sed -e '$ d' -e 's:\s*$::' {input} > {output}"


use rule collect_filtering_metadata_for_hmp2_spgc_strains as collect_filtering_metadata_for_een_spgc_strains with:
    output:
        "data/group/een/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
        # "test.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/een/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        sample_to_strain="data/group/een/{stem}.spgc_ss-{ss}.strain_samples.tsv",
        mgen="meta/een-mgen/sample.list",
    params:
        min_species_genes_frac=90 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=1 / 1000,
        min_geno_positions=100,


rule plot_een_strain_composition:
    output:
        "fig/group/een/species/sp-{species}/r.proc.gtpro.{sfacts}.gene{pang}.spgc_specgene-{specgene}.een_strain_plot.pdf",
    input:
        script="scripts/plot_een_strain_composition.py",
        sample="meta/een-mgen/sample.tsv",
        species_depths="data/group/een/r.proc.gene{pang}.spgc_specgene-{specgene}.species_depth.tsv",
        sfacts="data/group/een/species/sp-{species}/r.proc.gtpro.{sfacts}.world.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.sample} {input.species_depths} {input.sfacts} {wildcards.species} {output}"


# TODO: Move this rule to a more obvious location.
rule construct_group_files_for_een_select_species:
    output:
        touch("{stemA}/group/een/{stemB}.SELECT_SPECIES.flag"),
    input:
        lambda w: [
            f"{w.stemA}/group/een/species/sp-{species}/{w.stemB}"
            for species in config["species_group"]["een"]
        ],
