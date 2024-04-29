# EEN Metagenomes


# TODO: rule concatenate_species_depths:


rule process_raw_zotus_table:
    output:
        "data/group/een/a.proc.zotu_counts.tsv",
    input:
        "raw/een-mgen/2023-09-25_deborah.haecker@tum-create.edu.sg/zOTUs-Table_byron_final2.tab",
    shell:
        "sed -e '$ d' -e 's:\s*$::' {input} > {output}"


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
