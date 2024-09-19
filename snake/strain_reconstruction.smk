# FIXME: Rename a bunch of files to match the */midasdb_v15.gene75.* format.
rule select_species_core_genes_from_reference_by_filtered_set_prevalence:
    output:
        species_gene="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.spgc_specgene-ref-filt-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes2.py",
        prevalence="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.prevalence.tsv",
    params:
        threshold=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.prevalence} {params.threshold} {output}"


# NOTE: (2023-12-01) This rule is the first step in implementing the
# new packaged SPGC pipeline.
rule export_gene_depth_table_from_netcdf:
    output:
        "{stem}.depth2.tsv.gz",
    input:
        script="scripts/export_gene_depth_table_from_netcdf.py",
        depth="{stem}.depth2.nc",
    shell:
        "{input.script} {input.depth} {output}"


rule calculate_species_depth_directly:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
    input:
        depth="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.depth2.tsv.gz",
        species_genes="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.spgc_specgene-{specgene_params}.species_gene.list",
    params:
        trim_frac_species_genes=0.15,
    conda:
        "conda/toolz4.yaml"
    shell:
        """
        spgc estimate_species_depth --trim-frac-species-genes {params.trim_frac_species_genes} {input.depth} {input.species_genes} {output}
        """


# See matching rule in snake/gtpro.smk
rule estimate_all_species_depths_in_group_spgc:
    output:
        "data/group/{group}/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.spgc_specgene-{specgene_params}.all_species_depth.tsv",
    input:
        species=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/{w.stem}.gene{w.centroidA}_{w.dbv}-{w.btp}-agg{w.centroidB}.spgc_specgene-{w.specgene_params}.species_depth.tsv"
            for species in config["species_group"][w.group]
        ],
    params:
        header="sample	species_id	depth",
        species_list=lambda w: config["species_group"][w.group],
        species_pattern="data/group/{group}/species/sp-$species/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                file={params.species_pattern}
                echo $file >&2
                awk -v species=$species -v OFS='\t' '{{print $1,species,$2}}' $file
            done
        ) > {output}
        """
