rule combine_midasdb_reference_genome_gtpro_data_loadable:  # Hub-rule?
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gtpro.tsv.bz2",
    input:
        genome=lambda w: [
            f"data/species/sp-{w.species}/genome/midasdb_{w.dbv}/{genome}.norm.tiles-l500-o31.gtpro_parse.tsv.bz2"
            for genome in config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species]
        ],
    params:
        genome_list=lambda w: config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species],
        pattern="data/species/sp-{species}/genome/midasdb_{dbv}/$genome.norm.tiles-l500-o31.gtpro_parse.tsv.bz2",
    shell:
        """
        for genome in {params.genome_list}
        do
            path={params.pattern}
            ( \
                bzip2 -dc $path \
                | awk -v OFS='\t' -v strain=$genome -v species={wildcards.species} '$1==species {{print strain,$0}}' \
                | bzip2 -zc \
            )
        done > {output}
        """


rule convert_midasdb_species_gene_info_to_reference_genome_table_v20:
    output:
        "data/species/sp-{species}/gene{centroid}_v20.reference_copy_number.nc",
    input:
        script="scripts/convert_gene_info_to_genome_table.py",
        genes="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    shell:
        "{input.script} {input.genes} centroid_{wildcards.centroid} {output}"


# TODO: Move this into a new snakefile with all of the reference database
# work.
rule ref_gene_copy_number_to_presence_table:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
    input:
        script="scripts/gene_copy_number_nc_to_strain_gene_tsv.py",
        copies="data/species/sp-{species}/gene{centroid}_{dbv}.reference_copy_number.nc",
    shell:
        "{input.script} {input.copies} {output}"
