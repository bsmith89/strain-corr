rule link_project_reference_genome:
    output:
        "data/species/sp-{species}/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


rule link_project_reference_genome_no_species:
    output:
        "data/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


rule link_midasdb_reference_genome:
    output:
        "data/species/sp-{species}/genome/midasdb_{dbv}/{genome}.fn",
    input:
        "ref/midasdb_uhgg_{dbv}/mags/{species}/{genome}.fa",
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


rule normalize_genome_sequence:
    output:
        "{stem}.norm.fn",
    input:
        "{stem}.fn",
    shell:
        "sed '/^>/!s/[^ACGT]/N/g' {input} > {output}"


rule tile_reference_genome:
    output:
        "{stem}.tiles-l{length}-o{overlap}.fn",
    input:
        script="scripts/tile_fasta.py",
        fn="{stem}.fn",
    wildcard_constraints:
        genome=noperiod_wc,
        length=integer_wc,
        overlap=integer_wc,
    params:
        length=lambda w: int(w.length),
        overlap=lambda w: int(w.overlap),
    shell:
        "{input.script} {params.length} {params.overlap} {input.fn} > {output}"


rule genome_fasta_to_fastq:
    """
    Convert a FASTA formatted file into FASTQ.
    \
    Input/output patterns are limited to files found in */genome/* in order to
    prevent circular dependencies.
    \
    """
    output:
        "{stemA}/genome/{stemB}.fq.gz",
    input:
        "{stemA}/genome/{stemB}.fn",
    conda:
        "conda/seqtk.yaml"
    shell:
        "seqtk seq -F '#' {input} | gzip -c > {output}"


rule gtpro_species_lines_counts:
    output:
        "{stemA}/genome/{genome}.{stemB}.gtpro_species_tally.tsv",
    input:
        "{stemA}/genome/{genome}.{stemB}.gtpro_parse.tsv.bz2",
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        """
        bzcat {input} \
                | sed '1,1d' \
                | cut -f1 \
                | sort \
                | uniq -c \
                | awk -v OFS='\t' -v genome={wildcards.genome} '{{print genome,$2,$1}}' \
            > {output}
        """


rule combine_project_genome_gtpro_species_lines_counts_no_species:
    output:
        "data/group/{group}/strain_genomes.{stem}.gtpro_species_tally.tsv"
    input:
        strain=lambda w: [
            f"data/genome/{genome}.{w.stem}.gtpro_species_tally.tsv"
            for genome, species in group_genomes(w.group)
        ],
    shell:
        """
        cat {input.strain} > {output}
        """


rule combine_strain_genome_gtpro_data_loadable:  # Hub-rule?
    output:
        "data/group/{group}/species/sp-{species}/strain_genomes.gtpro.tsv.bz2",
    input:
        strain=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.norm.tiles-l500-o31.gtpro_parse.tsv.bz2"
            for strain in species_group_genomes(w.species, w.group)
        ],
    params:
        strain_list=lambda w: species_group_genomes(w.species, w.group),
        pattern="data/species/sp-{species}/genome/$strain.norm.tiles-l500-o31.gtpro_parse.tsv.bz2",
    shell:
        """
        for strain in {params.strain_list}
        do
            path={params.pattern}
            ( \
                bzip2 -dc $path \
                | awk -v OFS='\t' -v strain=$strain -v species={wildcards.species} '$1==species {{print strain,$0}}' \
                | bzip2 -zc \
            )
        done > {output}
        """


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


rule alias_midas_uhgg_pangenome_cds_new:
    output:
        "data/species/sp-{species}/midasdb.gene99_{dbv}.centroids.fn",
    input:
        fasta="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/centroids.ffn",
    shell:
        alias_recipe


localrules:
    alias_midas_uhgg_pangenome_cds_new,


rule blastn_genome_against_midasdb_uhgg_new:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.tsv",
    input:
        query="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
        subject="data/species/sp-{species}/midasdb.gene99_{dbv}.centroids.fn",
    threads: 1
    shell:
        """
        blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
        """


rule assign_matching_genes_based_on_best_blastn_hit:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best.tsv",
    input:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.tsv",
    shell:
        """
        sort -k1,1 -k12,12rn {input} | sort -k1,1 -u | cut -f1,2 > {output}
        """


rule aggreggate_top_blastn_hits_by_midasdb_centroid_v20:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_v20-blastn.gene_matching-best-c{centroid}.uhggtop-strain_gene.tsv",
    input:
        script="scripts/identify_strain_genes_from_top_blastn_hits.py",
        hit="{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_v20-blastn.gene_matching-best.tsv",
        gene_clust="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    params:
        aggregate_genes_by=lambda w: {
            "99": "centroid_99",
            "95": "centroid_95",
            "90": "centroid_90",
            "85": "centroid_85",
            "80": "centroid_80",
            "75": "centroid_75",
        }[w.centroid],
    shell:
        """
        {input.script} {input.gene_clust} {params.aggregate_genes_by} {input.hit} {output}
        """


use rule diamond_search_fa as blastp_midasdb_uhgg_against_genome with:
    output:
        "{stemA}/species/sp-{species}/genome/midas_uhgg_pangenome.{stemB}-blastp.tsv",
    input:
        fasta="data/species/sp-{species}/pangenome.centroids.tran.fa",
        db="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.dmnd",


use rule diamond_search_fa as blastp_genome_against_midasdb_uhgg with:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastp.tsv",
    input:
        fasta="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.fa",
        db="data/species/sp-{species}/pangenome.centroids.tran.dmnd",
    params:
        db=lambda w, input: parse_diamond_db_from_path(input.db),
        extra_diamond_blastp_args="--ultra-sensitive --max-hsps 10000",


use rule diamond_search_fa as reciprocal_blastp_genome with:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.{stemB}-blastp.tsv",
    input:
        fasta="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.fa",
        db="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.dmnd",


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
