rule link_reference_genome:
    output:
        "data/species/sp-{species}/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    shell:
        alias_recipe


rule tile_reference_genome:
    output:
        "{stem}.tiles-l{length}-o{overlap}.fn",
    input:
        script="scripts/tile_fasta.py",
        fn="{stem}.fn",
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
        "{stemA}/genome/{stemB}.fq",
    input:
        "{stemA}/genome/{stemB}.fn",
    conda:
        "conda/seqtk.yaml"
    shell:
        "seqtk seq -F '#' {input} > {output}"


rule alias_midas_uhgg_pangenome_cds:
    output:
        "data/species/sp-{species}/pangenome.centroids.fn",
    input:
        midas_download_flag="data/species/sp-{species}/download_species_midasdb_uhgg.flag",
    params:
        fasta="ref/midasdb_uhgg/pangenomes/{species}/centroids.ffn",
    shell:
        """
        ln -rs {params.fasta} {output}
        """


use rule diamond_search_fa as blastp_midasdb_uhgg_against_genome with:
    output:
        "{stemA}/species/sp-{species}/genome/midas_uhgg_pangenome.{stemB}-blastp.tsv",
    input:
        fasta="data/species/sp-{species}/pangenome.centroids.tran.fa",
        db="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.dmnd",
    # params:
    #     db=lambda w, input: parse_diamond_db_from_path(input.db),
    #     extra_diamond_blastp_args="--mid-sensitive",


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
    # params:
    #     db=lambda w, input: parse_diamond_db_from_path(input.db),
    #     extra_diamond_blastp_args="--mid-sensitive",


# rule blast_midasdb_uhgg_against_genome:
#     output:
#         "{stemA}/species/sp-{species}/genome/midas_uhgg_pangenome.{stemB}-blastx.tsv",
#     input:
#         db="{stemA}/species/sp-{species}/genome/{stemB}.cds.tran.dmnd",
#         fa="data/species/sp-{species}/pangenome.centroids.tran.fa",
#     params:
#         db_stem="{stemA}/species/sp-{species}/genome/{stemB}.cds.tran",
#     threads: 4
#     shell:
#         """
#         diamond blastp --threads {threads} --db {params.db_stem} --query {input.fa} > {output}
#         """


# rule blast_genome_against_midasdb_uhgg:
#     output:
#         "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastx.tsv",
#     input:
#         fa="{stemA}/species/sp-{species}/genome/{stemB}.cds.tran.fa",
#         db="data/species/sp-{species}/pangenome.centroids.tran.dmnd",
#     params:
#         db_stem="data/species/sp-{species}/pangenome.centroids.tran.dmnd",
#     threads: 4
#     shell:
#         """
#         diamond blastp --threads {threads} --db {params.db_stem} --query {input.fa} > {output}
#         """

def species_genomes(species):
    strain_list = idxwhere(config['genome'].species_id == species)
    assert len(strain_list) > 0
    return strain_list
