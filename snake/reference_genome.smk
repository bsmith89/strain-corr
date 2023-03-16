rule link_reference_genome:
    output:
        "data/species/sp-{species}/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


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


rule combine_strain_genome_gtpro_data_loadable:
    output:
        "data/species/sp-{species}/strain_genomes.gtpro.tsv.bz2",
    input:
        strain=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.tiles-l100-o99.gtpro_parse.tsv.bz2"
            for strain in species_genomes(w.species)
        ],
    params:
        strain_list=lambda w: species_genomes(w.species),
        pattern="data/species/sp-{species}/genome/$strain.tiles-l100-o99.gtpro_parse.tsv.bz2",
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


# rule blastn_midasdb_uhgg_against_genome:
#     output:
#         "{stemA}/species/sp-{species}/genome/midas_uhgg_pangenome.{stemB}-blastn.tsv",
#     input:
#         query="data/species/sp-{species}/pangenome.centroids.fn",
#         subject="{stemA}/species/sp-{species}/genome/{stemB}.fn",
#     shell:
#         """
#         blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
#         """


rule blastn_genome_against_midasdb_uhgg:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastn.tsv",
    input:
        query="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
        subject="data/species/sp-{species}/pangenome.centroids.fn",
    threads: 1
    shell:
        """
        blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
        """


rule blastn_genome_against_genome:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.{stemB}-blastn.tsv",
    input:
        query="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
        subject="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
    threads: 1
    shell:
        """
        blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
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
    strain_list = idxwhere(config["genome"].species_id == species)
    # assert len(strain_list) > 0
    return strain_list


rule debug_species_genomes:
    output:
        "debug_{species}_genomes.flag",
    params:
        test=lambda w: species_genomes(w.species),
    shell:
        "echo {params.test}; false"


rule calculate_bitscore_ratio_of_orfs_and_pangenome_genes:
    output:
        "data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastn.bitscore_ratio-c{centroid}.tsv",
    input:
        script="scripts/calculate_bitscore_ratio_for_gene_matching.py",
        orf_x_orf="data/species/sp-{species}/genome/{stemB}.{stemB}-blastn.tsv",
        orf_x_midas="data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastn.tsv",
        midasdb=ancient("ref/midasdb_uhgg"),
    params:
        aggregate_genes_by=lambda w: {
            "99": "centroid_99",
            "95": "centroid_95",
            "90": "centroid_90",
            "85": "centroid_85",
            "80": "centroid_80",
            "75": "centroid_75",
        }[w.centroid],
        gene_clust="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        """
        {input.script} \
                {input.orf_x_orf} \
                {input.orf_x_midas} \
                {params.gene_clust} \
                {params.aggregate_genes_by} \
                {output}
        """


rule assign_matching_genes:
    output:
        "data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastn.gene_matching-c{centroid}-t{thresh}.tsv",
    input:
        ratio="data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastn.bitscore_ratio-c{centroid}.tsv",
    params:
        thresh=lambda w: int(w.thresh) / 100,
    shell:
        """
        awk -v OFS='\t' -v thresh='{params.thresh}' '(NR > 1) && ($3 > thresh) {{print $1,$2}}' {input} > {output}
        """


use rule run_bowtie_multi_species_dereplicated_pangenome as run_bowtie_multi_species_pangenome_on_reference_genome_tiles with:
    output:
        "data/group/{group}/species/sp-{species}/genome/{stem}.pangenomes{centroid}.bam",
    input:
        bt2_dir="data/group/{group}/r.proc.pangenomes{centroid}.bt2.d",  # Assumes `proc` infix.
        r1="data/species/sp-{species}/genome/{stem}.fq.gz",
        r2="data/species/sp-{species}/genome/{stem}.fq.gz",


# FIXME: This is also a hub rule, and extends another hub rule. Consider
# commenting it out.
# NOTE: The "samples" here are actually reference strains, not samples.
# Then the sample_pattern+sample_list trick is doing something extra tricky
# where the species and strain are combined together.
use rule build_new_pangenome_db as build_pangenome_db_on_reference_genome_tiles with:
    output:
        protected("data/group/{group}/ALL_STRAINS.{stem}.pangenomes{centroid}-mapq{q}.db"),
    input:
        samples=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/genome/{strain}.{w.stem}.pangenomes{w.centroid}.gene_mapping_tally-mapq{q}.tsv.lz4"
            for strain, species in config["genome"].species_id.items()
        ],
        genes="data/group/{group}/r.proc.pangenomes.gene_info.tsv.lz4",
        nlength="data/group/{group}/r.proc.pangenomes99.nlength.tsv",
    params:
        sample_pattern="data/group/{group}/species/$sample.{stem}.pangenomes{centroid}.gene_mapping_tally-mapq{q}.tsv.lz4",
        sample_list=lambda w: [
            f"sp-{species}/genome/{strain}"
            for strain, species in config["genome"].species_id.items()
        ],


