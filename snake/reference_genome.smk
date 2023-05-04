use rule start_shell as start_shell_eggnog with:
    conda:
        "conda/eggnog.yaml"


use rule start_shell as start_shell_dbcan with:
    conda:
        "conda/dbcan.yaml"


rule link_reference_genome:
    output:
        "data/species/sp-{species}/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


rule eggnog_mapper_translated_orfs:
    output:
        "{stem}.emapper.d/proteins.emapper.annotations",
        "{stem}.emapper.d/proteins.emapper.hits",
        "{stem}.emapper.d/proteins.emapper.seed_orthologs",
    input:
        fasta="{stem}.tran.fa",
        db="ref/eggnog_mapper_db",
    params:
        outdir="{stem}.emapper.d",
        tax_scope="auto",
        sensmode="more-sensitive",
        mapper="diamond",
    conda:
        "conda/eggnog.yaml"
    threads: 24
    resources:
        walltime_hr=240,
        mem_mb=20_000,
        pmem=20_000 // 24,
    shell:
        """
        tmpdir=$(mktemp -d)
        export EGGNOG_DATA_DIR={input.db}
        rm -rf {params.outdir}.temp {params.outdir}
        mkdir -p {params.outdir}.temp
        emapper.py \
                -m {params.mapper} \
                -i {input.fasta} \
                --itype proteins \
                --sensmode {params.sensmode} \
                --go_evidence all \
                --dbmem \
                --tax_scope {params.tax_scope} \
                --temp_dir $tmpdir \
                --override \
                --cpu {threads} \
                --output_dir {params.outdir}.temp \
                --output 'proteins'
        mv {params.outdir}.temp {params.outdir}
        # TODO  # Test on 100035 because it has very few genes
        """


rule dbCAN_annotate_translated_orfs:
    output:
        dir=directory("{stem}.dbcan.d"),
    input:
        fasta="{stem}.tran.fa",
        db="ref/dbcan",
    conda:
        "conda/dbcan.yaml"
    threads: 4
    resources:
        walltime_hr=24,
        mem_mb=20_000,
        pmem=20_000 // 4,
    shell:
        """
        run_dbcan \
                {input.fasta} \
                protein \
                --db_dir {input.db} \
                --tools hmmer diamond \
                --tf_cpu {threads} --stp_cpu {threads} --dia_cpu {threads} --hmm_cpu {threads} --dbcan_thread {threads} \
                --out_dir {output.dir}
        """


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
        "data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-{blastn_or_p}.bitscore_ratio-c{centroid}.tsv",
    input:
        script="scripts/calculate_bitscore_ratio_for_gene_matching.py",
        orf_x_orf="data/species/sp-{species}/genome/{stemB}.{stemB}-{blastn_or_p}.tsv",
        orf_x_midas="data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-{blastn_or_p}.tsv",
        midasdb=ancient("ref/midasdb_uhgg"),
    wildcard_constraints:
        blastn_or_p="blastn|blastp",
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


rule assign_matching_genes_based_on_bitscore_ratio:
    output:
        "data/species/sp-{species}/genome/{stemB}.gene_matching-c{centroid}-t{thresh}.tsv",
    input:
        ratio="data/species/sp-{species}/genome/{stemB}.bitscore_ratio-c{centroid}.tsv",
    params:
        thresh=lambda w: int(w.thresh) / 100,
    shell:
        """
        awk -v OFS='\t' -v thresh='{params.thresh}' '(NR > 1) && ($3 > thresh) {{print $1,$2}}' {input} > {output}
        """


# NOTE: So-as to re-use the "gene matching" instrument designed for
# blastn-based accuracy assessment, here we "label" each gene with arbitrary,
# ascending integers.
rule assign_matching_genes_based_on_tile_depth:
    output:
        "{stemA}/species/sp-{species}/genome/{strain}.tiles-{tile_params}.gene{centroidA}-{params}-agg{centroidB}.gene_matching-t{thresh}.tsv",
    input:
        script="scripts/identify_strain_genes_from_tiling_depth.py",
        depth="{stemA}/species/sp-{species}/ALL_STRAINS.tiles-{tile_params}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    params:
        thresh=lambda w: int(w.thresh),
        strain_id=lambda w: f"sp-{w.species}/genome/{w.strain}",
    shell:
        """
        {input.script} {input.depth} {params.strain_id} {params.thresh} {output}
        """


use rule run_bowtie_multispecies_pangenome_v22 as run_bowtie_multispecies_pangenome_on_reference_genome_tiles_v22 with:
    output:
        "data/group/{group}/species/sp-{species}/genome/{genome}.tiles-{tile_params}.pangenomes{centroid}-v22.{bam_or_cram}",
    input:
        db="data/group/{group}/r.proc.pangenomes{centroid}.bt2.d/centroids.bt2db",
        r1="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",
        r2="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",


use rule run_bowtie_species_pangenome_v22 as run_bowtie_species_pangenome_on_reference_genome_tiles_v22 with:
    output:
        "data/group/{group}/species/sp-{species}/genome/{genome}.tiles-{tile_params}.pangenome{centroid}-v22.{bam_or_cram}",
    input:
        db="data/species/sp-{species}/pangenome{centroid}.bt2.d/centroids.bt2db",
        r1="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",
        r2="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",


# rule run_bowtie_multispecies_pangenome_on_reference_genome_tiles_v10:
#     output:
#         "data/group/{group}/species/sp-{species}/genome/{stem}.pangenomes{centroid}-v10.bam",
#     input:
#         bt2_dir="data/group/{group}/r.proc.pangenomes{centroid}.bt2.d",  # Assumes `proc` infix.
#         tiles="data/species/sp-{species}/genome/{stem}.fq.gz",
#     wildcard_constraints:
#         centroid="99|95|90|85|80|75"
#     params:
#         extra_flags="",
#         seed=0,
#     conda:
#         "conda/midas.yaml"
#     threads: 24
#     resources:
#         walltime_hr=24,
#         mem_mb=100_000,
#         pmem=100_000 // 24,
#     shell:
#         """
#         bowtie2 --no-unal --local --very-sensitive-local \
#             -x {input.bt2_dir}/centroids \
#             --threads {threads} --mm -q \
#             -U {input.tiles} \
#             --seed {params.seed} \
#             {params.extra_flags} \
#             | samtools view --threads 2 -u \
#             | samtools sort --threads {threads} -u -T {output} \
#             | samtools view --threads 2 -O BAM --reference {input.bt2_dir}/centroids.fn -t {input.bt2_dir}/centroids.fn.fai -o {output}
#         """
#
#
# use rule run_bowtie_multispecies_pangenome_on_reference_genome_tiles_v10_cram as run_bowtie_multispecies_pangenome_on_reference_genome_tiles_v13_cram with:
#     output:
#         "data/group/{group}/species/sp-{species}/genome/{stem}.pangenomes{centroid}-v13.cram",
#     input:
#         bt2_dir="data/group/{group}/r.proc.pangenomes{centroid}.bt2.d",  # Assumes `proc` infix.
#         tiles="data/species/sp-{species}/genome/{stem}.fq.gz",
#     params:
#         extra_flags="--all",
#         seed=0,
#     threads: 24


# FIXME: This is also a hub rule, and extends another hub rule. Consider
# commenting it out.
# NOTE: The "samples" here are actually reference strains, not samples.
# Then the sample_pattern+sample_list trick is doing something extra tricky
# where the species and strain are combined together.
use rule build_new_pangenome_profiling_db as build_new_pangenome_profiling_db_on_reference_genome_tiles with:
    output:
        protected("data/group/{group}/ALL_STRAINS.{stem}.pangenome{bowtie_params}.db"),
    input:
        samples=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/genome/{strain}.{w.stem}.pangenome{w.bowtie_params}.gene_mapping_tally.tsv.lz4"
            for strain, species in config["genome"].species_id.items()
        ],
        genes="data/group/{group}/r.proc.pangenomes.gene_info.tsv.lz4",
        nlength="data/group/{group}/r.proc.pangenomes99.nlength.tsv",
    params:
        sample_pattern="data/group/{group}/species/$sample.{stem}.pangenome{bowtie_params}.gene_mapping_tally.tsv.lz4",
        sample_list=lambda w: [
            f"sp-{species}/genome/{strain}"
            for strain, species in config["genome"].species_id.items()
        ],


use rule profile_pangenome_mapping_tally_aggregated_by_gene as profile_pangenome_mapping_tally_aggregated_by_gene_for_reference_genome_tiles with:
    output:
        "{stemA}/species/sp-{species}/genome/{genome}.{stem}.pangenome{mapping_params}.gene_mapping_tally.tsv.lz4",
    input:
        "{stemA}/species/sp-{species}/genome/{genome}.{stem}.pangenome{mapping_params}.bam",
