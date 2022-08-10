rule start_shell_midas:
    conda:
        "conda/midas.yaml"
    shell:
        "bash"


rule download_midasdb_uhgg_all_group_species:
    output:
        "data_temp/{group}.a.{stem}.gtpro.download_species_midasdb_uhgg.flag",
    input:
        midasdb=ancient("ref_temp/midasdb_uhgg"),
        species_list="data_temp/hmp2.a.r.proc.gtpro.horizontal_coverage.filt-h20-n1.list",
    conda:
        "conda/midas.yaml"
    threads: 64
    shell:
        """
        midas2 database --download \
                --debug --num_cores {threads} \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species_list {input.species_list}
        touch {output}
        """


rule download_midasdb_species_gene_annotations_all_tsv:
    output:
        "data_temp/sp-{species}.download_uhgg_gene_annotations_all_tsv.flag",
    params:
        url="s3://microbiome-pollardlab/uhgg_v1/gene_annotations/{species}",
        outdir="ref_temp/midasdb_uhgg/gene_annotations/{species}",
    conda:
        "conda/midas.yaml"
    shell:
        """
        aws s3 cp --no-sign-request --recursive --exclude "*" --include "*.tsv.lz4" {params.url} {params.outdir}
        touch {output}
        """


rule download_midasdb_species_pangenome_gene_list:
    output:
        "ref_temp/midasdb_uhgg/pangenomes/{species}/gene_info.txt.lz4",
    params:
        url="s3://microbiome-pollardlab/uhgg_v1/pangenomes/{species}/gene_info.txt.lz4",
    conda:
        "conda/midas.yaml"
    shell:
        """
        aws s3 cp --no-sign-request {params.url} {output}
        """


rule build_midas_one_species_pangenome_index:
    output:
        directory("ref_temp/midasdb_uhgg_indexes/{species}/pangenomes"),
    input:
        midasdb=ancient("ref_temp/midasdb_uhgg"),
        # # FIXME: This means that the species must be in the list for {group}.
        # # It would be better if we could ensure that each individual species
        # # had been downloaded.
        # midasdb_downloaded_flag="data_temp/{group}.a.{stem}.gtpro.download_species_midasdb_uhgg.flag",
    conda:
        "conda/midas.yaml"
    threads: 12
    shell:
        """
        midas2 build_bowtie2db \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species_list {wildcards.species} --select_threshold=-1 \
                --bt2_indexes_name pangenomes --bt2_indexes_dir {output} \
                --num_cores {threads}
        """


rule run_midas_genes:
    output:
        directory("data_temp/{group}.a.r.{stem}.midas_output/{mgen}/genes"),
    input:
        midasdb=ancient("ref_temp/midasdb_uhgg"),
        midasdb_downloaded_flag="data_temp/{group}.a.{stem}.gtpro.download_species_midasdb_uhgg.flag",
        r1="data/{mgen}.r1.{stem}.fq.gz",
        r2="data/{mgen}.r2.{stem}.fq.gz",
        species_coverage="data/{group}.a.r.proc.gtpro.horizontal_coverage.tsv",
    params:
        outdir=lambda w: f"data_temp/{w.group}.a.r.{w.stem}.midas_output",
        thresh=0.2,
        min_reads=0,
        min_mapq=0,
    conda:
        "conda/midas.yaml"
    threads: 8
    resources:
        walltime_hr=48,
    shell:
        """
        tmp=$(mktemp)
        awk \
                -v mgen={wildcards.mgen} \
                -v thresh={params.thresh} \
                '$1==mgen && $3 >= thresh {{print $2}}' \
                < {input.species_coverage} \
            | tr '\\n' ',' | sed 's:,$::'\
            > $tmp
        [ -s "$tmp" ] || echo 102506 > $tmp
        cat $tmp
        midas2 run_genes --sample_name {wildcards.mgen} \
                -1 {input.r1} -2 {input.r2} \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species_list $(cat $tmp) --select_threshold=-1 \
                --read_depth {params.min_reads} --aln_mapq {params.min_mapq} \
                --num_cores {threads} {params.outdir}
        """


rule run_midas_genes_one_species:
    output:
        directory("data_temp/sp-{species}.{group}.a.r.{stem}.midas_output/{mgen}/genes"),
    input:
        midasdb=ancient("ref_temp/midasdb_uhgg"),
        bt2_dir="ref_temp/midasdb_uhgg_indexes/{species}/pangenomes",
        r1="data/{mgen}.r1.{stem}.fq.gz",
        r2="data/{mgen}.r2.{stem}.fq.gz",
    params:
        outdir=lambda w: f"data_temp/sp-{w.species}.{w.group}.a.r.{w.stem}.midas_output",
        min_reads=0,
        min_mapq=0,
    conda:
        "conda/midas.yaml"
    threads: 4
    resources:
        walltime_hr=24,
        mem_mb=2_000,
        pmem=2_000,
    shell:
        """
        midas2 run_genes --sample_name {wildcards.mgen} \
                -1 {input.r1} -2 {input.r2} \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --prebuilt_bowtie2_indexes {input.bt2_dir}/pangenomes --prebuilt_bowtie2_species {input.bt2_dir}/pangenomes.species \
                --select_threshold=-1 \
                --read_depth {params.min_reads} --aln_mapq {params.min_mapq} \
                --num_cores {threads} {params.outdir}
        """


rule build_mgen_group_midas_manifest:
    output:
        "{stemA}{group}.a.r.{stemB}.midas_manifest.tsv",
    input:
        genes=lambda w: [
            f"{w.stemA}{w.group}.a.r.{w.stemB}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    wildcard_constraints:
        stemA=endswith_period_or_slash_wc,
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        outdir=lambda w: f"{w.stemA}{w.group}.a.r.{w.stemB}.midas_output",
    run:
        with open(output[0], "w") as f:
            print("sample_name", "midas_outdir", sep="\t", file=f)
            for mgen in params.mgen_list:
                print(mgen, params.outdir, sep="\t", file=f)


rule merge_midas_genes_one_species:
    output:
        directory("data_temp/sp-{species}.{group}.a.r.{stem}.midas_merge/genes"),
    input:
        manifest="data_temp/sp-{species}.{group}.a.r.{stem}.midas_manifest.tsv",
        genes=lambda w: [
            f"data_temp/sp-{w.species}.{w.group}.a.r.{w.stem}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        outdir="data_temp/sp-{species}.{group}.a.r.{stem}.midas_merge",
        midasdb="ref_temp/midasdb_uhgg",
    conda:
        "conda/midas.yaml"
    threads: 24
    shell:
        """
        midas2 merge_genes \
                --num_cores {threads} \
                --samples_list {input.manifest} \
                --species_list {wildcards.species} \
                --midasdb_name uhgg \
                --midasdb_dir {params.midasdb} \
                --genome_depth 1e-3 \
                --cluster_pid 99 \
                {params.outdir}
        """


rule merge_midas_genes:
    output:
        directory("data_temp/{group}.a.r.{stem}.midas_merge/genes"),
    input:
        manifest="data_temp/{group}.a.r.{stem}.midas_manifest.tsv",
        genes=lambda w: [
            f"data_temp/{w.group}.a.r.{w.stem}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        outdir="data_temp/{group}.a.r.{stem}.midas_merge",
        midasdb="ref_temp/midasdb_uhgg",
        species=lambda w: ','.join(
            checkpoint_select_species_with_greater_max_coverage_gtpro(
                group=w.group,
                stem=w.stem,
                cvrg_thresh=0.2,
                num_samples=1,
                require_in_species_group=True,
            )
        )
    conda:
        "conda/midas.yaml"
    threads: 24
    shell:
        """
        midas2 merge_genes \
                --num_cores {threads} \
                --samples_list {input.manifest} \
                --species_list {params.species} \
                --midasdb_name uhgg \
                --midasdb_dir {params.midasdb} \
                --genome_depth 1e-3 \
                --cluster_pid 99 \
                {params.outdir}
        """


# TODO: Make this decompress/move the species file from the traditional MIDAS
# workflow (not the one-species-only workflow).
rule unzip_and_relocate_midas_merge_genes_output:
    output:
        "data_temp/sp-{species}.{stem}.midas_genes.{out_type}.tsv",
    input:
        "data_temp/sp-{species}.{stem}.midas_merge/genes/{species}/{species}.genes_{out_type}.tsv.lz4",
    shell:
        "lz4cat {input} > {output}"
