use rule start_shell as start_shell_midas with:
    conda:
        "conda/midas.yaml"


rule download_midasdb_uhgg_species:
    output:
        "data/species/sp-{species}/download_species_midasdb_uhgg.flag",
    input:
        midasdb=ancient("ref/midasdb_uhgg"),
    conda:
        "conda/midas.yaml"
    threads: 2
    shell:
        """
        midas2 database --download \
                --debug --num_cores {threads} \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species {wildcards.species}
        touch {output}
        """


localrules:
    download_midasdb_uhgg_species,


rule download_midasdb_species_gene_annotations_all_tsv:
    output:
        directory("ref/midasdb_uhgg_gene_annotations/{species}"),
    params:
        url="s3://microbiome-pollardlab/uhgg_v1/gene_annotations/{species}",
    conda:
        "conda/midas.yaml"
    shell:
        """
        aws s3 cp --quiet --no-sign-request --recursive --exclude "*" --include "*.tsv.lz4" {params.url} {output}
        """


localrules:
    download_midasdb_species_gene_annotations_all_tsv,


rule download_midasdb_species_pangenome_gene_list:
    output:
        "ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
    params:
        url="s3://microbiome-pollardlab/uhgg_v1/pangenomes/{species}/gene_info.txt.lz4",
    conda:
        "conda/midas.yaml"
    shell:
        """
        aws s3 cp --quiet --no-sign-request {params.url} {output}
        """


localrules:
    download_midasdb_species_pangenome_gene_list,


rule build_midas_one_species_pangenome_index:
    output:
        directory("ref/midasdb_uhgg_indexes/{species}/pangenomes"),
    input:
        midasdb=ancient("ref/midasdb_uhgg"),
        species_downloaded_flag="data/species/sp-{species}/download_species_midasdb_uhgg.flag",
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


rule build_midas_multi_species_pangenome_index:
    output:
        directory("data/group/{group}/r.{proc}.pangenomes"),
    input:
        midasdb=ancient("ref/midasdb_uhgg"),
        species_downloaded_flags=lambda w: [
            f"data/species/sp-{species}/download_species_midasdb_uhgg.flag"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
        checkpoint_species="data/group/{group}/r.{proc}.gtpro.horizontal_coverage.tsv",
    params:
        species_list=lambda w: ",".join(
            checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ),
    conda:
        "conda/midas.yaml"
    threads: 100
    shell:
        """
        midas2 build_bowtie2db \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species_list {params.species_list} --select_threshold=-1 \
                --bt2_indexes_name pangenomes --bt2_indexes_dir {output} \
                --num_cores {threads}
        """


rule run_midas_genes_one_species:
    output:
        directory(
            "data/group/{group}/species/sp-{species}/r.{stem}.midas_output/{mgen}/genes"
        ),
    input:
        midasdb=ancient("ref/midasdb_uhgg"),
        bt2_dir="ref/midasdb_uhgg_indexes/{species}/pangenomes",
        r1="data/reads/{mgen}/r1.{stem}.fq.gz",
        r2="data/reads/{mgen}/r2.{stem}.fq.gz",
    params:
        outdir=lambda w: f"data/group/{w.group}/species/sp-{w.species}/r.{w.stem}.midas_output",
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


rule run_midas_genes_multi_species:
    output:
        directory("data/group/{group}/r.{proc}.midas_output/{mgen}/genes"),
    input:
        midasdb=ancient("ref/midasdb_uhgg"),
        bt2_dir="data/group/{group}/r.{proc}.pangenomes",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    params:
        outdir="data/group/{group}/r.{proc}.midas_output",
        min_reads=0,
        min_mapq=0,
    conda:
        "conda/midas.yaml"
    threads: 5
    resources:
        walltime_hr=24,
        mem_mb=100_000,
        pmem=20_000,
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
        "data/group/{group}/{stem}.midas_manifest.tsv",
    wildcard_constraints:
        stemA=endswith_period_or_slash_wc,
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        outdir=lambda w: f"data/group/{w.group}/{w.stem}.midas_output",
    run:
        with open(output[0], "w") as f:
            print("sample_name", "midas_outdir", sep="\t", file=f)
            for mgen in params.mgen_list:
                print(mgen, params.outdir, sep="\t", file=f)


rule merge_midas_genes_one_species:
    output:
        directory("data/group/{group}/species/sp-{species}/r.{stem}.midas_merge/genes"),
    input:
        manifest="data/group/{group}/species/sp-{species}/r.{stem}.midas_manifest.tsv",
        genes=lambda w: [
            f"data/group/{w.group}/species/sp-{w.species}/r.{w.stem}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        outdir="data/group/{group}/species/sp-{species}/r.{stem}.midas_merge",
        midasdb="ref/midasdb_uhgg",
    conda:
        "conda/midas.yaml"
    threads: 24
    resources:
        walltime_hr=48,
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


rule merge_midas_genes_from_multi_species:
    output:
        directory(
            "data/group/{group}/species/sp-{species}/r.{proc}.midas_merge_from_multi/genes"
        ),
    input:
        manifest="data/group/{group}/r.{proc}.midas_manifest.tsv",
        genes=lambda w: [
            f"data/group/{w.group}/r.{w.proc}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        outdir="data/group/{group}/species/sp-{species}/r.{proc}.midas_merge_from_multi",
        midasdb="ref/midasdb_uhgg",
    conda:
        "conda/midas.yaml"
    threads: 4
    resources:
        mem_mb=100_000,
        pmem=100_000 // 4,
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


