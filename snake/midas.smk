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
        species_list="data/hmp2.a.r.proc.gtpro.horizontal_coverage.filt-20.list",
    conda:
        "conda/midas.yaml"
    threads: 64
    shell:
        """
        midas2 database --download --debug --num_cores {threads} --midasdb_name uhgg --midasdb_dir {input.midasdb} --species_list {input.species_list}
        touch {output}
        """


rule build_midas_one_species_pangenome_index:
    output: directory("ref_temp/midasdb_uhgg_indexes/{species}/pangenomes"),
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
    conda:
        "conda/midas.yaml"
    threads: 1
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
                --num_cores {threads} {params.outdir}
        """


rule build_mgen_group_midas_manifest:
    output:
        "data_temp/{group}.a.r.{stem}.midas_manifest.tsv",
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        outdir=lambda w: f"data_temp/{w.group}.a.r.{w.stem}.midas_output",
    run:
        with open(output[0], "w") as f:
            print("sample_name", "midas_outdir", sep="\t", file=f)
            for mgen in params.mgen_list:
                print(mgen, params.outdir, sep="\t", file=f)


rule merge_midas_genes:
    output:
        directory("data_temp/sp-{species}.{group}.a.r.{stem}.midas_merge/genes"),
    input:
        manifest="data_temp/{group}.a.r.{stem}.midas_manifest.tsv",
        genes=lambda w: [
            f"data_temp/sp-{w.species}.{w.group}.a.r.{w.stem}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    shell:
        """
        midas2 merge_species --samples_list {input.manifest} {output}
        """
