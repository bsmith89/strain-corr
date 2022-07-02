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
    conda:
        "conda/midas.yaml"
    threads: 4
    resources:
        walltime_hr=24,
        pmem=int(3.2e3 / 1),
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
        midas2 run_genes --sample_name {wildcards.mgen} \
                -1 {input.r1} -2 {input.r2} \
                --midasdb_name uhgg --midasdb_dir {input.midasdb} \
                --species_list $(cat $tmp) --select_threshold=-1 \
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
        directory("data_temp/{group}.a.r.{stem}.midas_merge/genes"),
    input:
        manifest="data_temp/{group}.a.r.{stem}.midas_manifest.tsv",
        genes=lambda w: [
            f"data_temp/{w.group}.a.r.{w.stem}.midas_output/{mgen}/genes"
            for mgen in config["mgen_group"][w.group]
        ],
    shell:
        """
        midas2 merge_species --samples_list {input.manifest} {output}
        """
