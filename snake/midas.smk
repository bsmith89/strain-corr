rule run_midas_species:
    output:
        directory("data/reads/{mgen}/r.{proc}.species-v3.midas.d"),
    input:
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
        midasdb_dir="ref/midasdb_uhgg_v20_all",
    conda:
        "conda/midas3.yaml"
    threads: 4
    resources:
        walltime_hr=3,
        mem_mb=1_000,
        pmem=1_000 // 4,
    shell:
        """
        midas2 run_species \
                --midasdb_dir {input.midasdb_dir} \
                --midasdb_name newdb \
                -1 {input.r1} -2 {input.r2} \
                --num_cores {threads} \
                --sample_name {wildcards.mgen} \
                {output}
        """


rule write_species_list:
    output:
        "data/hash/{hash}/species.list",
    params:
        species_list=lambda w: "\n".join(config["hash_to_species_set"][w.hash]),
    shell:
        dd(
            """
        cat > {output} <<EOF
        {params.species_list}
        EOF
        """
        )


rule build_midas3_pangenomes_bowtie_index:
    output:
        directory("data/hash/{hash}/pangenomes99_v3.bt2.d"),
    input:
        species_list="data/hash/{hash}/species.list",
        midasdb_dir="ref/midasdb_uhgg_v20_all",
    conda:
        "conda/midas3.yaml"
    threads: 64
    resources:
        walltime_hr=72,
        mem_mb=640_000,
    shell:
        """
        midas2 build_bowtie2db \
                --bt2_indexes_dir {output} \
                --bt2_indexes_name pangenomes \
                --midasdb_name newdb \
                --midasdb_dir {input.midasdb_dir} \
                --select_threshold=-1 \
                --species_list {input.species_list} \
                --num_cores {threads} \
                --prune_centroids \
                --remove_singleton
        """


rule run_midas_genes:
    output:
        directory("data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes99-v3.midas.d"),
    input:
        species_list="data/hash/{hash}/species.list",
        midasdb_dir="ref/midasdb_uhgg_v20_all",
        bowtie_indexes="data/hash/{hash}/pangenomes99_v3.bt2.d",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    conda:
        "conda/midas3.yaml"
    threads: 8
    resources:
        walltime_hr=48,
        mem_mb=160_000,
        pmem=160_000 // 8,
    shell:
        """
        midas2 run_genes --num_cores {threads} \
                -1 {input.r1} -2 {input.r2} \
                --sample_name {wildcards.mgen} \
                --midasdb_name newdb \
                --midasdb_dir {input.midasdb_dir} \
                --prebuilt_bowtie2_indexes {input.bowtie_indexes}/pangenomes \
                --prebuilt_bowtie2_species {input.species_list} \
                --select_threshold=-1 \
                --aln_speed sensitive \
                --aln_extra_flags '--mm --ignore-quals' \
                --total_depth 0 \
                --cluster_level 75 \
                {output}
        """


rule load_one_species_pangenome3_depth_into_netcdf:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene99_v30-v23-agg75.depth2.nc",
    input:
        script="scripts/merge_midas_pangenomes_depth.py",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenomes99-v3.midas.d".format(
                w=w, mgen=mgen, _hash=config["species_group_to_hash"][w.group]
            )
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        args=lambda w: [
            "{mgen}=data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenomes99-v3.midas.d/{mgen}/genes/{w.species}.genes.tsv.lz4".format(
                w=w, mgen=mgen, _hash=config["species_group_to_hash"][w.group]
            )
            for mgen in config["mgen_group"][w.group]
        ],
    conda:
        "conda/toolz.yaml"
    threads: 1
    resources:
        walltime_hr=24,
        mem_mb=20_000,
        pmem=20_000 // 1,
    shell:
        """
        {input.script} {output} {params.args}
        """


# This comment is only needed to get the last rule off the bottom of the file.
