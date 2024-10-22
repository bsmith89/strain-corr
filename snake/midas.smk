rule run_midas_species:
    output:
        directory("data/reads/{mgen}/r.{proc}.species-v20.midas.d"),
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


rule build_midas3_pangenomes_bowtie_index:  # Hub-rule
    output:
        directory("data/hash/{hash}/pangenomes99_v20.bt2.d"),
    input:
        species_list="data/hash/{hash}/species.list",
        midasdb_dir=ancient("ref/midasdb_uhgg_v20_all"),
    conda:
        "conda/midas3.yaml"
    threads: 96
    resources:
        walltime_hr=72,
        mem_mb=480_000,
        pmem=lambda w, threads: 480_000 // threads,
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


rule run_midas_genes_align_only:  # Hub-rule
    output:
        bam="data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes99_v20-v23.bam",
    input:
        species_list="data/hash/{hash}/species.list",
        midasdb_dir="ref/midasdb_uhgg_v20_all",
        bowtie_indexes="data/hash/{hash}/pangenomes99_v20.bt2.d",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    params:
        outdir="data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes99_v20.midas.d",
        outbam="data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes99_v20.midas.d/{mgen}/genes/{mgen}.bam",
    conda:
        "conda/midas3.yaml"
    threads: 4
    resources:
        walltime_hr=72,
        mem_mb=80_000,
        pmem=lambda w, threads: 80_000 // threads,
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
                --alignment_only \
                {params.outdir}
        ln {params.outbam} {output.bam}
        """


# This comment is only needed to get the last rule off the bottom of the file.
