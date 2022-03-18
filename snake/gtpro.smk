rule download_gtpro_reference_core_snps:
    output:
        "raw/gtpro_refs/variation_in_species/{species_id}/core_snps.vcf.gz",
    params:
        url=lambda w: f"https://fileshare.czbiohub.org/s/waXQzQ9PRZPwTdk/download?path=%2Fvariation_in_species%2F{w.species_id}&files=core_snps.vcf.gz",
    container:
        None
    shell:
        curl_recipe


rule run_gtpro:
    output:
        "{stem}.gtpro_raw.gz",
    input:
        "{stem}.fq.gz",
    params:
        db_l=32,
        db_m=36,
        db_name="ref/gtpro/20190723_881species",
    threads: 8
    resources:
        mem_mb=60000,
        pmem=60000 // 8,
        walltime_hr=4,
    container:
        config["container"]["default"]  # FIXME: Why does gtpro not work here?? File writing is just fine in the other rule...
    shell:
        dd(
            """
        GT_Pro genotype -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_name} -f -o {output} {input}
        mv {output}.tsv.gz {output}
        """
        )


# NOTE: Comment-out this rule after files have been completed to
# save DAG processing time.
rule gtpro_finish_processing_reads:
    output:
        "{stem}.gtpro_parse.tsv.bz2",
    input:
        gtpro="{stem}.gtpro_raw.gz",
    params:
        db="ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
        script="/src/gt-pro/scripts/gtp_parse.py",
    container:
        config["container"]["gtpro"]
    shell:
        dd(
            """
        python3 {params.script} --dict {params.db} --in <(zcat {input.gtpro}) \
                | bzip2 -c \
            > {output}
        """
        )


# NOTE: Comment out this rule to speed up DAG evaluation.
# Selects a single species from every file and concatenates.
rule concatenate_mgen_group_one_read_count_data_from_one_species:
    output:
        "data/{group}.a.{stem}.sp-{species}.gtpro_combine.tsv.bz2",
    input:
        script="scripts/select_gtpro_species_lines.sh",
        gtpro=lambda w: [
            f"data/{mgen}.{{stem}}.gtpro_parse.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        species=lambda w: w.species,
        args=lambda w: "\n".join(
            [
                f"{mgen}\tdata/{mgen}.{w.stem}.gtpro_parse.tsv.bz2"
                for mgen in config["mgen_group"][w.group]
            ]
        ),
    threads: 6
    shell:
        dd(
            """
        tmp=$(mktemp)
        cat >$tmp <<EOF
        {params.args}
        EOF
        parallel --colsep='\t' --bar -j {threads} \
                {input.script} {params.species} :::: $tmp \
            | bzip2 -c \
            > {output}
        """
        )


rule merge_both_reads_species_count_data:
    output:
        "data/{group}.a.{stem}.gtpro.sp-{species}.metagenotype.nc",
    input:
        script="scripts/merge_both_gtpro_reads_to_netcdf.py",
        r1="data/{group}.a.r1.{stem}.sp-{species}.gtpro_combine.tsv.bz2",
        r2="data/{group}.a.r2.{stem}.sp-{species}.gtpro_combine.tsv.bz2",
    threads: 4
    resources:
        mem_mb=100000,
        pmem=lambda w, threads: 100000 // threads,
    shell:
        "{input.script} {input.r1} {input.r2} {output}"


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group.
rule estimate_all_species_coverage_from_metagenotype:
    output:
        touch("data/{stem}.species_cvrg.tsv"),
    input:
        script="scripts/estimate_species_coverage_from_metagenotype.py",
        mgt=[
            f"data/{{stem}}.sp-{species}.metagenotype.nc"
            for species in config["species_list"]
        ],
    params:
        trim=0.05,
        mgt=[
            f"{species}=data/{{stem}}.sp-{species}.metagenotype.nc"
            for species in config["species_list"]
        ],
    shell:
        "{input.script} {params.trim} {params.mgt} > {output}"


rule gather_mgen_group_for_all_species:
    output:
        touch("data/{group}.a.{stemA}.ALL_SPECIES.{stemB}.flag"),
    input:
        [
            f"data/{{group}}.a.{{stemA}}.sp-{species}.{{stemB}}"
            for species in config["species_list"]
        ],
    shell:
        "touch {output}"
