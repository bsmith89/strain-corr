rule start_gtpro_shell:
    container:
        config["container"]["gtpro"]
    shell:
        "bash"


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
    threads: 4
    resources:
        mem_mb=60000,
        pmem=60000 // 4,
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


rule load_gtpro_snp_dict:
    output:
        "ref/gtpro.snp_dict.db",
    input:
        "ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
    shell:
        dd(
            """
        rm -f {output}.tmp
        sqlite3 {output}.tmp <<EOF
        CREATE TABLE snp (
          species TEXT
          , global_pos INT
          , contig TEXT
          , local_pos INT
          , ref_allele VARCHAR(1)
          , alt_allele VARCHAR(1)
          , PRIMARY KEY (species_id, snp_index)
        );
        EOF
        cat {input} \
            | tqdm --unit-scale 1 \
            | sqlite3 -separator '\t' {output}.tmp '.import /dev/stdin snp'
        mv {output}.tmp {output}
        """
        )


# Helper rule that pre-formats paths from library_id to r1 and r2 paths.
rule count_species_lines_from_both_reads_helper:
    output:
        temp("data/{group}.a.r.{stem}.gtpro_species_tally.tsv.args"),
    run:
        with open(output[0], "w") as f:
            for mgen in config["mgen_group"][wildcards.group]:
                print(
                    mgen,
                    f"data/{mgen}.r1.{wildcards.stem}.gtpro_parse.tsv.bz2",
                    f"data/{mgen}.r2.{wildcards.stem}.gtpro_parse.tsv.bz2",
                    sep="\t",
                    file=f,
                )


rule count_species_lines_from_both_reads:
    output:
        "data/{group}.a.r.{stem}.gtpro_species_tally.tsv",
    input:
        script="scripts/tally_gtpro_species_lines.sh",
        r1=lambda w: [
            f"data/{mgen}.r1.{{stem}}.gtpro_parse.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/{mgen}.r2.{{stem}}.gtpro_parse.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
        helper="data/{group}.a.r.{stem}.gtpro_species_tally.tsv.args",
    threads: 24
    shell:
        r"""
        parallel --colsep='\t' --bar -j {threads} \
                bash {input.script} :::: {input.helper} \
            > {output}

        """


rule estimate_all_species_horizontal_coverage:
    output:
        "data/{stem}.gtpro.horizontal_coverage.tsv",
    input:
        script="scripts/estimate_all_species_horizontal_coverage_from_position_tally.py",
        snps="ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
        r="data/{stem}.gtpro_species_tally.tsv",
    shell:
        "{input.script} {input.snps} {input.r} {output}"


# Helper rule that pre-formats paths from library_id *.gtpro_parse.tsv.bz2 files.
rule concatenate_mgen_group_one_read_count_data_from_one_species_helper:
    output:
        temp("data/{group}.a.{stem}.gtpro_combine.tsv.bz2.args"),
    run:
        with open(output[0], "w") as f:
            for mgen in config["mgen_group"][wildcards.group]:
                print(
                    mgen,
                    f"data/{mgen}.{wildcards.stem}.gtpro_parse.tsv.bz2",
                    sep="\t",
                    file=f,
                )


# NOTE: Comment out this rule to speed up DAG evaluation.
# Selects a single species from every file and concatenates.
rule concatenate_mgen_group_one_read_count_data_from_one_species:
    output:
        "data/sp-{species}.{group}.a.{stem}.gtpro_combine.tsv.bz2",
    input:
        script="scripts/select_gtpro_species_lines.sh",
        gtpro=lambda w: [
            f"data/{mgen}.{{stem}}.gtpro_parse.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
        helper="data/{group}.a.{stem}.gtpro_combine.tsv.bz2.args",
    params:
        species=lambda w: w.species,
    threads: 6
    shell:
        dd(
            """
        parallel --colsep='\t' --bar -j {threads} \
                {input.script} {params.species} :::: {input.helper} \
            | bzip2 -c \
            > {output}
        """
        )


rule merge_both_reads_species_count_data:
    output:
        "data/{group_stem}.a.r.{stem}.gtpro_combine.tsv.bz2",
    input:
        script="scripts/sum_merged_gtpro_tables.py",
        r1="data/{group_stem}.a.r1.{stem}.gtpro_combine.tsv.bz2",
        r2="data/{group_stem}.a.r2.{stem}.gtpro_combine.tsv.bz2",
    resources:
        mem_mb=100000,
        pmem=lambda w, threads: 100000 // threads,
    shell:
        """
        {input.script} {input.r1} {input.r2} {output}
        """


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group.
rule estimate_all_species_depth_from_metagenotype:
    output:
        "data/{stem}.species_depth.tsv",
    input:
        script="scripts/estimate_species_depth_from_metagenotype.py",
        mgen=[
            f"data/sp-{species}.{{stem}}.mgen.nc" for species in config["species_list"]
        ],
    params:
        trim=0.05,
        mgen=[
            f"{species}=data/sp-{species}.{{stem}}.mgen.nc"
            for species in config["species_list"]
        ],
    shell:
        "{input.script} {params.trim} {params.mgen} > {output}"


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
