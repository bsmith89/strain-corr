use rule start_shell as start_shell_gtpro with:
    container:
        config["container"]["gtpro"]


rule download_gtpro_reference_core_snps:
    output:
        "raw/gtpro_refs/variation_in_species/{species_id}/core_snps.vcf.gz",
    params:
        url=lambda w: f"https://fileshare.czbiohub.org/s/waXQzQ9PRZPwTdk/download?path=%2Fvariation_in_species%2F{w.species_id}&files=core_snps.vcf.gz",
    shell:
        curl_recipe


localrules:
    download_gtpro_reference_core_snps,


rule run_gtpro:
    output:
        temp("{stem}.gtpro_raw.gz"),
    input:
        r="{stem}.fq.gz",
        db="ref/gtpro",
    params:
        db_l=32,
        db_m=36,
        db_name="ref/gtpro/20190723_881species",
    threads: 8
    resources:
        mem_mb=10000,
        pmem=10000 // 8,
        walltime_hr=8,
    container:
        config["container"]["gtpro"]
    shell:
        """
        GT_Pro genotype -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_name} -o {output}.temp {input.r}
        mv {output}.temp.tsv.gz {output}
        """


# # This rule seems to work where the other fails.
# # Test cases:
# # data/species/sp-103694/genome/midasdb_v20/GUT_GENOME037857.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-100084/genome/midasdb_v20/GUT_GENOME036376.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-101359/genome/midasdb_v20/GUT_GENOME251829.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-103899/genome/midasdb_v20/GUT_GENOME254885.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-100189/genome/midasdb_v20/GUT_GENOME243557.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-101338/genome/midasdb_v20/GUT_GENOME058802.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-102492/genome/midasdb_v20/GUT_GENOME066406.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-100060/genome/midasdb_v20/GUT_GENOME122534.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-101346/genome/midasdb_v20/GUT_GENOME054037.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-102512/genome/midasdb_v20/GUT_GENOME193443.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-102545/genome/midasdb_v20/GUT_GENOME125056.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# # data/species/sp-102478/genome/midasdb_v20/GUT_GENOME038804.norm.tiles-l500-o31.gtpro_parse.tsv.bz2
# rule run_gtpro_stopgap:
#     output:
#         temp("{stem}.gtpro_raw.gz"),
#     input:
#         r="{stem}.fq.gz",
#         db="ref/gtpro",
#     params:
#         db_l=32,
#         db_m=36,
#         db_name="ref/gtpro/20190723_881species",
#     threads: 8
#     resources:
#         mem_mb=10000,
#         pmem=10000 // 8,
#         walltime_hr=8,
#     container:
#         config["container"]["gtpro"]
#     shell:
#         """
#         gzip -cd {input.r} | GT_Pro genotype -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_name} | gzip -c > {output}.temp
#         mv {output}.temp {output}
#         """


rule load_gtpro_snp_dict:
    output:
        "ref/gtpro.snp_dict.db",
    input:
        "ref/gtpro.snp_dict.tsv",
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
          , PRIMARY KEY (species, global_pos)
        );
        EOF
        cat {input} \
            | tqdm --unit-scale 1 \
            | sqlite3 -separator '\t' {output}.tmp '.import /dev/stdin snp'
        mv {output}.tmp {output}
        """
        )


rule gtpro_finish_processing_reads:
    output:
        "{stem}.gtpro_parse.tsv.bz2",
    input:
        raw="{stem}.gtpro_raw.gz",
        db="ref/gtpro.snp_dict.db",
    shell:
        dd(
            """
        rm -f {output}.tmp
        sqlite3 {output}.tmp <<EOF

        CREATE TABLE _gtpro_tally (
            snv_id TEXT
          , tally INT
        );

        EOF
        zcat {input.raw} \
            | tqdm --unit-scale 1 \
            | sqlite3 -separator '\t' {output}.tmp '.import /dev/stdin _gtpro_tally'
        (
        sqlite3 -header -separator '\t' {output}.tmp <<EOF

        ATTACH DATABASE '{input.db}' AS ref;

        CREATE TEMPORARY VIEW gtpro_tally AS
        SELECT
          snv_id
        , substr(snv_id, 1, 6) AS species
        , substr(snv_id, 7, 1) AS snv_type
        , substr(snv_id, 8) AS global_pos
        , tally
        FROM _gtpro_tally
        ;

        CREATE TEMPORARY VIEW snp_hit AS
        SELECT *
        , CASE snv_type
            WHEN '0' THEN tally
            WHEN '1' THEN 0
        END AS ref_count
        , CASE snv_type
            WHEN '0' THEN 0
            WHEN '1' THEN tally
        END AS alt_count
        FROM gtpro_tally
        JOIN ref.snp USING (species, global_pos)
        ;

        SELECT
            species
            , global_pos
            , contig
            , local_pos
            , ref_allele
            , alt_allele
            , SUM(ref_count) AS ref_count
            , SUM(alt_count) AS alt_count
        FROM snp_hit
        GROUP BY species, global_pos
        ORDER BY CAST(species AS INT), CAST(global_pos AS INT)
        ;

        EOF
        ) | bzip2 -c > {output}.tmp2
        rm {output}.tmp
        mv {output}.tmp2 {output}
        """
        )


rule count_species_lines_from_both_reads_helper:  # Hub-rule
    output:
        "data/group/{group}/r.{proc}.gtpro_species_tally.tsv.args",
    input:
        reads=lambda w: [
            f"data/reads/{mgen}/{r}.{w.proc}.gtpro_parse.tsv.bz2"
            for mgen, r in product(
                config["mgen_group"][w.group],
                ["r1", "r2"],
            )
        ],
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        pattern=lambda w: f"data/reads/{{mgen}}/{{r}}.{w.proc}.gtpro_parse.tsv.bz2",
    run:
        with open(output[0], "w") as f:
            for mgen in params.mgen_list:
                print(
                    mgen,
                    params.pattern.format(mgen=mgen, r="r1"),
                    params.pattern.format(mgen=mgen, r="r2"),
                    sep="\t",
                    file=f,
                )


rule count_species_lines_from_both_reads:
    output:
        "{stem}.gtpro_species_tally.tsv",
    input:
        script="scripts/tally_gtpro_species_lines.sh",
        helper="{stem}.gtpro_species_tally.tsv.args",
    threads: 12
    shell:
        r"""
        parallel --colsep='\t' --bar -j {threads} \
                bash {input.script} :::: {input.helper} \
            > {output}.tmp
        mv {output}.tmp {output}

        """


checkpoint estimate_all_species_horizontal_coverage_with_gtpro:
    output:
        "data/group/{group}/r.{proc}.gtpro.horizontal_coverage.tsv",
    input:
        script="scripts/estimate_all_species_horizontal_coverage_from_position_tally.py",
        snps="ref/gtpro.snp_dict.tsv",
        r="data/group/{group}/r.{proc}.gtpro_species_tally.tsv",
    shell:
        "{input.script} {input.snps} {input.r} {output}"


def checkpoint_estimate_all_species_horizontal_coverage_with_gtpro(path):
    chkpt, w = get_checkpoint_by_path(
        checkpoints.estimate_all_species_horizontal_coverage_with_gtpro, path
    )
    return chkpt, w


# FIXME: This is too convoluted of a checkpointing mechanism. I think I should
# go back to the checkpoint rule being a species list and the checkpoint function
# taking that path, using get_checkpoint_by_path(path) to get the checkpoint reference
# and then simply reading the species list out of that file.
# Parameterization by horizontal coverage and sample count should be included
# in the file path directly.
def checkpoint_select_species(
    path, cvrg_thresh, num_samples, require_in_species_group=False
):
    chkpt, w = checkpoint_estimate_all_species_horizontal_coverage_with_gtpro(path)
    d = (
        pd.read_table(
            chkpt.output[0],
            names=["sample", "species", "max_coverage"],
            dtype={"sample": str, "species": str, "max_coverage": float},
            index_col=["sample", "species"],
        )
        .squeeze()
        .unstack("species")
        .astype(float)
        .fillna(0)
    )
    species_with_sufficient_coverage = idxwhere((d >= cvrg_thresh).sum() >= num_samples)
    if require_in_species_group:
        out = sorted(
            set(species_with_sufficient_coverage)
            & set(config["species_group"][w["group"]])
        )
    else:
        out = sorted(species_with_sufficient_coverage)
    return out


rule concatenate_mgen_group_one_read_count_data_from_one_species_helper:
    output:
        "data/group/{group}/{stem}.gtpro.tsv.bz2.args",
    input:
        reads=lambda w: [
            f"data/reads/{mgen}/{w.stem}.gtpro_parse.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        mgen_list=lambda w: config["mgen_group"][w.group],
        pattern=lambda w: f"data/reads/{{mgen}}/{w.stem}.gtpro_parse.tsv.bz2",
    run:
        with open(output[0], "w") as f:
            for mgen in params.mgen_list:
                print(
                    mgen,
                    params.pattern.format(mgen=mgen),
                    sep="\t",
                    file=f,
                )


rule concatenate_mgen_group_one_read_count_data_from_one_species:  # Hub-rule?
    output:
        "{stemA}/species/sp-{species}/{r12}.{stemB}.gtpro.tsv.bz2",
    input:
        script="scripts/select_gtpro_species_lines.sh",
        helper="{stemA}/{r12}.{stemB}.gtpro.tsv.bz2.args",
    wildcard_constraints:
        r12="r[12]",
    params:
        species=lambda w: w.species,
    threads: 6
    shell:
        dd(
            """
        parallel --colsep='\t' --bar -j {threads} \
                {input.script} {params.species} :::: {input.helper} \
            | bzip2 -c \
            > {output}.tmp
        mv {output}.tmp {output}
        """
        )


rule merge_both_reads_species_count_data:
    output:
        "{stemA}/r.{stemB}.gtpro.tsv.bz2",
    input:
        script="scripts/sum_merged_gtpro_tables.py",
        r1="{stemA}/r1.{stemB}.gtpro.tsv.bz2",
        r2="{stemA}/r2.{stemB}.gtpro.tsv.bz2",
    resources:
        mem_mb=100000,
        pmem=lambda w, threads: 100000 // threads,
        walltime_hr=4,
    shell:
        """
        {input.script} {input.r1} {input.r2} {output}
        """


rule load_metagenotype_from_merged_gtpro:  # Hub-rule
    output:
        "{stem}.gtpro.mgtp.nc",
    input:
        "{stem}.gtpro.tsv.bz2",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts load --gtpro-metagenotype {input} {output}
        """


# FIXME: Change script to output only a mapping sample\tdepth; no header, no second index.
rule estimate_species_depth_from_metagenotype:
    output:
        "{stemA}/species/sp-{species}/r.{stemB}.gtpro.species_depth.tsv",
    input:
        script="scripts/estimate_species_depth_from_metagenotype.py",
        mgen="{stemA}/species/sp-{species}/r.{stemB}.gtpro.mgtp.nc",
    params:
        trim=0.05,
    shell:
        "{input.script} {params.trim} {input.mgen} {output}"


rule estimate_all_species_depths_in_group_gtpro:  # Hub-rule
    output:
        "data/group/{group}/r.{proc}.gtpro.all_species_depth.tsv",
    input:
        species=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/r.{w.proc}.gtpro.species_depth.tsv"
            for species in config["species_group"][w.group]
        ],
    params:
        header="sample	species_id	depth",
        species_list=lambda w: config["species_group"][w.group],
        species_pattern="data/group/{group}/species/sp-$species/r.{proc}.gtpro.species_depth.tsv",
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                file={params.species_pattern}
                echo $file >&2
                awk -v species=$species -v OFS='\t' '{{print $1,species,$2}}' $file
            done
        ) > {output}
        """


rule gather_mgen_group_for_all_species:
    output:
        touch("{stemA}/ALL_SPECIES/{stemB}.flag"),
    input:
        lambda w: [
            f"{w.stemA}/species/sp-{species}/{w.stemB}"
            for species in config["species_group"][w.group]
        ],
    shell:
        "touch {output}"


# TODO: Move this rule to a more obvious location.
rule construct_group_files_for_all_select_species:
    output:
        touch("data/group/{group}/{stem}.SELECT_SPECIES.flag"),
    input:
        lambda w: [
            f"data/group/{w.group}/species/sp-{species}/{w.stem}"
            for species in config["species_group"][w.group]
        ],


rule construct_groupfree_files_for_all_select_species:
    output:
        touch("data/group/{group}/{stem}.SELECT_SPECIES.flag"),
    input:
        lambda w: [
            f"data/species/sp-{species}/{w.stem}"
            for species in config["species_group"][w.group]
        ],


ruleorder: construct_group_files_for_all_select_species > construct_groupfree_files_for_all_select_species
