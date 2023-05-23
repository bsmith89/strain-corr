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
    threads: 4
    resources:
        mem_mb=60000,
        pmem=60000 // 4,
        walltime_hr=4,
    container:
        config["container"]["gtpro"]
    shell:
        dd(
            """
        GT_Pro genotype -t {threads} -l {params.db_l} -m {params.db_m} -d {params.db_name} -f -o {output} {input.r}
        mv {output}.tsv.gz {output}
        """
        )


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


# NOTE: Comment-out this rule after files have been completed to
# save DAG processing time.
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


rule count_species_lines_from_both_reads:
    output:
        "{stem}.gtpro_species_tally.tsv",
    input:
        script="scripts/tally_gtpro_species_lines.sh",
        helper="{stem}.gtpro_species_tally.tsv.args",
    threads: 24
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


rule list_checkpoint_select_species:
    output:
        "data/group/{group}/r.{proc}.gtpro.horizontal_coverage.select_species.list",
    input:
        "data/group/{group}/r.{proc}.gtpro.horizontal_coverage.tsv",
    params:
        obj=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
    shell:
        "for species in {params.obj}; do echo $species; done > {output}"


# NOTE: Comment out this rule to speed up DAG evaluation.
# Selects a single species from every file and concatenates.
rule concatenate_mgen_group_one_read_count_data_from_one_species:
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


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group/species.
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


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group.
rule estimate_all_species_depths:
    output:
        "data/group/{group}/r.{proc}.gtpro.species_depth.tsv",
    input:
        species=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/r.{w.proc}.gtpro.species_depth.tsv"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        header="sample	species_id	depth",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
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


rule construct_files_for_all_select_species:
    output:
        touch("data/group/{group}/r.{proc}.gtpro.{suffix}.SELECT_SPECIES.flag"),
    input:
        lambda w: [
            f"data/group/{w.group}/species/sp-{species}/r.{w.proc}.gtpro.{w.suffix}"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
