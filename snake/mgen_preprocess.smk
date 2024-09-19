# {{{1 Download and organize reference data


rule alias_raw_read_unsafe_r1:
    output:
        "data/reads/{mgen}/r1.fq.gz",
    input:
        lambda w: config["mgen"]["r1_path"][w.mgen],
    shell:
        alias_recipe


localrules:
    alias_raw_read_unsafe_r1,


rule alias_raw_read_unsafe_r2:
    output:
        "data/reads/{mgen}/r2.fq.gz",
    input:
        lambda w: config["mgen"]["r2_path"][w.mgen],
    shell:
        alias_recipe


localrules:
    alias_raw_read_unsafe_r2,


# {{{1 Process data
# {{{2 Metagenomic reads

# Useful for when no additional processing is necessary.
rule dummy_operation_on_reads:
    output:
        temp("{stem}.noop.fq.gz"),
    input:
        "{stem}.fq.gz",
    shell:
        alias_recipe


localrules:
    dummy_operation_on_reads,


rule alias_cleaned_reads:
    output:
        "data/reads/{mgen}/{r}.proc.fq.gz",
    input:
        lambda w: f"data/reads/{w.mgen}/{w.r}."
        + config["mgen"]["preprocessing"][w.mgen]
        + ".fq.gz",
    shell:
        alias_recipe


localrules:
    alias_cleaned_reads,


# {{{1 Checkpoint rules
# NOTE: These may be useful for other parts of the workflow besides
# just pre-processing.


rule gather_all_mgen_read_pairs_from_mgen_group:
    output:
        touch("data/group/{group}/r.{stem}.ALL_MGEN_PAIRS.flag"),
    input:
        r1=lambda w: [
            f"data/reads/{mgen}/r1.{{stem}}" for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"data/reads/{mgen}/r2.{{stem}}" for mgen in config["mgen_group"][w.group]
        ],
    shell:
        "touch {output}"


localrules:
    gather_all_mgen_read_pairs_from_mgen_group,


rule gather_all_mgen_from_mgen_group:
    output:
        touch("data/group/{group}/r.{stem}.ALL_MGEN.flag"),
    input:
        lambda w: [
            f"data/reads/{mgen}/r.{{stem}}" for mgen in config["mgen_group"][w.group]
        ],
    shell:
        "touch {output}"


localrules:
    gather_all_mgen_from_mgen_group,
