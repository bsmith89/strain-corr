# {{{1 Download and organize reference data


rule download_illumina_adapters:
    output:
        "raw/ref/illumina_adapters.fn",
    params:
        url="https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa",
    resources:
        network_connections=1,
    shell:
        curl_recipe


rule link_illumina_adapters:
    output:
        "ref/illumina_adapters.fn",
    input:
        "raw/ref/illumina_adapters.fn",
    shell:
        alias_recipe


rule download_GRCh38_index:
    output:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    resources:
        network_connections=1,
    shell:
        curl_recipe


rule unpack_GRCh38_index:
    output:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.1.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.2.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.3.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.4.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.1.bt2",
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.rev.2.bt2",
    input:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz",
    shell:
        "tar -C raw/ref/ -xzf {input}"


rule alias_GRCh38_index_file:
    output:
        "ref/GRCh38.{suffix}",
    input:
        "raw/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.{suffix}",
    shell:
        alias_recipe


ruleorder: alias_GRCh38_index_file > bowtie_index_build


# {{{1 Organize raw data


rule alias_raw_read_r1:
    output:
        "sdata/reads/{mgen}/r1.fq.gz",
    input:
        lambda w: config["mgen"]["r1_path"][w.mgen],
    shell:
        alias_recipe


localrules:
    alias_raw_read_r1,


rule alias_raw_read_r2:
    output:
        "sdata/reads/{mgen}/r2.fq.gz",
    input:
        lambda w: config["mgen"]["r2_path"][w.mgen],
    shell:
        alias_recipe


localrules:
    alias_raw_read_r2,


# {{{1 Process data
# {{{2 Metagenomic reads


rule qc_reads:
    output:
        directory("{stemA}/group/{group}/a.r.{stemB}.fastqc.d"),
    input:
        r1=lambda w: [
            f"{{stemA}}/reads/{mgen}/r1.{{stemB}}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
        r2=lambda w: [
            f"{{stemA}}/reads/{mgen}/r2.{{stemB}}.fq.gz"
            for mgen in config["mgen_group"][w.group]
        ],
    container:
        config["container"]["toolz"]
    threads: config["MAX_THREADS"]
    shell:
        dd(
            """
        mkdir -p {output}
        fastqc -t {threads} -o {output} {input}
        """
        )


# Useful for when no additional processing is necessary.
rule dummy_operation_on_reads:
    output:
        temp("{stem}.noop.fq.gz"),
    input:
        "{stem}.fq.gz",
    shell:
        alias_recipe


rule deduplicate_reads:
    output:
        r1=temp("{stemA}/r1{stemB}dedup.fq.gz"),
        r2=temp("{stemA}/r2{stemB}dedup.fq.gz"),
    input:
        script="scripts/fastuniq_wrapper.sh",
        r1="{stemA}/r1{stemB}fq.gz",
        r2="{stemA}/r2{stemB}fq.gz",
    resources:
        mem_mb=10_000,
        walltime_min=600,
    container:
        config["container"]["toolz"]
    shell:
        "{input.script} {input.r1} {input.r2} {output.r1} {output.r2}"


rule trim_adapters:
    output:
        fq=temp("{stem}.deadapt.fq.gz"),
    input:
        adapters="ref/illumina_adapters.fn",
        fq="{stem}.fq.gz",
    log:
        "log/{stem}.scythe.log",
    container:
        config["container"]["toolz"]
    threads: 2
    resources:
        walltime_hr=2,
    shell:
        dd(
            """
        scythe -a {input.adapters} {input.fq} 2>{log} | gzip -c > {output.fq}
        ! grep -Fxq 'Blank FASTA header or sequence in adapters file.' {log}
        """
        )


rule quality_trim_reads:
    output:
        r1=temp("{stemA}/r1.{stemB}.qtrim.fq.gz"),
        r2=temp("{stemA}/r2.{stemB}.qtrim.fq.gz"),
        r3=temp("{stemA}/r3.{stemB}.qtrim.fq.gz"),
    input:
        r1="{stemA}/r1.{stemB}.fq.gz",
        r2="{stemA}/r2.{stemB}.fq.gz",
    params:
        qual_type="sanger",
        qual_thresh=20,
    container:
        config["container"]["toolz"]
    resources:
        walltime_hr=3,
    shell:
        dd(
            """
        sickle pe -t {params.qual_type} -q {params.qual_thresh} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """
        )


# NOTE: Input from sdata/ output to data/
# NOTE: Had to drop periods on either side of {stem} to account for
# NULL case.
rule filter_out_host:
    output:
        r1="data/reads/{mgen}/r1{stem}hfilt.fq.gz",
        r2="data/reads/{mgen}/r2{stem}hfilt.fq.gz",
    input:
        script="scripts/filter_out_mapping_reads.sh",
        r1="sdata/reads/{mgen}/r1{stem}fq.gz",
        r2="sdata/reads/{mgen}/r2{stem}fq.gz",
        index=[
            "ref/GRCh38.1.bt2",
            "ref/GRCh38.2.bt2",
            "ref/GRCh38.3.bt2",
            "ref/GRCh38.4.bt2",
            "ref/GRCh38.rev.1.bt2",
            "ref/GRCh38.rev.2.bt2",
        ],
    params:
        index="ref/GRCh38",
    container:
        config["container"]["toolz"]
    threads: 8
    shell:
        dd(
            """
        {input.script} {threads} {params.index} {input.r1} {input.r2} {output.r1} {output.r2}
        """
        )


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once all libraries have been processed
rule alias_cleaned_reads:
    output:
        "data/{stem}.proc.fq.gz",
    input:
        f"data/{{stem}}.{config['proc']}.fq.gz",
    shell:
        alias_recipe


localrules:
    alias_cleaned_reads,


# {{{1 Checkpoint rules
# NOTE: These may be useful for other parts of the workflow besides
# just pre-processing.


rule gather_all_mgen_read_pairs_from_mgen_group:
    output:
        touch("data/group/{group}/a.{stem}.ALL_MGEN_PAIRS.flag"),
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
        touch("data/group/{group}/a.{stem}.ALL_MGEN.flag"),
    input:
        lambda w: [f"data/reads/{mgen}/r.{{stem}}" for mgen in config["mgen_group"][w.group]],
    shell:
        "touch {output}"


localrules:
    gather_all_mgen_from_mgen_group,
