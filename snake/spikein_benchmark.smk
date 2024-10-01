rule download_gtdb_metadata:
    output:
        "raw/gtdb/bac120_metadata_r220.tsv.gz",
    params:
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz",
    shell:
        curl_recipe


rule unzip_gtdb_metadata:
    output:
        "ref/gtdb/metadata.tsv",
    input:
        "raw/gtdb/bac120_metadata_r220.tsv.gz",
    shell:
        "gzip -dc {input} > {output}"


rule subset_ecoli_metadata:
    output:
        "ref/gtdb/species/102506/metadata.tsv",
    input:
        "ref/gtdb/metadata.tsv",
    params:
        taxonomy_string="d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli",
    shell:
        """
        awk 'NR==1 || $0~/{params.taxonomy_string}/' {input} > {output}
        """


rule download_gtdb_genome:
    output:
        "raw/genomes/gtdb/{genome}/assembly.fa",
    params:
        url=lambda w: config["genome"].loc[w.genome].assembly_download_url,
    shell:
        curl_unzip_recipe


rule download_sra_raw_fastq:
    output:
        r1="raw/sra/{srr}_1.fastq",
        r2="raw/sra/{srr}_2.fastq",
    conda:
        "conda/sra_tools.yaml"
    shell:
        """
        fasterq-dump --outdir raw/sra {wildcards.srr}
        """


rule gzip_sra_fastq:
    output:
        "raw/sra/{stem}.fq.gz",
    input:
        "raw/sra/{stem}.fastq",
    shell:
        "gzip -c {input} > {output}"


rule link_sra_reads_to_genome_dir:
    output:
        "raw/genomes/gtdb/{genome}/r{read}.fq.gz",
    input:
        lambda w: "raw/sra/{srr}_{{read}}.fq.gz".format(
            srr=config["genome"].loc[w.genome].sra_run_accession
        ),
    shell:
        alias_recipe
