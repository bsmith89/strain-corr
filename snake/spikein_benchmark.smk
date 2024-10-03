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


rule download_gtdb_genome_from_assembly_accession:
    output:
        fasta="raw/genomes/gtdb/{genome}/assembly.fa",
    log:
        genbank="raw/genomes/gtdb/{genome}/genbank.log"
    params:
        ncbi_assembly_name=lambda w: config["genome"].loc[w.genome].ncbi_assembly_name,
    conda:
        "conda/ncbi_datasets.yaml"
    shell:
        """
        genbank_accession=$(esearch -db assembly -query "{params.ncbi_assembly_name}" | esummary | xmllint --xpath "string(//Synonym/RefSeq)" -)
        echo $genbank_accession | tee {log.genbank}
        datasets download genome accession ${{genbank_accession}} --include genome --filename raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip
        unzip -p raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip ncbi_dataset/data/${{genbank_accession}}/${{genbank_accession}}_{params.ncbi_assembly_name}_genomic.fna > {output}
        rm raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip
        """


rule download_sra_run_for_biosample_as_raw_fastq:
    output:
        r1="raw/genomes/gtdb/{genome}/r1.fq.gz",
        r2="raw/genomes/gtdb/{genome}/r2.fq.gz",
    log:
        sra="raw/genomes/gtdb/{genome}/sra.log"
    conda:
        "conda/sra_tools.yaml"
    params:
        biosample=lambda w: config["genome"].loc[w.genome].ncbi_assembly_biosample,
    shell:
        """
        outdir=$(dirname {output.r1})
        srr=$(esearch -db sra -query "{params.biosample}" | esummary | xmllint --xpath 'string(//Runs/Run/@acc)' -)
        echo $srr > {log.sra}
        fasterq-dump --outdir ${{outdir}} ${{srr}}
        gzip -c ${{outdir}}/${{srr}}_1.fastq > {output.r1}
        gzip -c ${{outdir}}/${{srr}}_2.fastq > {output.r2}
        rm ${{outdir}}/${{srr}}_1.fastq ${{outdir}}/${{srr}}_2.fastq
        """
