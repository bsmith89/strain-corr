rule download_gtdb_metadata:
    output:
        "raw/gtdb/bac120_metadata_r220.tsv.gz",
    params:
        url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz",
    resources:
        network_connections=1,
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


rule subset_fprausnitzii_e_metadata:
    output:
        "ref/gtdb/species/100195/metadata.tsv",
    input:
        "ref/gtdb/metadata.tsv",
    params:
        taxonomy_string="d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium prausnitzii_E",
    shell:
        """
        awk 'NR==1 || $0~/{params.taxonomy_string}/' {input} > {output}
        """


rule download_gtdb_genome_from_assembly_accession:
    output:
        fasta="raw/genomes/gtdb/{genome}/assembly.fa",
    log:
        genbank="raw/genomes/gtdb/{genome}/genbank.log",
    params:
        ncbi_assembly_name=lambda w: config["genome"].loc[w.genome].ncbi_assembly_name,
    resources:
        network_connections=1,
    conda:
        "conda/ncbi_datasets.yaml"
    shell:
        """
        genbank_accession=$(esearch -db assembly -query "{params.ncbi_assembly_name}" | esummary | xmllint --xpath "string(//Synonym/RefSeq)" -)
        echo $genbank_accession | tee {log.genbank}
        datasets download genome accession ${{genbank_accession}} --include genome --filename raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip
        unzip -Z raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip "ncbi_dataset/data/${{genbank_accession}}/*.fna" > {log.genbank}
        unzip -p raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip "ncbi_dataset/data/${{genbank_accession}}/*.fna" > {output}
        rm raw/genomes/gtdb/{wildcards.genome}/ncbi_dataset.zip
        """


rule download_genbank_genome_from_assembly_accession:
    output:
        fasta="raw/genomes/ncbi/{genome}/assembly.fa",
    log:
        esearch="raw/genomes/ncbi/{genome}/esearch_assembly.log",
    params:
        ncbi_assembly_name=lambda w: config["genome"].loc[w.genome].ncbi_assembly_name,
    resources:
        network_connections=1,
    conda:
        "conda/ncbi_datasets.yaml"
    shell:
        """
        esearch -db assembly -query "{params.ncbi_assembly_name}" | esummary > {log.esearch}
        ftp=$(xmllint --xpath "string(//FtpPath_GenBank)" {log.esearch})
        asm=$(xmllint --xpath "string(//AssemblyName)" {log.esearch})
        genbank=$(xmllint --xpath "string(//Synonym/Genbank)" {log.esearch})
        curl ${{ftp}}/${{genbank}}_${{asm}}_genomic.fna.gz | zcat > {output}
        """


rule simulate_reads_from_ncbi_genome:
    output:
        r1="data/reads/sim_{genome}_len150_seed{seed}_cov{cov}/r1.fq.gz",
        r2="data/reads/sim_{genome}_len150_seed{seed}_cov{cov}/r2.fq.gz",
    input:
        fasta="data/genome/{genome}.fn",
    params:
        outdir="data/reads/sim_{genome}_len150_seed{seed}_cov{cov}",
    conda:
        "conda/art_read_sim.yaml"
    shell:
        """
        art_illumina --rndSeed {wildcards.seed} \
                --seqSys HS25 \
                --len 150 --paired --mflen 450 --sdev 30 \
                --fcov {wildcards.cov} \
                --noALN \
                --in {input.fasta} --out {params.outdir}/
        gzip -c {params.outdir}/1.fq > {output.r1} && rm {params.outdir}/1.fq
        gzip -c {params.outdir}/2.fq > {output.r2} && rm {params.outdir}/2.fq
        """


# # Rules below are meant to more efficiently implement this rule by merging downstream.
# # This way we don't need to re-run all of the read mapping and GT-Pro counting.
# rule spikein_simulated_benchmark_isolate_to_mgen_fastq:
#     output:
#         r1="data/reads/{mgen}_spikein_{genome}_seed{seed}_cov{cov}/r1.proc.fq.gz",
#         r2="data/reads/{mgen}_spikein_{genome}_seed{seed}_cov{cov}/r2.proc.fq.gz"
#     input:
#         r1_mgen="data/reads/{mgen}/r1.proc.fq.gz",
#         r1_sim="data/reads/sim_{genome}_len150_seed{seed}_cov{cov}/r1.proc.fq.gz",
#         r2_mgen="data/reads/{mgen}/r2.proc.fq.gz",
#         r2_sim="data/reads/sim_{genome}_len150_seed{seed}_cov{cov}/r2.proc.fq.gz",
#     shell:
#         """
#         cat {input.r1_mgen} {input.r1_sim} > {output.r1}
#         cat {input.r2_mgen} {input.r2_sim} > {output.r2}
#         """


# rule spikein_simulated_benchmark_isolate_to_mgen_bam:
#     output: "{stemA}/reads/{mgen}_spikein_{genome}_seed{seed}_cov{cov}/{stemB}.bam"
#     input:
#         mgen="{stemA}/reads/{mgen}/{stemB}.bam",
#         sim="{stemA}/reads/sim_{genome}_len150_seed{seed}_cov{cov}/{stemB}.bam",
#     conda:
#         "conda/midas3.yaml"
#     shell:
#         """
#         samtools merge -o {output} {input.mgen} {input.sim}
#         """


# # TODO: Replace spikein_simulated_benchmark_isolate_to_mgen_bam with this for
# # efficiency. BAM mergeing is _expensive_.
rule spikein_simulated_benchmark_isolate_to_mgen_read_tally:
    output:
        "{stemA}/reads/{mgen}_spikein_{genome}_seed{seed}_cov{cov}/{stemB}.gene_mapping_tally.tsv.lz4",
    input:
        mgen="{stemA}/reads/{mgen}/{stemB}.gene_mapping_tally.tsv.lz4",
        sim="{stemA}/reads/sim_{genome}_len150_seed{seed}_cov{cov}/{stemB}.gene_mapping_tally.tsv.lz4",
    conda:
        "conda/toolz5.yaml"
    shell:
        """
        (
            echo "gene_id\ttally" ;\
            sort -m \
                <(lz4cat {input.mgen} | sed '1,1d' | sort) \
                <(lz4cat {input.sim}  | sed '1,1d' | sort) \
            | awk -v OFS="\t" '(NR==1){{gene=$1; tally=$2}} (NR>1 && $1==gene){{tally+=$2}} (NR>1 && $1!=gene) {{print gene,tally; gene=$1; tally=$2}}  END{{print gene,tally}}' \
        ) | lz4 -c > {output}
        """


rule spikein_simulated_benchmark_isolate_to_mgen_gtpro:
    output:
        "{stemA}/reads/{mgen}_spikein_{genome}_seed{seed}_cov{cov}/{stemB}.gtpro_parse.tsv.bz2",
    input:
        script="scripts/sum_merged_gtpro_parsed_samples.py",
        mgen="{stemA}/reads/{mgen}/{stemB}.gtpro_parse.tsv.bz2",
        sim="{stemA}/reads/sim_{genome}_len150_seed{seed}_cov{cov}/{stemB}.gtpro_parse.tsv.bz2",
    shell:
        """
        {input.script} {input.mgen} {input.sim} {output}
        """


rule write_hmp2_spikein_benchmark_known_samples:
    output:
        partition="data/group/hmp2_spikein_benchmark/species/sp-102506/ecoli-spiked.strain_samples.tsv",
    params:
        table="""
CSM79HGX_G116737_spikein_Escherichia-coli-GCF_030205875-1_seed0_cov01	Escherichia-coli-GCF_030205875-1
CSM79HGZ_G116733_spikein_Escherichia-coli-GCF_030205875-1_seed1_cov02	Escherichia-coli-GCF_030205875-1
CSM79HOV_G110595_spikein_Escherichia-coli-GCF_030205875-1_seed2_cov04	Escherichia-coli-GCF_030205875-1
CSM79HOX_G110583_spikein_Escherichia-coli-GCF_030205875-1_seed3_cov08	Escherichia-coli-GCF_030205875-1
CSM79HOZ_G116676_spikein_Escherichia-coli-GCF_030205875-1_seed4_cov16	Escherichia-coli-GCF_030205875-1
CSM79HGX_G116737_spikein_Escherichia-coli-GCF_030205145-1_seed0_cov01	Escherichia-coli-GCF_030205145-1
CSM79HGZ_G116733_spikein_Escherichia-coli-GCF_030205145-1_seed1_cov02	Escherichia-coli-GCF_030205145-1
CSM79HOV_G110595_spikein_Escherichia-coli-GCF_030205145-1_seed2_cov04	Escherichia-coli-GCF_030205145-1
CSM79HOX_G110583_spikein_Escherichia-coli-GCF_030205145-1_seed3_cov08	Escherichia-coli-GCF_030205145-1
CSM79HOZ_G116676_spikein_Escherichia-coli-GCF_030205145-1_seed4_cov16	Escherichia-coli-GCF_030205145-1
CSM79HGX_G116737_spikein_Escherichia-coli-GCF_030202075-1_seed0_cov01	Escherichia-coli-GCF_030202075-1
CSM79HGZ_G116733_spikein_Escherichia-coli-GCF_030202075-1_seed1_cov02	Escherichia-coli-GCF_030202075-1
CSM79HOV_G110595_spikein_Escherichia-coli-GCF_030202075-1_seed2_cov04	Escherichia-coli-GCF_030202075-1
CSM79HOX_G110583_spikein_Escherichia-coli-GCF_030202075-1_seed3_cov08	Escherichia-coli-GCF_030202075-1
CSM79HOZ_G116676_spikein_Escherichia-coli-GCF_030202075-1_seed4_cov16	Escherichia-coli-GCF_030202075-1
CSM79HGX_G116737_spikein_Escherichia-coli-GCF_030198905-1_seed0_cov01	Escherichia-coli-GCF_030198905-1
CSM79HGZ_G116733_spikein_Escherichia-coli-GCF_030198905-1_seed1_cov02	Escherichia-coli-GCF_030198905-1
CSM79HOV_G110595_spikein_Escherichia-coli-GCF_030198905-1_seed2_cov04	Escherichia-coli-GCF_030198905-1
CSM79HOX_G110583_spikein_Escherichia-coli-GCF_030198905-1_seed3_cov08	Escherichia-coli-GCF_030198905-1
CSM79HOZ_G116676_spikein_Escherichia-coli-GCF_030198905-1_seed4_cov16	Escherichia-coli-GCF_030198905-1
CSM79HGX_G116737_spikein_Escherichia-coli-GCF_030204715-1_seed0_cov01	Escherichia-coli-GCF_030204715-1
CSM79HGZ_G116733_spikein_Escherichia-coli-GCF_030204715-1_seed1_cov02	Escherichia-coli-GCF_030204715-1
CSM79HOV_G110595_spikein_Escherichia-coli-GCF_030204715-1_seed2_cov04	Escherichia-coli-GCF_030204715-1
CSM79HOX_G110583_spikein_Escherichia-coli-GCF_030204715-1_seed3_cov08	Escherichia-coli-GCF_030204715-1
CSM79HOZ_G116676_spikein_Escherichia-coli-GCF_030204715-1_seed4_cov16	Escherichia-coli-GCF_030204715-1
"""
    shell:
        """
cat <<EOF > {output}
{params.table}
EOF
        """


use rule run_spgc as run_spgc_with_known_strain_samples with:
    output:
        "data/group/hmp2_spikein_benchmark/species/sp-102506/{proc_stem}.gtpro.gene{centroidA}_{dbv}-{pang_stem}-agg{centroidB}.spgc_specgene-{specgene}_ss-ecoli-spiked_t-10_thresh-corr{cthresh}-depth{dthresh}.nc",
    input:
        depth="data/group/hmp2_spikein_benchmark/species/sp-102506/{proc_stem}.gene{centroidA}_{dbv}-{pang_stem}-agg{centroidB}.depth2.tsv.gz",
        partition="data/group/hmp2_spikein_benchmark/species/sp-102506/ecoli-spiked.strain_samples.tsv",
        species_genes="data/species/sp-102506/midasdb.gene{centroidB}_{dbv}.spgc_specgene-{specgene}.species_gene.list",


use rule aggregate_strain_metagenotype as aggregate_strain_metagenotype_with_known_strain_samples with:
    output:
        "data/group/hmp2_spikein_benchmark/species/sp-102506/{stem}.gtpro.spgc_ss-ecoli-spiked.mgtp.nc",
    input:
        script="scripts/aggregate_strain_metagenotypes_across_strain_samples.py",
        mgtp="data/group/hmp2_spikein_benchmark/species/sp-102506/{stem}.gtpro.mgtp.nc",
        mapping="data/group/hmp2_spikein_benchmark/species/sp-102506/ecoli-spiked.strain_samples.tsv",


use rule compile_spgc_strain_metadata as compile_spgc_strain_metadata_with_known_strain_samples with:
    output:
        "data/group/hmp2_spikein_benchmark/species/sp-102506/{stemA}.gtpro.gene{gene_params}.spgc_specgene-{specgene_params}_ss-ecoli-spiked_t-{t}_thresh-{thresh_params}.strain_meta.tsv",
    input:
        script="scripts/compile_spgc_results_metadata.py",
        agg_mgtp="data/group/hmp2_spikein_benchmark/species/sp-102506/{stemA}.gtpro.spgc_ss-ecoli-spiked.mgtp.nc",
        spgc="data/group/hmp2_spikein_benchmark/species/sp-102506/{stemA}.gtpro.gene{gene_params}.spgc_specgene-{specgene_params}_ss-ecoli-spiked_t-{t}_thresh-{thresh_params}.nc",


use rule collect_filtering_metadata as collect_filtering_metadata_with_known_strain_samples with:
    output:
        "data/group/hmp2_spikein_benchmark/species/sp-102506/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-ecoli-spiked_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/hmp2_spikein_benchmark/species/sp-102506/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-ecoli-spiked_t-{t}_thresh-{thresh}.strain_meta.tsv",


ruleorder: collect_filtering_metadata_with_known_strain_samples > collect_filtering_metadata


use rule assess_infered_strain_accuracy_uhgg_tiles as assess_infered_strain_accuracy_emapper_unit_with_known_strain_samples with:
    output:
        "data/group/hmp2_spikein_benchmark/species/sp-102506/{stemA}.gene{pangenome_params}.{stemB}.{strain}.{unit}-reconstruction_accuracy.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/hmp2_spikein_benchmark/species/sp-102506/{stemA}.gene{pangenome_params}.{stemB}.{unit}-strain_gene.tsv",
        truth="data/species/sp-102506/genome/{strain}.prodigal-single.cds.emapper.{unit}-strain_gene.tsv",


rule download_sra_run_for_biosample_as_raw_fastq:
    output:
        r1="raw/genomes/{db_name}/{genome}/r1.fq.gz",
        r2="raw/genomes/{db_name}/{genome}/r2.fq.gz",
    wildcard_constraints:
        db_name="ncbi|gtdb",
    log:
        esearch="raw/genomes/{db_name}/{genome}/esearch_sra.log",
    conda:
        "conda/sra_tools.yaml"
    resources:
        network_connections=1,
    params:
        biosample=lambda w: config["genome"].loc[w.genome].ncbi_assembly_biosample,
    shell:
        """
        outdir=$(dirname {output.r1})
        esearch -db sra -query "{params.biosample}" | esummary > {log.esearch}
        srr=$(xmllint --xpath 'string(//Runs/Run/@acc)' {log.esearch})
        fasterq-dump --outdir ${{outdir}} ${{srr}}
        gzip -c ${{outdir}}/${{srr}}_1.fastq > {output.r1}
        gzip -c ${{outdir}}/${{srr}}_2.fastq > {output.r2}
        rm ${{outdir}}/${{srr}}_1.fastq ${{outdir}}/${{srr}}_2.fastq
        """
