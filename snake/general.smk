# {{{1 Generalizable rules


rule find_genes:
    output:
        "data/{stem}.prodigal.gff",
    input:
        "data/{stem}.fn",
    threads: 12
    shell:
        dd(
            """
        cat {input} | parallel --gnu --plain -j {threads} --recstart '>' --pipe \
                prodigal -m -q -g 11 -p meta -f gff \
            > {output}
        """
        )


rule convert_prodigal_gff_to_bed:
    output:
        "{stem}.prodigal.bed",
    input:
        "{stem}.prodigal.gff",
    shell:
        dd(
            r"""
        cat {input} \
                | awk -v OFS='\t' \
                        '/^#/ {{}} \
                         !/^#/ {{seq[$1]++; print $1, $4 - 1, $5, $1"_"seq[$1], 0, $7}} \
                         ' \
                | sort -k1,1n -k2,2n \
            > {output}
        """
        )


# Trivial conversion of lengths into a BED file format describing the
# full span of contigs as a feature.
rule nlength_to_bed:
    output:
        "{stem}.bed",
    input:
        "{stem}.nlength.tsv",
    shell:
        dd(
            """
        awk -v OFS='\\t' '{{print $1, "0", $2, $1, "0", "+"}}' {input} > {output}
        """
        )


rule fetch_prodigal_cds:
    output:
        "{stem}.cds.fn",
    input:
        bed="{stem}.prodigal.bed",
        fn="{stem}.fn",
        fai="{stem}.fn.fai",
    shell:
        dd(
            """
        bedtools getfasta -name -fi {input.fn} -bed {input.bed} -s \
                | awk -v FS=':' '/^>/{{print $1}} !/^>/{{print $0}}' \
            > {output}
        """
        )


rule translate_nucleotide_to_protein:
    output:
        "{stem}.tran.fa",
    input:
        "{stem}.fn",
    shell:
        dd(
            r"""
        cat {input} \
                | transeq -auto -filter -sformat1 fasta -table 11 -frame 1 -trim -stdout \
                | sed 's:\(^>.*\)_1$:\1:' \
                | sed 's:\*:x:g' \
            > {output}
        """
        )


rule denovo_align_fa:
    output:
        "{stem}.muscle.afa",
    input:
        "{stem}.fa",
    shell:
        "muscle < {input} > {output}"


rule hmmbuild:
    output:
        "{stem}.hmm",
    input:
        "{stem}.afa",
    shell:
        dd(
            """
        hmmbuild --amino --informat afa {output} {input}
        """
        )


rule infer_aa_phylogeny:
    output:
        "{stem}.prot.nwk",
    input:
        "{stem}.afa",
    shell:
        "fasttree < {input} > {output}"


rule bowtie_index_build:
    output:
        "{stem}.1.bt2",
        "{stem}.2.bt2",
        "{stem}.3.bt2",
        "{stem}.4.bt2",
        "{stem}.rev.1.bt2",
        "{stem}.rev.2.bt2",
    input:
        "{stem}.fn",
    threads: 12
    resources:
        mem_mb=int(200e3),
        pmem=int(200e3) // 12,
    shell:
        dd(
            """
        bowtie2-build --threads {threads} {input} {wildcards.stem}
        """
        )


ruleorder: bowtie_index_build > bowtie_index_build_from_gzipped


rule bowtie_index_build_from_gzipped:
    output:
        "{stem}.1.bt2",
        "{stem}.2.bt2",
        "{stem}.3.bt2",
        "{stem}.4.bt2",
        "{stem}.rev.1.bt2",
        "{stem}.rev.2.bt2",
    input:
        "{stem}.fn.gz",
    params:
        db=lambda w: f"{w.stem}",
    threads: 12
    shell:
        dd(
            """
        tmp=$(mktemp)
        echo $tmp
        zcat {input} > $tmp
        bowtie2-build --threads {threads} $tmp {params.db}
        """
        )


rule build_samtools_index:
    output:
        "{stem}.fn.fai",
    input:
        "{stem}.fn",
    shell:
        "samtools faidx {input}"


rule index_bam:
    output:
        "data/{stem}.sort.bam.bai",
    input:
        "data/{stem}.sort.bam",
    shell:
        "samtools index {input} {output}"


rule index_cram:
    output:
        "data/{stem}.sort.cram.crai",
    input:
        "data/{stem}.sort.cram",
    shell:
        "samtools index {input} {output}"


rule make_diamond_db:
    output:
        "ref/{db}.dmnd",
    input:
        "ref/{db}.fa",
    threads: config['MAX_THREADS']
    shell:
        dd(
            """
        diamond makedb --threads {threads} --db ref/{wildcards.db} --in {input}
        """
        )


rule diamond_search_fn:
    output:
        "data/{query}.{db}-blastx.tsv",
    input:
        fasta="data/{query}.fn",
        db="ref/{db}.dmnd",
    params:
        db="ref/{db}",
    threads: config['MAX_THREADS']
    shell:
        dd(
            """
        diamond blastx --threads {threads} --db {params.db} --query {input.fasta} > {output}
        """
        )


rule diamond_search_fa:
    output:
        "data/{query}.{db}-blastp.tsv",
    input:
        fasta="data/{query}.fa",
        db="ref/{db}.dmnd",
    params:
        db="ref/{db}",
    threads: config['MAX_THREADS']
    shell:
        dd(
            """
        diamond blastp --threads {threads} --db {params.db} --query {input.fasta} > {output}
        """
        )


# TODO: Do I actually have to feed seqtk into diamond?  Surely diamond can read
# a .fq.gz, right?
rule diamond_search_fastq:
    output:
        "data/{query}.{db}-blastqx.tsv.gz",
    input:
        fasta="data/{query}.fq.gz",
        db="ref/{db}.dmnd",
    params:
        db="ref/{db}",
    threads: 4
    resources:
        mem_mb=int(1.6e4),
        pmem=int(1.6e4) // 4,
        walltime_hr=12,
    shell:
        dd(
            """
        tmpdir=$(mktemp -d)
        echo Using temporary directory $tmpdir for {output}
        seqtk seq -A {input.fasta} \
                | diamond blastx --tmpdir $tmpdir --threads {threads} \
                      --max-target-seqs 1 --db {params.db} \
                | gzip -c \
            > {output}
        """
        )


rule hmmpress:
    output:
        h3m="ref/hmm/{model}.hmm.h3m",
        h3i="ref/hmm/{model}.hmm.h3i",
        h3f="ref/hmm/{model}.hmm.h3f",
        h3p="ref/hmm/{model}.hmm.h3p",
    input:
        hmm="ref/hmm/{model}.hmm",
    shell:
        "hmmpress -f {input}"


rule parse_hmmsearch_tblout:
    output:
        "data/{stem}.{model}-hmmer-{hmm_cutoff}.tsv",
    input:
        "data/{stem}.{model}-hmmer-{hmm_cutoff}.tblout",
    shell:
        dd(
            """
        grep -v '^#' {input} | sed 's:\s\+:\t:g' | cut -f1,3,6 > {output}
        """
        )


rule squeeze_alignment:
    output:
        "{stem}.sqz.afa",
    input:
        script="scripts/squeeze_alignment.py",
        seq="{stem}.afa",
    shell:
        "{input.script} '-.*abcdefghijklmnopqrstuvwxyz' < {input.seq} > {output}"


rule gblocks_afa:
    output:
        "{stem}.gb.afa",
    input:
        script="scripts/Gblocks.py",
        seq="{stem}.afa",
    shell:
        dd(
            """
        {input.script} < {input.seq} > {output}
        """
        )


rule count_seq_lengths_nucl:
    output:
        "{stem}.nlength.tsv",
    input:
        script="scripts/count_seq_lengths.py",
        seqs="{stem}.fn",
    shell:
        dd(
            """
        {input.script} {input.seqs} > {output}
        """
        )


rule count_seq_lengths_prot:
    output:
        "{stem}.alength.tsv",
    input:
        script="scripts/count_seq_lengths.py",
        seqs="{stem}.fa",
    shell:
        dd(
            """
        {input.script} {input.seqs} > {output}
        """
        )


# TODO: Figure out if 0M is the correct amount of overlap between edges.
# No!  It's not.  It should be k_max
rule fastg_to_gfa:
    output:
        "{stem}.a-k{k}.gfa",
    input:
        "{stem}.a-k{k}.fg",
    params:
        k=lambda w: int(w.k),
    shell:
        "fastg2gfa {input} | sed 's:\<NODE_\([0-9]\+\)_[^\\t]*\>:\\1:g' | sed 's:\\<0M\\>:{params.k}M:' > {output}"


rule gfa_to_fn:
    output:
        "{stem}.a-k{k}.fn",
    input:
        "{stem}.a-k{k}.gfa",
    shell:
        dd(
            """
        awk '/^S/{{print ">"$2"\\n"$3}}' < {input} > {output}
        """
        )


rule select_top_blastp_hits:
    output:
        "data/{stem}.{ref}-blastp.top.tsv",
    input:
        "data/{stem}.{ref}-blastp.tsv",
    shell:
        dd(
            """
        sort -k1,1 -k12,12rn {input} | sort -k1,1 -u | sort -k2,2 > {output}
        """
        )


rule select_top_blastx_hits:
    output:
        "data/{stem}.{ref}-blastx.top.tsv",
    input:
        "data/{stem}.{ref}-blastx.tsv",
    shell:
        dd(
            """
        sort -k1,1 -k12,12rn {input} | sort -k1,1 -u | sort -k2,2 > {output}
        """
        )


rule fetch_top_blastx_hits:
    output:
        "data/{stem}.{ref}-blastx.top.fn",
    input:
        tsv="data/{stem}.{ref}-blastx.top.tsv",
        fn="data/{stem}.fn",
    shell:
        dd(
            """
        seqtk subseq {input.fn} <(cut -f1 {input.tsv}) > {output}
        """
        )


rule fetch_top_blastp_hits:
    output:
        "data/{stem}.{ref}-blastp.top.fa",
    input:
        tsv="data/{stem}.{ref}-blastp.top.tsv",
        fa="data/{stem}.fa",
    shell:
        dd(
            """
        seqtk subseq {input.fa} <(cut -f1 {input.tsv}) > {output}
        """
        )
