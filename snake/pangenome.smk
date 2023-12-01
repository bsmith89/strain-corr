rule collect_species_pangenome_centroids_new:
    output:
        fn="data/species/sp-{species}/pangenome{centroid}_{dbv}.bt2.d/centroids.fn",
    input:
        fn="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/centroids.ffn",
        gene_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
    params:
        col=lambda w: {"99": 2, "95": 3, "90": 4, "85": 5, "80": 6, "75": 7}[w.centroid],
    conda:
        "conda/seqtk.yaml"
    shell:
        """
        seqtk subseq \
                {input.fn} \
                <(cat {input.gene_info} | cut -f{params.col} | sort | uniq) \
            > {output.fn}
        """


rule collect_multispecies_pangenome_centroids_to_hash:
    output:
        fn="data/hash/{_hash}/pangenomes{centroid}_{dbv}.bt2.d/centroids.fn",
    input:
        fasta=lambda w: [
            f"ref/midasdb_uhgg_{w.dbv}/pangenomes/{species}/centroids.ffn"
            for species in config["hash_to_species_set"][w._hash]
        ],
        gene_info=lambda w: [
            f"ref/midasdb_uhgg_{w.dbv}/pangenomes/{species}/gene_info.txt"
            for species in config["hash_to_species_set"][w._hash]
        ],
    params:
        fasta_pattern="ref/midasdb_uhgg_{dbv}/pangenomes/$species/centroids.ffn",
        gene_info_pattern="ref/midasdb_uhgg_{dbv}/pangenomes/$species/gene_info.txt",
        species_list=lambda w: config["hash_to_species_set"][w._hash],
        col=lambda w: {"99": 2, "95": 3, "90": 4, "85": 5, "80": 6, "75": 7}[w.centroid],
    resources:
        walltime_hr=12,
    conda:
        "conda/seqtk.yaml"
    shell:
        """
        echo "Collecting representative sequences."
        for species in {params.species_list}
        do
            fasta_path={params.fasta_pattern}
            gene_info_path={params.gene_info_pattern}
            seqtk subseq \
                    $fasta_path \
                    <(cat $gene_info_path | cut -f{params.col} | sort | uniq)
            echo -n "$species " >&2
        done \
            > {output.fn}
        echo "" >&2
        """


# NOTE: Hub-rule?
# TODO: Move this to "general.smk".
rule build_large_bowtie_index:
    output:
        db="data/hash/{stem}.bt2db",
        fn="data/hash/{stem}.bt2db.fn",
        fai="data/hash/{stem}.bt2db.fn.fai",
        db1="data/hash/{stem}.bt2db.1.bt2l",
        db2="data/hash/{stem}.bt2db.2.bt2l",
        db3="data/hash/{stem}.bt2db.3.bt2l",
        db4="data/hash/{stem}.bt2db.4.bt2l",
        dbrev1="data/hash/{stem}.bt2db.rev.1.bt2l",
        dbrev2="data/hash/{stem}.bt2db.rev.2.bt2l",
        md5="data/hash/{stem}.bt2db.checksum.md5",
    input:
        fn="data/hash/{stem}.fn",
    conda:
        "conda/midas.yaml"
    threads: 48
    shell:
        """
        echo "Linking input reference." >&2
        ln -rs {input.fn} {output.fn}

        echo "Building bowtie2 reference index." >&2
        bowtie2-build  \
                --large-index \
                --threads {threads} \
                --seed 0 \
            {output.fn} {output.db}

        echo "Indexing reference sequences." >&2
        samtools faidx {output.fn}

        echo "This file is used as a flag and has the same name as the bowtie2 database overall." > {output.db}

        echo "Collecting md5 checksums." >&2
        md5sum {output.db}* > {output.md5}
        """


# NOTE: I no longer use this original set of parameters
# but base other rules on this template.
rule run_bowtie_multispecies_pangenome_v0:
    output:
        "data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes{centroid}-v0.{bam_or_cram}",
    input:
        db="data/hash/{hash}/pangenomes{centroid}.bt2.d/centroids.bt2db",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    wildcard_constraints:
        centroid="99|95|90|85|80|75",
        bam_or_cram="bam|cram",
    params:
        extra_flags="--local --very-sensitive-local ",
        seed=0,
    conda:
        "conda/bowtie2.yaml"  # FIXME
    threads: 12
    resources:
        walltime_hr=48,
        mem_mb=100_000,
        pmem=100_000 // 12,
    shell:
        """
        bowtie2 --no-unal \
            -x {input.db} \
            --threads {threads} --mm -q \
            -U {input.r1} \
            -U {input.r2} \
            --seed {params.seed} \
            {params.extra_flags} \
            | samtools view --threads 2 -u \
            | samtools sort --threads {threads} -u -T {output} \
            | samtools view --threads 2 -O {wildcards.bam_or_cram} --reference {input.db}.fn -t {input.db}.fn.fai -o {output}
        """


use rule run_bowtie_multispecies_pangenome_v0 as run_bowtie_multispecies_pangenome_v22_new with:
    output:
        "data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes{centroid}_{dbv}-v22.{bam_or_cram}",
    input:
        db="data/hash/{hash}/pangenomes{centroid}_{dbv}.bt2.d/centroids.bt2db",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    benchmark:
        "data/hash/{hash}/reads/{mgen}/r.{proc}.pangenomes{centroid}_{dbv}-v22.{bam_or_cram}.benchmark"
    params:
        extra_flags="--ignore-quals --end-to-end --very-sensitive",
        seed=0,


rule alias_to_pangenome_sorted:
    output:
        "{stem}.pangenomes.sort.bam",
    input:
        "{stem}.pangenomes.bam",
    shell:
        alias_recipe


# TODO: Move this to "general.smk".
# TODO: Double check that --min-MQ 0 has no effect.
rule profile_mapping_depths:
    output:
        "{stem}.position_depth.tsv.lz4",
    input:
        "{stem}.bam",
    conda:
        "conda/midas.yaml"
    shell:
        """
        samtools depth -g SECONDARY -a {input} | lz4 -9zc > {output}
        """


rule profile_mapping_depths_mapq2:
    output:
        "{stem}.position_depth-mapq2.tsv.lz4",
    input:
        "{stem}.bam",
    conda:
        "conda/midas.yaml"
    shell:
        """
        samtools depth -g SECONDARY -a {input} --min-MQ 2 | lz4 -9zc > {output}
        """


# NOTE: "gene_mapping_tally" here means counting the number of bases aligned.
rule profile_pangenome_mapping_tally_aggregated_by_gene:
    output:
        "{stem}.pangenome{mapping_params}.gene_mapping_tally.tsv.lz4",
    input:
        "{stem}.pangenome{mapping_params}.bam",
    conda:
        "conda/bowtie2.yaml"
    threads: 1
    resources:
        walltime_hr=2,
    shell:
        """
        samtools depth -@ {threads} -g SECONDARY --min-MQ 0 {input} \
            | awk -v OFS='\t' '\\
                BEGIN {{gene_id="__START__"; depth_tally=0}} \\
                $1==gene_id {{depth_tally+=$3}} \\
                $1!=gene_id {{print gene_id, depth_tally; gene_id=$1; depth_tally=0}} \\
                END {{print gene_id,depth_tally}} \\
                ' \
            | (echo 'gene_id\ttally'; sed '1,1d') \
            | lz4 -9 -zc > {output}.temp
        mv {output}.temp {output}
        """


# NOTE: Hub-rule
# NOTE: Paired with it's child rule in reference_genome.smk.
rule load_one_species_pangenome2_depth_into_netcdf_new:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.depth2.nc",
    input:
        script="scripts/merge_pangenomes_depth.py",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenomes{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w, mgen=mgen, _hash=config["species_group_to_hash"][w.group]
            )
            for mgen in config["mgen_group"][w.group]
        ],
        gene_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
        gene_length="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/genes.len",
    wildcard_constraints:
        centroidA="99|95|90|85|80|75",
        centroidB="99|95|90|85|80|75",
    params:
        args=lambda w: [
            "{mgen}=data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenomes{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w, mgen=mgen, _hash=config["species_group_to_hash"][w.group]
            )
            for mgen in config["mgen_group"][w.group]
        ],
        centroidB_col=lambda w: f"centroid_{w.centroidB}",
    conda:
        "conda/toolz.yaml"
    threads: 1
    resources:
        walltime_hr=24,
        mem_mb=20_000,
        pmem=20_000 // 1,
    shell:
        """
        {input.script} {input.gene_length} {input.gene_info} {params.centroidB_col} {output} {params.args}
        """
