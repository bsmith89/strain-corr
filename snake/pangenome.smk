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



rule load_one_species_pangenome2_depth_into_netcdf_v20:  # Hub-rule (also note child rule in reference_genome.smk)
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.depth2.nc",
    input:
        script="scripts/merge_pangenomes_depth_v20.py",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenomes{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w, mgen=mgen, _hash=config["species_group_to_hash"][w.group]
            )
            for mgen in config["mgen_group"][w.group]
        ],
        gene_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/genes_info.tsv",
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
        {input.script} {input.gene_info} {params.centroidB_col} {output} {params.args}
        """


# This comment is only needed to get the last rule off the bottom of the file.
