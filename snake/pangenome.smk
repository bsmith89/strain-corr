rule collect_species_pangenome_centroids_new:
    output:
        fn="data/species/sp-{species}/pangenome{centroid}_new.bt2.d/centroids.fn",
    input:
        fn="ref/midasdb_uhgg_new/pangenomes/{species}/centroids.ffn",
        gene_info="ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt",
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


rule collect_multispecies_pangenome_centroids_new:
    output:
        fn="data/group/{group}/r.{proc}.pangenomes{centroid}_new.bt2.d/centroids.fn",
    input:
        fasta=lambda w: [
            f"ref/midasdb_uhgg_new/pangenomes/{species}/centroids.ffn"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
        gene_info=lambda w: [
            f"ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        fasta_pattern="ref/midasdb_uhgg_new/pangenomes/$species/centroids.ffn",
        gene_info_pattern="ref/midasdb_uhgg_new/pangenomes/$species/gene_info.txt",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.proc.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
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


# TODO: Move this to "general.smk".
rule build_large_bowtie_index:
    output:
        db="{stem}.bt2db",
        fn="{stem}.bt2db.fn",
        fai="{stem}.bt2db.fn.fai",
        db1="{stem}.bt2db.1.bt2l",
        db2="{stem}.bt2db.2.bt2l",
        db3="{stem}.bt2db.3.bt2l",
        db4="{stem}.bt2db.4.bt2l",
        dbrev1="{stem}.bt2db.rev.1.bt2l",
        dbrev2="{stem}.bt2db.rev.2.bt2l",
        md5="{stem}.bt2db.checksum.md5",
    input:
        fn="{stem}.fn",
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
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes{centroid}-v0.{bam_or_cram}",
    input:
        db="data/group/{group}/r.{proc}.pangenomes{centroid}.bt2.d/centroids.bt2db",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    wildcard_constraints:
        centroid="99|95|90|85|80|75",
        bam_or_cram="bam|cram",
    params:
        extra_flags="--local --very-sensitive-local ",
        seed=0,
    conda:
        "conda/midas.yaml"
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
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes{centroid}_new-v22.{bam_or_cram}",
    input:
        db="data/group/{group}/r.{proc}.pangenomes{centroid}_new.bt2.d/centroids.bt2db",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    params:
        extra_flags="--ignore-quals --end-to-end --very-sensitive",
        seed=0,


use rule run_bowtie_multispecies_pangenome_v22_new as run_bowtie_species_pangenome_v22_new with:
    output:
        "data/reads/{mgen}/r.{proc}.pangenome{centroid}_new-{species}-v0.{bam_or_cram}",
    input:
        db="data/species/sp-{species}/pangenome{centroid}_new.bt2.d/centroids.bt2db",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    resources:
        walltime_hr=10,
        mem_mb=30_000,
        pmem=30_000 // 8,


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
        "conda/midas.yaml"
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


rule concatenate_pangenome_gene_info:
    output:
        "data/group/{group}/r.{proc}.pangenomes.gene_info.tsv.lz4",
    input:
        genes=lambda w: [
            f"ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        gene_pattern="ref/midasdb_uhgg_pangenomes/$species/gene_info.txt.lz4",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
        header="species_id	centroid_99	centroid_95	centroid_90	centroid_85	centroid_80	centroid_75",
    resources:
        walltime_hr=12,
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                path="{params.gene_pattern}"
                echo -n "$species " >&2
                lz4 -dc $path | sed '1,1d' | awk -v OFS='\t' -v species=$species '{{print species,$2,$3,$4,$5,$6,$7}}' | sort -k2,2 | uniq
            done
        ) \
                | lz4 -9 -zc \
            > {output}
        echo "" >&2
        """


# TODO: Figure out what downstram rules use pangenomes.gene_info.tsv.lz4 and
# make a *_new version taking this input.
rule concatenate_pangenome_gene_info_new:
    output:
        "data/group/{group}/r.{proc}.pangenomes_new.gene_info.tsv.lz4",
    input:
        genes=lambda w: [
            f"ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        gene_pattern="ref/midasdb_uhgg_new/pangenomes/$species/gene_info.txt",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
        header="species_id	centroid_99	centroid_95	centroid_90	centroid_85	centroid_80	centroid_75",
    resources:
        walltime_hr=12,
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                path="{params.gene_pattern}"
                echo -n "$species " >&2
                cat $path | sed '1,1d' | awk -v OFS='\t' -v species=$species '{{print species,$2,$3,$4,$5,$6,$7}}' | sort -k2,2 | uniq
            done
        ) \
                | lz4 -9 -zc \
            > {output}
        echo "" >&2
        """


# TODO: Generalize to other centroid definitions.
rule concatenate_reference_gene_lengths:
    output:
        "data/group/{group}/r.{proc}.pangenomes99.nlength.tsv",
    input:
        genes=lambda w: [
            f"ref/midasdb_uhgg_gene_annotations/sp-{species}.gene99_annotations.tsv"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        gene_pattern="ref/midasdb_uhgg_gene_annotations/sp-$species.gene99_annotations.tsv",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
        header="gene_id	nlength",
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                path="{params.gene_pattern}"
                echo -n "$species " >&2
                cut -f1,3 $path
            done
        ) \
            > {output}
        echo "" >&2
        """


# TODO: Generalize to other centroid definitions.
# TODO: DB-UPDATE: Take from new DB
rule concatenate_reference_gene_lengths_new:
    output:
        "data/group/{group}/r.{proc}.pangenomes99_new.nlength.tsv",
    input:
        genes=lambda w: [
            f"ref/midasdb_uhgg_new/pangenomes/{species}/genes.len"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        gene_pattern="ref/midasdb_uhgg_new/pangenomes/$species/genes.len",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
        header="gene_id	nlength",
    shell:
        """
        (
            echo "{params.header}"
            for species in {params.species_list}
            do
                path="{params.gene_pattern}"
                echo -n "$species " >&2
                cut -f1,3 $path
            done
        ) \
            > {output}
        echo "" >&2
        """


ruleorder: concatenate_reference_gene_lengths > count_seq_lengths_nucl


rule load_one_species_pangenome2_depth_into_netcdf_new:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_new-{bowtie_params}-agg{centroidB}.depth2.nc",
    input:
        script="scripts/merge_pangenomes_depth.py",
        samples=lambda w: [
            f"data/group/{w.group}/reads/{mgen}/{w.stem}.pangenomes{w.centroidA}_new-{w.bowtie_params}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
        gene_info="ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt",
        gene_length="ref/midasdb_uhgg_new/pangenomes/{species}/genes.len",
    wildcard_constraints:
        centroidA="99|95|90|85|80|75",
        centroidB="99|95|90|85|80|75",
    params:
        args=lambda w: [
            f"{mgen}=data/group/{w.group}/reads/{mgen}/{w.stem}.pangenomes{w.centroidA}_new-{w.bowtie_params}.gene_mapping_tally.tsv.lz4"
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


rule concatenate_pangenome_depth_from_samples:
    output:
        "data/group/{group}/species/sp-{species}/r.{proc}.pangenome{centroid}-{mapping_params}.gene_depth.tsv.lz4",
    input:
        samples=lambda w: [
            f"data/reads/{mgen}/r.{w.proc}.pangenome{w.centroid}-{w.species}-{w.mapping_params}.gene_depth.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        header="sample\tgene_id\tdepth",
        sample_pattern="data/reads/$sample/r.{proc}.pangenome{centroid}-{species}-{mapping_params}.gene_depth.tsv.lz4",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
    shell:
        """
        (
            echo "{params.header}"
            for sample in {params.sample_list}
            do
                path="{params.sample_pattern}"
                echo {wildcards.species} $sample >&2
                lz4 -dc $path | awk -v OFS='\t' -v sample=$sample 'NR>1 {{print sample,$1,$2}}'
            done
        ) | lz4 -9zc > {output}
        """
