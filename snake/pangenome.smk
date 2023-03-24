rule build_multispecies_dereplicated_pangenome_bowtie_index:
    output:
        directory("data/group/{group}/r.{proc}.pangenomes{centroid}.bt2.d"),
    input:
        download_flag=lambda w: [
            f"data/species/sp-{species}/download_species_midasdb_uhgg.flag"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
        gene_info=lambda w: [
            f"ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4"
            for species in checkpoint_select_species(
                f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        fasta_pattern="ref/midasdb_uhgg/pangenomes/$species/centroids.ffn",
        gene_info_pattern="ref/midasdb_uhgg_pangenomes/$species/gene_info.txt.lz4",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.proc.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
        col=lambda w: {"99": 2, "95": 3, "90": 4, "85": 5, "80": 6, "75": 7}[w.centroid],
    conda:
        "conda/midas.yaml"
    threads: 48
    shell:
        """
        rm -rf {output}.temp
        mkdir -p {output}.temp
        for species in {params.species_list}
        do
            fasta_path={params.fasta_pattern}
            gene_info_path={params.gene_info_pattern}
            seqtk subseq \
                    $fasta_path \
                    <(lz4 -dc $gene_info_path | cut -f{params.col} | sort | uniq)
        done \
            > {output}.temp/centroids.fn
        bowtie2-build  \
                --large-index \
                --threads {threads} \
                --seed 0 \
            {output}.temp/centroids.fn {output}.temp/centroids
        mv {output}.temp {output}
        """


rule run_bowtie_multi_species_dereplicated_pangenome:
    output:
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes{centroid}.bam",
    input:
        bt2_dir="data/group/{group}/r.{proc}.pangenomes{centroid}.bt2.d",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    wildcard_constraints:
        centroid="99|95|90|85|80|75"
    params:
        extra_flags="-X 5000",
        seed=0,
    conda:
        "conda/midas.yaml"
    threads: 5
    resources:
        walltime_hr=24,
        mem_mb=100_000,
        pmem=100_000 // 5,
    shell:
        """
        bowtie2 --no-unal --local --very-sensitive-local \
            -x {input.bt2_dir}/centroids \
            --threads {threads} --mm -q \
            -1 {input.r1} \
            -2 {input.r2} \
            --seed {params.seed} \
            {params.extra_flags} \
            | samtools view --threads 1 -b - \
            | samtools sort --threads {threads} -o {output}.temp
        mv {output}.temp {output}
        """
        # TODO: Assign a fixed seed for bowtie2


# rule build_single_species_dereplicated_pangenome_bowtie_index:
#     output:
#         directory("data/species/sp-{species}/pangenome{centroid}.bt2.d"),
#     input:
#         download_flag="data/species/sp-{species}/download_species_midasdb_uhgg.flag",
#         gene_info="ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
#     params:
#         fasta="ref/midasdb_uhgg/pangenomes/{species}/centroids.ffn",
#         col=lambda w: {"99": 2, "95": 3, "90": 4, "85": 5, "80": 6, "75": 7}[w.centroid],
#     conda:
#         "conda/midas.yaml"
#     threads: 18
#     shell:
#         """
#         rm -rf {output}.temp
#         mkdir -p {output}.temp
#         seqtk subseq \
#                 {params.fasta} \
#                 <(lz4 -dc {input.gene_info} | cut -f{params.col} | sort | uniq) \
#             > {output}.temp/centroids.fn
#         bowtie2-build  \
#                 --large-index \
#                 --threads {threads} \
#                 --seed 0 \
#             {output}.temp/centroids.fn {output}.temp/centroids
#         mv {output}.temp {output}
#         """


# rule run_bowtie_single_species_pangenome:
#     output:
#         "data/species/sp-{species}/reads/{mgen}/r.{proc}.pangenome{centroid}.bam",
#     input:
#         bt2_dir="data/species/sp-{species}/pangenome{centroid}.bt2.d",
#         r1="data/reads/{mgen}/r1.{proc}.fq.gz",
#         r2="data/reads/{mgen}/r2.{proc}.fq.gz",
#     conda:
#         "conda/midas.yaml"
#     threads: 4
#     resources:
#         walltime_hr=24,
#         mem_mb=10_000,
#         pmem=10_000 // 4,
#     shell:
#         """
#         bowtie2 --no-unal -X 5000  --local --very-sensitive-local \
#             -x {input.bt2_dir}/centroids \
#             --threads {threads} --mm -q \
#             -1 {input.r1} \
#             -2 {input.r2} \
#             | samtools view --threads 1 -b - \
#             | samtools sort --threads {threads} -o {output}.temp
#         mv {output}.temp {output}
#         """


# rule profile_single_species_pangenome_depth_aggregated_by_gene:
#     output:
#         "{stem}.pangenome{centroid}.gene_depth.tsv.lz4",
#     input:
#         "{stem}.pangenome{centroid}.bam",
#     params:
#         mapq_thresh=0,  # MAPQ didn't have much effect in one test run. Also, I can't have such a high threshold if centroid=99, since multi-mapping is expected...
#     conda:
#         "conda/midas.yaml"
#     threads: 1
#     resources:
#         walltime_hr=2,
#     shell:
#         """
#         samtools depth -@ {threads} --min-MQ {params.mapq_thresh} -a {input} \
#             | awk -v OFS='\t' '\\
#                 BEGIN {{gene_id="__START__"; position_tally=0; depth_sum=0}} \\
#                 $1==gene_id {{position_tally+=1; depth_sum+=$3}} \\
#                 $1!=gene_id {{print gene_id,depth_sum / position_tally; gene_id=$1; position_tally=0; depth_sum=0}} \\
#                 END {{print gene_id,depth_sum / position_tally}} \\
#                 ' \
#             | (echo 'gene_id\tdepth'; sed '1,1d') \
#             | lz4 -9 -zc > {output}.temp
#         mv {output}.temp {output}
#         """


# rule concatenate_single_species_pangenome_depth_aggregated_by_gene:
#     output:
#         "data/group/{group}/species/sp-{species}/r.{proc}.pangenome{clust}.gene_depth.tsv.lz4",
#     input:
#         samples=lambda w: [
#             f"data/species/sp-{w.species}/reads/{mgen}/r.{w.proc}.pangenome{w.clust}.gene_depth.tsv.lz4"
#             for mgen in config["mgen_group"][w.group]
#         ],
#     params:
#         sample_pattern="data/species/sp-{species}/reads/$sample/r.{proc}.pangenome{clust}.gene_depth.tsv.lz4",
#         sample_list=lambda w: list(config["mgen_group"][w.group]),
#         header="sample	gene_id	depth",
#     shell:
#         """
#         (
#             echo "{params.header}"
#             for sample in {params.sample_list}
#             do
#                 path="{params.sample_pattern}"
#                 echo -n "." >&2
#                 lz4 -dc $path | sed '1,1d' | awk -v OFS='\t' -v sample=$sample '{{print sample,$1,$2}}'
#             done
#         ) \
#                 | lz4 -9 -zc \
#             > {output}
#         echo "" >&2
#         """


# rule load_single_species_pangenome_depth_into_netcdf:
#     output:
#         "data/group/{group}/species/sp-{species}/r.{proc}.pangenome{clust}.gene_depth.nc",
#     input:
#         script="scripts/generic_netcdf_loader.py",
#         table="data/group/{group}/species/sp-{species}/r.{proc}.pangenome{clust}.gene_depth.tsv.lz4",
#     params:
#         index="sample,gene_id",
#     shell:
#         """
#         lz4 -dc {input.table} | {input.script} {params.index} {output}
#         """


# rule profile_pangenome_depth:
#     output:
#         "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes.position_depth.tsv.bz2",
#     input:
#         "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes.bam",
#     conda:
#         "conda/midas.yaml"
#     shell:
#         """
#         samtools depth -a {input} | bzip2 -zc > {output}
#         """


rule alias_to_pangenome_sorted:
    output:
        "{stem}.pangenomes.sort.bam",
    input:
        "{stem}.pangenomes.bam",
    shell:
        alias_recipe


rule profile_pangenome_position_depths:
    output:
        "{stem}.pangenomes{bowtie_params}.position_depth-mapq{q}.tsv.lz4",
    input:
        "{stem}.pangenomes{bowtie_params}.bam",
    params:
        mapq_thresh=lambda w: int(w.q),
    conda:
        "conda/midas.yaml"
    shell:
        """
        samtools depth --min-MQ {params.mapq_thresh} -a {input} | lz4 -9zc > {output}
        """


rule profile_pangenome_depth_aggregated_by_gene:
    output:
        "{stem}.pangenomes{centroid}.gene_depth.tsv.lz4",
    input:
        "{stem}.pangenomes{centroid}.bam",
    params:
        mapq_thresh=0,
    conda:
        "conda/midas.yaml"
    threads: 1
    resources:
        walltime_hr=2,
    shell:
        """
        samtools depth -@ {threads} --min-MQ {params.mapq_thresh} -a {input} \
            | awk -v OFS='\t' '\\
                BEGIN {{gene_id="__START__"; position_tally=0; depth_sum=0}} \\
                $1==gene_id {{position_tally+=1; depth_sum+=$3}} \\
                $1!=gene_id {{print gene_id,depth_sum / position_tally; gene_id=$1; position_tally=0; depth_sum=0}} \\
                END {{print gene_id,depth_sum / position_tally}} \\
                ' \
            | (echo 'gene_id\tdepth'; sed '1,1d') \
            | lz4 -9 -zc > {output}.temp
        mv {output}.temp {output}
        """


rule profile_pangenome_mapping_tally_aggregated_by_gene:
    output:
        "{stem}.pangenomes{centroid}.gene_mapping_tally-mapq{q}.tsv.lz4",
    input:
        "{stem}.pangenomes{centroid}.bam",
    params:
        mapq_thresh=lambda w: int(w.q),
    conda:
        "conda/midas.yaml"
    threads: 1
    resources:
        walltime_hr=2,
    shell:
        """
        samtools depth -@ {threads} --min-MQ {params.mapq_thresh} {input} \
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

rule profile_pangenome_mapping_coverage_aggregated_by_gene:
    output:
        "{stem}.pangenomes{centroid}.gene_coverage_tally-mapq{q}.tsv.lz4",
    input:
        "{stem}.pangenomes{centroid}.bam",
    params:
        mapq_thresh=lambda w: int(w.q),
    conda:
        "conda/midas.yaml"
    threads: 1
    resources:
        walltime_hr=2,
    shell:
        """
        samtools depth -@ {threads} --min-MQ {params.mapq_thresh} {input} \
            | awk -v OFS='\t' '\\
                BEGIN {{gene_id="__START__"; position_tally=0}} \\
                $1==gene_id {{position_tally+=1}} \\
                $1!=gene_id {{print gene_id, position_tally; gene_id=$1; position_tally=0}} \\
                END {{print gene_id,position_tally}} \\
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


ruleorder: concatenate_reference_gene_lengths > count_seq_lengths_nucl


# FIXME: Hub rule. Commenting this out can greatly speed up
# startup time for downstream tasks.
rule build_new_pangenome_db:
    output:
        protected("data/group/{group}/r.{proc}.pangenomes{centroid}-mapq{q}.db"),
    input:
        samples=lambda w: [
            f"data/group/{w.group}/reads/{mgen}/r.{w.proc}.pangenomes{w.centroid}.gene_mapping_tally-mapq{w.q}.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
        genes="data/group/{group}/r.{proc}.pangenomes.gene_info.tsv.lz4",
        nlength="data/group/{group}/r.{proc}.pangenomes99.nlength.tsv",
    params:
        sample_pattern="data/group/{group}/reads/$sample/r.{proc}.pangenomes{centroid}.gene_mapping_tally-mapq{q}.tsv.lz4",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
    conda:
        "conda/toolz.yaml"
    threads: 1
    resources:
        walltime_hr=24,
        mem_mb=100_000,
        pmem=100_000 // 1,
    shell:
        dd(
            """
        db={output}.tmp
        echo Writing to temporary db path: $db
        rm -rf $db

        sqlite3 $db <<EOF
        CREATE TABLE _gene (
            species TEXT
          , gene_id TEXT PRIMARY KEY
          , centroid_95 TEXT
          , centroid_90 TEXT
          , centroid_85 TEXT
          , centroid_80 TEXT
          , centroid_75 TEXT
        );
        CREATE INDEX gene__species ON _gene (species);

        CREATE TABLE _gene_length (
            gene_id TEXT PRIMARY KEY REFERENCES _gene(gene_id)
          , nlength INT
        );

        CREATE VIEW gene AS
        SELECT
            species
          , gene_id
          , nlength
          , gene_id AS centroid_99
          , centroid_95
          , centroid_90
          , centroid_85
          , centroid_80
          , centroid_75
        FROM _gene
        JOIN _gene_length USING (gene_id);

        CREATE TABLE sample_x_gene (
            sample TEXT
          , gene_id TEXT REFERENCES _gene(gene_id)
          , tally INT
          , PRIMARY KEY (sample, gene_id)
        );
        EOF

        echo "Loading gene lists."
        lz4 -dc {input.genes} \
                | sed '1,1d' \
                | tqdm \
            | sqlite3 -separator '\t' $db '.import /dev/stdin _gene'

        echo "Loading gene lengths."
        cat {input.nlength} \
            | sqlite3 -separator '\t' $db '.import /dev/stdin _gene_length'

        echo "Loading mapping tally data."
        for sample in {params.sample_list}
        do
            path="{params.sample_pattern}"
            echo -n '.' >&2
            lz4 -dc $path | awk -v OFS='\t' -v sample=$sample 'NR>1 {{print sample,$1,$2}}'
        done | sqlite3 -separator '\t' $db '.import /dev/stdin sample_x_gene'
        echo '' >&2
        echo '' >&2

        echo "DONE"
        mv $db {output}
        """
        )


rule load_one_species_pangenome2_depth_into_netcdf:
    output:
        "{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-mapq{q}-agg{centroidB}.depth2.nc",
    input:
        script="scripts/load_one_species_pangenome2_depth_into_netcdf.py",
        db="{stemA}/{stemB}.pangenomes{centroidA}-mapq{q}.db",
    wildcard_constraints:
        centroidA="99|95",
        centroidB="99|95|90|85|80|75",
    conda:
        "conda/toolz.yaml"
    threads: 1
    resources:
        walltime_hr=12,
        mem_mb=20_000,
        pmem=20_000 // 1,
    shell:
        """
        {input.script} {input.db} {wildcards.species} {wildcards.centroidB} {output}
        """


# rule load_one_species_pangenome_depth_into_netcdf:
#     output:
#         "data/group/{group}/species/sp-{species}/r.{proc}.gene99.depth.nc",
#     input:
#         script="scripts/load_one_species_pangenome_depth_into_netcdf.py",
#         db="data/group/{group}/r.{proc}.pangenomes.db",
#     conda:
#         "conda/toolz.yaml"
#     threads: 1
#     resources:
#         walltime_hr=12,
#         mem_mb=20_000,
#         pmem=20_000 // 1,
#     shell:
#         """
#         {input.script} {input.db} {wildcards.species} {output}
#         """


rule normalize_multispecies_gene_depth_profile:
    output:
        "{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-agg{centroidB}.normed_depth2.nc",
    input:
        script="scripts/normalize_gene_depth_by_sample.py",
        depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-agg{centroidB}.depth2.nc",
        norm="{stemA}/{stemB}.gtpro.total_species_depth.tsv",
    shell:
        """
        {input.script} {input.norm} {input.depth} {output}
        """
