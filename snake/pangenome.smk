rule run_bowtie_multi_species_pangenome:
    output:
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes.bam",
    input:
        bt2_dir="data/group/{group}/r.{proc}.pangenomes",
        r1="data/reads/{mgen}/r1.{proc}.fq.gz",
        r2="data/reads/{mgen}/r2.{proc}.fq.gz",
    conda:
        "conda/midas.yaml"
    threads: 5
    resources:
        walltime_hr=24,
        mem_mb=100_000,
        pmem=100_000 // 5,
    shell:
        """
        bowtie2 --no-unal -X 5000  --local --very-sensitive-local \
            -x {input.bt2_dir}/pangenomes \
            --threads {threads} --mm -q \
            -1 {input.r1} \
            -2 {input.r2} \
            | samtools view --threads 1 -b - \
            | samtools sort --threads {threads} -o {output}.temp
        mv {output}.temp {output}
        """


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


rule profile_pangenome_depth_aggregated_by_gene:
    output:
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes.gene_depth.tsv.bz2",
    input:
        "data/group/{group}/reads/{mgen}/r.{proc}.pangenomes.bam",
    conda:
        "conda/midas.yaml"
    resources:
        walltime_hr=2,
    shell:
        """
        samtools depth -a {input} \
            | awk -v OFS='\t' '\\
                BEGIN {{gene_id="__START__"; position_tally=0; depth_sum=0}} \\
                $1==gene_id {{position_tally+=1; depth_sum+=$3}} \\
                $1!=gene_id {{print gene_id,depth_sum / position_tally; gene_id=$1; position_tally=0; depth_sum=0}} \\
                END {{print gene_id,depth_sum / position_tally}} \\
                ' \
            | (echo 'gene_id\tdepth'; sed '1,1d') \
            | bzip2 -zc > {output}.temp
        mv {output}.temp {output}
        """


rule build_pangenome_depths_db:
    output:
        "data/group/{group}/r.{proc}.pangenomes.db",
    input:
        samples=lambda w: [
            f"data/group/{w.group}/reads/{mgen}/r.{w.proc}.pangenomes.gene_depth.tsv.bz2"
            for mgen in config["mgen_group"][w.group]
        ],
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
        sample_pattern="data/group/{group}/reads/$sample/r.{proc}.pangenomes.gene_depth.tsv.bz2",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
        gene_pattern="ref/midasdb_uhgg_pangenomes/$species/gene_info.txt.lz4",
        species_list=lambda w: checkpoint_select_species(
            f"data/group/{w.group}/r.{w.proc}.gtpro.horizontal_coverage.tsv",
            cvrg_thresh=0.2,
            num_samples=2,
            require_in_species_group=True,
        ),
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
        CREATE TABLE gene (
            gene_id TEXT PRIMARY KEY
          , species TEXT
        );
        CREATE TABLE sample_x_gene (
            sample TEXT
          , gene_id TEXT REFERENCES gene(gene_id)
          , depth FLOAT
          , PRIMARY KEY (sample, gene_id)
        );
        EOF

        echo "Compiling gene lists."
        for species in {params.species_list}
        do
            path="{params.gene_pattern}"
            echo -n '.' >&2
            lz4 -dc $path | sed '1,1d' | cut -f2 | sort | uniq | awk -v OFS='\t' -v species=$species '{{print $1,species}}'
        done | sqlite3 -separator '\t' $db '.import /dev/stdin gene'
        echo '' >&2
        echo '' >&2

        echo "Compiling depth data."
        for sample in {params.sample_list}
        do
            path="{params.sample_pattern}"
            echo -n '.' >&2
            bzip2 -dc $path | awk -v OFS='\t' -v sample=$sample 'NR>1 {{print sample,$1,$2}}'
        done | sqlite3 -separator '\t' $db '.import /dev/stdin sample_x_gene'
        echo '' >&2
        echo '' >&2

        echo "DONE"
        mv $db {output}
        """
        )


rule load_one_species_pangenome_depth_into_netcdf:
    output:
        "data/group/{group}/species/sp-{species}/r.{proc}.pangenomes.gene_depth.nc",
    input:
        script="scripts/load_one_species_pangenome_depth_into_netcdf.py",
        db="data/group/{group}/r.{proc}.pangenomes.db",
    conda:
        "conda/toolz.yaml"
    threads: 1
    shell:
        dd("""
        sqlite3 -separator '\t' -header {input.db} \
                ' \
                PRAGMA journal_mode=OFF; \
                PRAGMA synchronous=OFF; \
                PRAGMA locking_mode=EXCLUSIVE; \
                PRAGMA temp_store=MEMORY; \
                PRAGMA cache_size=1000000; \
                SELECT sample, gene_id, depth \
                FROM gene JOIN sample_x_gene USING (gene_id) \
                WHERE species = "{wildcards.species}"; \
                ' \
            | sed '1,4d' \
            | {input.script} {output}
        """)
