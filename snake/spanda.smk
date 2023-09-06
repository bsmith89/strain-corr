# StrainPanDA: https://github.com/xbiome/StrainPanDA


rule start_shell_spanda:
    singularity:
        config["container"]["spanda"]
    shell:
        "bash"


rule download_spanda_reference:
    output:
        "raw/ref/strainpanda/{species}.tar.gz",
    params:
        url=lambda w: "https://zenodo.org/record/6592017/files/"
        + config["species_to_spanda"][w.species]
        + ".tar.gz",
    shell:
        curl_recipe


rule extract_spanda_reference:
    output:
        directory("ref/strainpanda/{species}"),
    input:
        "raw/ref/strainpanda/{species}.tar.gz",
    params:
        dir="raw/ref/strainpanda",
        spanda_species=lambda w: config["species_to_spanda"][w.species],
    shell:
        """
        tar -xzf {input} --directory {params.dir}
        mv {params.dir}/{params.spanda_species} {output}
        """


# NOTE: (before 2023-06-13) I've hard-coded xjin_hmp2 -> XJIN_BENCHMARK so that
# this table does not require re-running the bowtie2 building and mapping
# steps.
# This means that I process all xjin samples but take their counts from
# xjin_hmp2 pangenome profiling.
rule construct_spanda_count_matrix_from_spgc_mapping_xjin_benchmark_new:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.spanda_counts.csv",
    input:
        script="scripts/merge_pangenomes_tallies_for_spanda.py",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                mgen=mgen,
                _hash=config["species_group_to_hash"]["xjin_ucfmt_hmp2"],
            )
            for mgen in config["mgen_group"]["xjin"]
        ],
        gene_info="ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt",
    params:
        sample_args=lambda w: [
            "{mgen}=data/hash/{_hash}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                mgen=mgen,
                _hash=config["species_group_to_hash"]["xjin_ucfmt_hmp2"],
            )
            for mgen in config["mgen_group"]["xjin"]
        ],
    shell:
        """
        {input.script} {input.gene_info} {output} {params.sample_args}
        """




rule run_spanda_decompose_new:
    output:
        gene="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}_new-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}.genefamily_strain.csv",
        sample="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}_new-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}.strain_sample.csv",
    input:
        data="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}_new-{bowtie_params}.spanda_counts.csv",
        pangenome="data/species/sp-{species}/midasdb_uhgg_pangenome{centroidB}_new.tsv",
    params:
        outstem="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}_new-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}",
        max_strains=lambda w: int(w.nstrain),
        expect_strains=lambda w: int(w.nstrain),
        libstrainpandar="include/StrainPanDA/src/strainpandar",
    singularity:
        config["container"]["spanda"]
    threads: 1
    resources:
        walltime_hr=24,
    shell:
        """
        tmpdir=$(mktemp -d)
        ln -rs {input.pangenome} $tmpdir/{wildcards.species}_pangenome.csv
        echo $tmpdir/{wildcards.species}_pangenome.csv
        Rscript include/StrainPanDA/bin/run_strainpandar.r \
                --counts {input.data} \
                --reference $tmpdir \
                --output {params.outstem} \
                --threads {threads} \
                --max_rank {params.max_strains} \
                --rank {params.expect_strains} \
                --libstrainpandar {params.libstrainpandar} \
                --noextras \
                --mincov 10 --minfrac 0.9 --minreads 1e6 \
                --minsamples 1
        """


# NOTE: Converting pangenomes -> gene.
rule convert_spanda_genes_table_to_strain_gene_format:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{pangenome_params}.spanda{spanda_params}.uhgg-strain_gene.tsv",
    input:
        "data/group/{group}/species/sp-{species}/{stem}.pangenomes{pangenome_params}.spanda{spanda_params}.genefamily_strain.csv",
    run:
        data = pd.read_csv(input[0]).rename_axis(index="gene_id")
        strain_names = data.columns.to_series()[
            lambda x: x.str.startswith("strain")
        ].values
        data = data[strain_names]
        data.to_csv(output[0], sep="\t")
