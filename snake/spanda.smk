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


# rule collect_panphlan_spanda_group:
#     output:
#         hit="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_spanda_hit.tsv",
#         depth="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_spanda_depth.tsv",
#         thresh="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_spanda_thresh.tsv",
#     input:
#         pangenome="ref/strainpanda/{species}",
#         samples=lambda w: [
#             f"data/species/sp-{w.species}/reads/{mgen}/r.{w.stem}.panphlan_spanda_map.tsv"
#             for mgen in config["mgen_group"][w.group]
#         ],
#     params:
#         spanda_species=lambda w: config["species_to_spanda"][w.species],
#         sample_pattern="data/species/sp-{species}/reads/$sample/r.{stem}.panphlan_spanda_map.tsv",
#         sample_list=lambda w: list(config["mgen_group"][w.group]),
#     conda:
#         "conda/panphlan.yaml"
#     threads: 2
#     shell:
#         """
#         tmpdir=$(mktemp -d)
#         echo linking samples into $tmpdir >&2
#         for sample in {params.sample_list}
#         do
#             ln -rs {params.sample_pattern} $tmpdir/$sample
#         done
#         panphlan_profiling.py -i $tmpdir \
#                 -p {input.pangenome}/panphlan_{params.spanda_species}_pangenome.csv \
#                 --o_matrix {output.hit} \
#                 --o_covmat {output.depth} \
#                 --o_idx {output.thresh}
#         rm -r $tmpdir
#         """


# NOTE: (before 2023-06-13) I've hard-coded xjin_hmp2 -> XJIN_BENCHMARK so that this table does not require
# re-running the bowtie2 building and mapping steps.
# This means that I process all xjin samples but take their counts from xjin_hmp2
# pangenome profiling.
rule construct_spanda_count_matrix_from_spgc_mapping_xjin_benchmark:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.spanda_counts.csv",
    input:
        script="scripts/construct_spanda_count_matrix.py",
        samples=lambda w: [
            f"data/group/xjin_hmp2/species/sp-{w.species}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"]["xjin"]
        ],
    params:
        sample_args=lambda w: [
            f"{mgen}=data/group/xjin_hmp2/species/sp-{w.species}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"]["xjin"]
        ],
    shell:
        """
        {input.script} {output} {params.sample_args}
        """


rule run_spanda_decompose:
    output:
        gene="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}.genefamily_strain.csv",
        sample="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}.strain_sample.csv",
    input:
        data="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}-{bowtie_params}.spanda_counts.csv",
        pangenome="ref/panphlan/{species}.midasdb_uhgg_pangenome{centroidB}.tsv",
    params:
        outstem="data/group/{group}/species/sp-{species}/{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.spanda-s{nstrain}",
        max_strains=lambda w: int(w.nstrain),
        expect_strains=lambda w: int(w.nstrain),
        libstrainpandar="include/StrainPanDA/src/strainpandar",
    singularity:
        config["container"]["spanda"]
    threads: 12
    resources:
        walltime_hr=12,
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
