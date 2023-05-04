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


rule run_panphlan_spanda:
    output:
        tally="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_spanda_map.tsv",
    input:
        pangenome="ref/strainpanda/{species}",
        r1="data/reads/{mgen}/r1.{stem}.fq.gz",
        r2="data/reads/{mgen}/r2.{stem}.fq.gz",
    params:
        panphlan_species=lambda w: config["species_to_spanda"][w.species],
    conda:
        "conda/panphlan.yaml"
    threads: 2
    resources:
        walltime_hr=20,
    shell:
        """
        fastq=$(mktemp)
        echo Unzipping {input.r1} and {input.r2} to $fastq >&2
        gzip -dc {input.r1} {input.r2} > $fastq
        echo Running panphlan_map on {input.r1} and {input.r2} as $fastq >&2
        panphlan_map.py --nproc {threads} \
                -p {input.pangenome}/panphlan_{params.panphlan_species}_pangenome.csv \
                --indexes {input.pangenome}/panphlan_{params.panphlan_species} \
                -i $fastq \
            -o {output.tally}
        rm $fastq
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


# FIXME: I've hard-coded xjin_hmp2 so that this table does not require
# re-running the bowtie2 building and mapping steps.
rule construct_spanda_count_matrix_from_spgc_mapping:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.spanda_counts.csv",
    input:
        script="scripts/construct_spanda_count_matrix.py",
        samples=lambda w: [
            f"data/group/xjin_hmp2/species/sp-{w.species}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        sample_args=lambda w: [
            f"{mgen}=data/group/xjin_hmp2/species/sp-{w.species}/reads/{mgen}/{w.stem}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
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
    singularity:
        config["container"]["spanda"]
    threads: 12
    shell:
        """
        tmpdir=$(mktemp -d)
        ln -rs {input.pangenome} $tmpdir/{wildcards.species}_pangenome.csv
        echo $tmpdir/{wildcards.species}_pangenome.csv
        Rscript include/StrainPanDA/bin/run_strainpandar.r \
                -c {input.data} \
                -r $tmpdir \
                -o {params.outstem} \
                -t {threads} \
                -m {params.max_strains} \
                -n {params.expect_strains}
        """
