# NOTE: Hub-rule
rule construct_panphlan_pangenome_metadata_from_midas_uhgg_new:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.panphlan_pangenome.tsv",
    input:
        script="scripts/construct_panphlan_pangenome_from_midas.py",
        gene_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
        cluster_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/cluster_info.txt",
    params:
        centroid=lambda w: int(w.centroid),
        fillna_gene_length=1000,
    shell:
        "{input.script} {input.gene_info} {input.cluster_info} {params.centroid} {params.fillna_gene_length} {output}"


# NOTE: This rule is equivilant to BOTH "construct_spanda_count_matrix" and
# "run_spanda_decompose" combined.
# NOTE: (before 2023-06-13) I've hard-coded xjin_hmp2 -> XJIN_BENCHMARK so that this table does not require
# re-running the bowtie2 building and mapping steps.
# This means that I process all xjin samples but take their counts from xjin_hmp2
# pangenome profiling.
# FIXME: I need a new way to build the species-specific mapping tallies for each species.
rule run_panphlan_on_spgc_mapping_xjin_benchmark_new:
    output:
        hit="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.pangenome{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.panphlan_hit.tsv",
        depth="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.pangenome{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.panphlan_depth.tsv",
        thresh="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.pangenome{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.panphlan_thresh.tsv",
    input:
        pangenome="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.panphlan_pangenome.tsv",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenome{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                mgen=mgen,
                _hash=config["species_group_to_hash"]["xjin_ucfmt_hmp2"],
            )
            for mgen in config["mgen_group"]["xjin"]
        ],
    params:
        sample_pattern=lambda w: "data/hash/{_hash}/reads/$sample/{w.stem}.pangenome{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
            w=w, _hash=config["species_group_to_hash"]["xjin_ucfmt_hmp2"]
        ),
        sample_list=lambda w: list(config["mgen_group"]["xjin"]),
        min_depth=0,
        left_max=1_000_000,
        right_min=0,
    conda:
        "conda/panphlan_dev.yaml"
    threads: 2
    resources:
        walltime_hr=12,
    shell:
        """
        tmpdir=$(mktemp -d)
        echo Reading sample depths from database into $tmpdir >&2
        for sample in {params.sample_list}
        do
            # echo $sample >&2
            lz4 -dc {params.sample_pattern} | sed '1,1d' | {{ grep -Ff <(cut -f2 {input.pangenome}) || [[ $? == 1 ]]; }} > $tmpdir/$sample
        done
        panphlan_profiling.py -i $tmpdir \
                -p {input.pangenome} \
                --left_max {params.left_max} --right_min {params.right_min} \
                --min_coverage {params.min_depth} \
                --o_matrix {output.hit} \
                --o_covmat {output.depth} \
                --o_idx {output.thresh}
        rm -r $tmpdir
        """


# NOTE: Converting pangenomes -> gene.
rule convert_panphlan_genes_table_to_strain_gene_format:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{pangenome_params}.panphlan.uhgg-strain_gene.tsv",
    input:
        "data/group/{group}/species/sp-{species}/{stem}.pangenomes{pangenome_params}.panphlan_hit.tsv",
    run:
        pd.read_table(input[0]).set_index("Unnamed: 0").rename_axis(
            index="gene_id"
        ).to_csv(output[0], sep="\t")
