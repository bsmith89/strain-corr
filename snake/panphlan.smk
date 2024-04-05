rule construct_panphlan_pangenome_metadata_from_midas_uhgg_v20:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_v20.panphlan_pangenome.tsv",
    input:
        script="scripts/construct_panphlan_pangenome_from_midas_v20.py",
        gene_info="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    params:
        centroid=lambda w: int(w.centroid),
    shell:
        "{input.script} {input.gene_info} {params.centroid} {output}"


# NOTE: This rule is equivilant to BOTH "construct_spanda_count_matrix" and
# "run_spanda_decompose" combined.
rule run_panphlan_on_spgc_mapping_xjin_benchmark_v20:  # Hub-rule
    output:
        hit="data/group/xjin/species/sp-{species}/{stem}.pangenome{centroidA}_v20-{bowtie_params}-agg{centroidB}.panphlan_hit.tsv",
        depth="data/group/xjin/species/sp-{species}/{stem}.pangenome{centroidA}_v20-{bowtie_params}-agg{centroidB}.panphlan_depth.tsv",
        thresh="data/group/xjin/species/sp-{species}/{stem}.pangenome{centroidA}_v20-{bowtie_params}-agg{centroidB}.panphlan_thresh.tsv",
    input:
        pangenome="data/species/sp-{species}/midasdb.gene{centroidB}_v20.panphlan_pangenome.tsv",
        samples=lambda w: [
            "data/hash/{_hash}/reads/{mgen}/{w.stem}.pangenome{w.centroidA}_v20-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                mgen=mgen,
                _hash=config["species_group_to_hash"]["xjin"],
            )
            for mgen in config["mgen_group"]["xjin"]
        ],
    params:
        sample_pattern=lambda w: "data/hash/{_hash}/reads/$sample/{w.stem}.pangenome{w.centroidA}_v20-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
            w=w, _hash=config["species_group_to_hash"]["xjin"]
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
