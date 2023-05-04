use rule start_shell as start_shell_panphlan with:
    conda:
        "conda/panphlan.yaml"


rule construct_panphlan_pangenome_metadata_from_midas_uhgg:
    output:
        "ref/panphlan/{species}.midasdb_uhgg_pangenome{centroid}.tsv",
    input:
        script="scripts/construct_panphlan_pangenome_from_midas.py",
        gene_info="ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    params:
        centroid=lambda w: int(w.centroid),
    shell:
        "{input.script} {input.gene_info} {input.cluster_info} {params.centroid} {output}"


# NOTE: This rule is equivilant to BOTH "construct_spanda_count_matrix" and
# "run_spanda_decompose" combined.
# FIXME: I've hard-coded xjin_hmp2 so that this table does not require
# re-running the bowtie2 building and mapping steps.
rule run_panphlan_on_spgc_mapping:
    output:
        hit="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan_hit.tsv",
        depth="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan_depth.tsv",
        thresh="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan_thresh.tsv",
    input:
        pangenome="ref/panphlan/{species}.midasdb_uhgg_pangenome{centroidB}.tsv",
        samples=lambda w: [
            f"data/group/xjin_hmp2/species/sp-{w.species}/reads/{mgen}/r.{w.stem}.pangenomes{w.centroidA}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        sample_pattern="data/group/xjin_hmp2/species/sp-{species}/reads/$sample/r.{stem}.pangenomes{centroidA}-{bowtie_params}.gene_mapping_tally.tsv.lz4",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
        min_depth=0.0,
    conda:
        "conda/panphlan.yaml"
    threads: 2
    shell:
        """
        tmpdir=$(mktemp -d)
        echo Reading sample depths from database into $tmpdir >&2
        for sample in {params.sample_list}
        do
            echo $sample >&2
            lz4 -dc {params.sample_pattern} | sed '1,1d' > $tmpdir/$sample
        done
        panphlan_profiling.py -i $tmpdir \
                -p {input.pangenome} \
                --min_coverage {params.min_depth} \
                --o_matrix {output.hit} \
                --o_covmat {output.depth} \
                --o_idx {output.thresh}
        rm -r $tmpdir
        """
