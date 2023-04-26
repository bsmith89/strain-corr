use rule start_shell as start_shell_panphlan with:
    conda:
        "conda/panphlan.yaml"


rule download_panphlan_reference:
    output:
        directory("raw/ref/panphlan/{panphlan_species}"),
    params:
        outdir="raw/ref/panphlan",
    conda:
        "conda/panphlan.yaml"
    shell:
        """
        panphlan_download_pangenome.py -i {wildcards.panphlan_species} -o {params.outdir}
        """


rule link_panphlan_reference:
    output:
        directory("ref/panphlan/{species}"),
    input:
        dir=lambda w: ancient(
            "raw/ref/panphlan/" + config["species_to_panphlan"][w.species]
        ),
    shell:
        alias_recipe


rule run_panphlan:
    output:
        tally="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_map.tsv",
    input:
        pangenome="ref/panphlan/{species}",
        r1="data/reads/{mgen}/r1.{stem}.fq.gz",
        r2="data/reads/{mgen}/r2.{stem}.fq.gz",
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
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
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --indexes {input.pangenome}/{params.panphlan_species} \
                -i $fastq \
            -o {output.tally}
        rm $fastq
        """


rule collect_panphlan_single:
    output:
        hit="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_hit.tsv",
        depth="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_depth.tsv",
        thresh="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_thresh.tsv",
    input:
        pangenome="ref/panphlan/{species}",
        tally="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_map.tsv",
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
    conda:
        "conda/panphlan.yaml"
    threads: 2
    resources:
        walltime_hr=20,
    shell:
        """
        tmpdir=$(mktemp -d)
        ln -s $(realpath {input.tally}) $tmpdir/{wildcards.mgen}
        echo Running panphlan_profiling on $tmpdir/{wildcards.mgen} >&2
        panphlan_profiling.py -i $tmpdir \
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --o_matrix {output.hit} \
                --o_covmat {output.depth} \
                --o_idx {output.thresh}
        rm -r $tmpdir
        """


rule collect_panphlan_group:
    output:
        hit="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_hit.tsv",
        depth="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_depth.tsv",
        thresh="data/group/{group}/species/sp-{species}/r.{stem}.panphlan_thresh.tsv",
    input:
        pangenome="ref/panphlan/{species}",
        samples=lambda w: [
            f"data/species/sp-{w.species}/reads/{mgen}/r.{w.stem}.panphlan_map.tsv"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
        sample_pattern="data/species/sp-{species}/reads/$sample/r.{stem}.panphlan_map.tsv",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
        min_depth=0.01,
    conda:
        "conda/panphlan.yaml"
    threads: 2
    shell:
        """
        tmpdir=$(mktemp -d)
        echo linking samples into $tmpdir >&2
        for sample in {params.sample_list}
        do
            ln -rs {params.sample_pattern} $tmpdir/$sample
        done
        panphlan_profiling.py -i $tmpdir \
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --min_coverage {params.min_depth} \
                --o_matrix {output.hit} \
                --o_covmat {output.depth} \
                --o_idx {output.thresh}
        rm -r $tmpdir
        """


rule panphlan_depth_to_spgc_format:
    output:
        "data/group/{group}/species/sp-{species}/r.{proc}.geneXX-panphlan-agg90.depth2.nc"
    input:
        script="scripts/construct_spgc_input_from_panphlan.py",
        depth="data/group/{group}/species/sp-{species}/r.{proc}.panphlan_depth.tsv",
    shell:
        """
        {input.script} {input.depth} {output}
        """


rule construct_panphlan_pangenome_from_midas_uhgg:
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


rule collect_panphlan_group_from_midas_uhgg_pangenome:
    output:
        hit="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan2_hit.tsv",
        depth="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan2_depth.tsv",
        thresh="data/group/{group}/species/sp-{species}/r.{stem}.pangenomes{centroidA}-{bowtie_params}-agg{centroidB}.panphlan2_thresh.tsv",
    input:
        pangenome="ref/panphlan/{species}.midasdb_uhgg_pangenome{centroidB}.tsv",
        # NOTE: Hard-coding xjin_hmp2 here temporarily, as it includes all of xjin_102395_subset.
        samples=lambda w: [
            f"data/group/xjin_hmp2/reads/{mgen}/r.{w.stem}.pangenomes{w.centroidA}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4"
            for mgen in config["mgen_group"][w.group]
        ],
    params:
        sample_pattern="data/group/xjin_hmp2/reads/$sample/r.{stem}.pangenomes{centroidA}-{bowtie_params}.gene_mapping_tally.tsv.lz4",
        sample_list=lambda w: list(config["mgen_group"][w.group]),
        aggregate_genes_by=lambda w: {
            "99": "centroid_99",
            "95": "centroid_95",
            "90": "centroid_90",
            "85": "centroid_85",
            "80": "centroid_80",
            "75": "centroid_75",
        }[w.centroidB],
        min_depth=0.25,
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
