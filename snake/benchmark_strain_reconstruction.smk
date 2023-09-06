# TODO: DB-UPDATE: For all XJIN_BENCHMARK rules, convert from xjin_ucfmt_hmp2.
# Also, add _new suffix where appropriate.


rule select_xjin_samples:
    output:
        "{stem}.spgc_ss-xjin-all.strain_samples.tsv",
    input:
        "{stem}.spgc_ss-all.strain_samples.tsv",
    params:
        pattern="^xjin_",
    shell:
        """
        awk 'NR==1 || $2~/{params.pattern}/ {{print $0}}' {input} > {output}
        """


# This will shuffle the order of lines in the samples map (reproducibly based on seed)
# and then take the top n_samples for each strain and write out a strain_samples
# file, but with fewer lines.
rule subsample_xjin_samples_for_benchmarking:
    output:
        "{stem}.spgc_ss-xjin-seed{seed}-n{n_samples}.strain_samples.tsv",
    input:
        script="scripts/subsample_strain_samples_for_benchmarking.py",
        samples="{stem}.spgc_ss-xjin-all.strain_samples.tsv",
    params:
        seed=lambda w: int(w.seed),
        n_samples=lambda w: int(w.n_samples),
    shell:
        """
        {input.script} {params.seed} {params.n_samples} {input.samples} {output}
        """


rule select_highest_depth_xjin_samples_for_benchmarking:
    output:
        "{stemA}.gtpro.{stemB}.spgc_ss-xjin-deepest-n{n_samples}.strain_samples.tsv",
    input:
        script="scripts/select_highest_depth_strain_samples_for_benchmarking.py",
        depth="{stemA}.gtpro.species_depth.tsv",
        samples="{stemA}.gtpro.{stemB}.spgc_ss-xjin-all.strain_samples.tsv",
    params:
        n_samples=lambda w: int(w.n_samples),
    shell:
        """
        {input.script} {params.n_samples} {input.depth} {input.samples} {output}
        """


# NOTE: This rule takes the super long filename and turns it into a much shorter one for benchmarking
rule alias_spgc_gene_hits_for_benchmarking:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gene{pangenome_params}.spgc-fit.uhgg-strain_gene.tsv",
    input:
        source=lambda w: (
            "data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gtpro.{sfacts_stem}.gene{pangenome_params}.spgc_{spgc_stem}.strain_gene.tsv".format(
                species=w.species,
                stem=w.stem,
                pangenome_params=w.pangenome_params,
                spgc_stem=config["species_group_to_spgc_stem"][
                    (w.species, "xjin_ucfmt_hmp2")
                ],
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, "xjin_ucfmt_hmp2")
                ],
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_spgc_gene_hits_for_benchmarking,


rule alias_spgc_depth_only_gene_hits_for_benchmarking:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gene{pangenome_params}.spgc-depth{thresh}.uhgg-strain_gene.tsv",
    input:
        source=lambda w: (
            "data/group/xjin_ucfmt_hmp2/species/sp-{w.species}/{w.stem}.gtpro.{sfacts_stem}.gene{w.pangenome_params}.spgc_specgene-ref-t25-p95_ss-all_t-30_thresh-corr0-depth{w.thresh}.strain_gene.tsv".format(
                w=w,
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, "xjin_ucfmt_hmp2")
                ],
            )
        ),
    shell:
        alias_recipe


# FIXME: Consider moving this to a different snakefile.
rule predict_sfacts_strain_gene_content_by_nearest_neighbor_matching:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.gene{centroid}.nnmatched-m{min_diss}.uhgg-strain_gene.tsv",
    input:
        script="scripts/predict_gene_content_by_nearest_neighbor.py",
        spgc_geno="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.strain_mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb.geno.nc",
        ref_gene="data/species/sp-{species}/gene{centroid}_new.reference_copy_number.nc",
        ref_meta="ref/uhgg_genomes_all_4644.tsv",
    params:
        min_geno_diss=lambda w: int(w.min_diss) / 1000,
        min_completeness=0.9,
        max_contamination=0.02,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.spgc_geno} {input.ref_geno} {params.min_geno_diss} {input.ref_gene} {input.ref_meta} {wildcards.species} {params.min_completeness} {params.max_contamination} {output}"


# NOTE: This rule takes the super long filename and turns it into a much shorter one for benchmarking
rule alias_nnmatched_predictions_for_benchmarking:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gene{centroidA}-{pangenome_params}-agg{centroidB}.nnmatched{nnmatch_params}.uhgg-strain_gene.tsv",
    input:
        lambda w: "data/group/xjin_ucfmt_hmp2/species/sp-{w.species}/{w.stem}.gtpro.{sfacts_stem}.gene{w.centroidB}.nnmatched{w.nnmatch_params}.uhgg-strain_gene.tsv".format(
            w=w,
            sfacts_stem=config["species_group_to_sfacts_stem"][
                (w.species, "xjin_ucfmt_hmp2")
            ],
        ),
    shell:
        alias_recipe


# ruleorder: alias_spgc_gene_hits_as_uhgg_strain_gene > aggregate_uhgg_strain_gene_by_annotation


rule aggregate_uhgg_strain_gene_by_annotation:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.{agg}-strain_gene.tsv",
    wildcard_constraints:
        agg="eggnog|top_eggnog|cog|ko",
    input:
        uhgg="data/group/{group}/species/sp-{species}/{stem}.uhgg-strain_gene.tsv",
        agg="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_{agg}.tsv",
    group:
        "assess_gene_inference_benchmark"
    run:
        gene_table = pd.read_table(input.uhgg, index_col="gene_id")
        agg = pd.read_table(input.agg, index_col="gene_id").squeeze()
        result = (
            gene_table.astype(bool).join(agg).groupby(wildcards.agg).any().astype(int)
        )
        result.to_csv(output[0], sep="\t")


rule alias_xjin_tiles_as_xjin_benchmark:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/genome/{strain}.tiles-{paramsA}.gene{pangenome_params}.gene_matching-{paramsB}.uhggtiles-strain_gene.tsv",
    input:
        "data/group/xjin_ucfmt_hmp2/species/sp-{species}/genome/{strain}.tiles-{paramsA}.gene{pangenome_params}.gene_matching-{paramsB}.uhggtiles-strain_gene.tsv",
    shell:
        alias_recipe


ruleorder: alias_xjin_tiles_as_xjin_benchmark > assign_matching_genes_based_on_tile_depth


localrules:
    alias_xjin_tiles_as_xjin_benchmark,


rule match_strains_to_genomes_based_on_genotype:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.gene{pangen_params}.spgc{spgc_params}.{strain}.geno_matching_stats.tsv",
    input:
        script="scripts/match_strains_to_genomes_in_sample_subset.py",
        strain_geno="data/species/sp-{species}/strain_genomes.gtpro.mgtp.nc",
        # TODO: Consider whether it's okay that I'm getting the spgc genotype from the full sample set, not just ss-xjin-all
        spgc_geno="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.strain_mgtp.nc",
        strain_sample="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-xjin-all.strain_samples.tsv",
        species_depth="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gene{pangen_params}.spgc{spgc_params}.species_depth.tsv",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.strain_geno} {wildcards.strain} {input.spgc_geno} {input.species_depth} {input.strain_sample} {output}"


# NOTE: (2023-06-20) UHGG accuracy gets its own rule, separate from the
# other units, because it's assigned based on tiling depth with a particular
# centroidA, centroidB, mapping strategy, etc.
# Also note that UHGG tiles as an assessment unit isn't great. I could try
# switching to "best hit" UHGG assignment, but I think the accuracy would look
# misleadingly low for all of the tools...
rule assess_infered_strain_accuracy_uhgg_tiles:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.uhggtiles-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/group/XJIN_BENCHMARK/species/sp-{species}/genome/{strain}.tiles-l100-o99.gene{pangenome_params}.gene_matching-t30.uhggtiles-strain_gene.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        """
        {input.script} \
                {input.truth} \
                {input.infer} \
                {output}
        """


use rule assess_infered_strain_accuracy_uhgg_tiles as assess_infered_strain_accuracy_emapper_unit with:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.{unit}-reconstruction_accuracy.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{unit}-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.{unit}-strain_gene.tsv",


use rule assess_infered_strain_accuracy_uhgg_tiles as assess_infered_strain_accuracy_uhgg_best_hit with:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{centroidA}_new-{bowtie_params}-agg{centroidB}.{stemB}.{strain}.uhggtop-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{centroidA}_new-{bowtie_params}-agg{centroidB}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.midas_uhgg_pangenome_new-blastn.gene_matching-best-c{centroidB}.uhggtop-strain_gene.tsv",  # FIXME: convert to a strain_gene.tsv




rule xjin_accuracy_grid:
    output:
        touch("data/group/XJIN_BENCHMARK/{stem}.ACCURACY_BENCHMARK_GRID.flag"),
    input:
        lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{species}/{w.stem}.{tool}.{genome}.{unit}-reconstruction_accuracy.tsv"
            for (genome, species), tool, unit in product(
                config["genome"].species_id.items(),
                [
                "spgc-fit",
                "spgc-depth250",
                "spanda-s2",
                "spanda-s3",
                "spanda-s4",
                "panphlan",
                "nnmatched-m0",
                "nnmatched-m1",
                "nnmatched-m10",
                "nnmatched-m50",
            ],
            [
                "uhggtop",
                "uhggtiles",
                "eggnog",
            ],
        )
        if species != "TODO"
        ],
    shell:
        "cat {input} > {output}"


localrules:
    xjin_accuracy_grid,
