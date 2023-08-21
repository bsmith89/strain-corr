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


# # TODO: Parameterize this rule (or a child of this rule) by the partition of
# # the strain_samples map.
# # NOTE: I could also parameterize the main rule, instead of splitting this one
# # out, but it will require handling filenames somehow.
# use rule calculate_strain_specific_correlation_and_depth_ratio_of_genes as calculate_strain_specific_correlation_and_depth_ratio_of_genes_with_subsampling with:
#     output:
#         corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
#         depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
#     input:
#         script="scripts/calculate_strain_partitioned_gene_stats.py",
#         nospecies_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_free_samples.list",
#         strain_partitions="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.strain_samples.tsv",
#         species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",
#         gene_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",


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


rule alias_spgc_gene_hits_as_uhgg_strain_gene:
    output:
        "{stemA}.spgc{stemB}.uhgg-strain_gene.tsv",
    input:
        "{stemA}.spgc{stemB}.strain_gene.tsv",
    shell:
        alias_recipe


ruleorder: alias_spgc_gene_hits_as_uhgg_strain_gene > aggregate_uhgg_strain_gene_by_annotation


localrules:
    alias_spgc_gene_hits_as_uhgg_strain_gene,


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


rule select_dominant_xjin_sfacts_strains_for_benchmark:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.strain_gene.tsv",
    input:
        script="scripts/select_dominant_strain_in_sample_subset.py",
        hits="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.strain_gene.tsv",
        comm="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.comm.tsv",
        sample_list="meta/XJIN_BENCHMARK/mgen.tsv",
    shell:
        "{input.script} {input.comm} {input.sample_list} {input.hits} > {output}"


# NOTE: (2023-06-13) This rule order should mean that
# alias_spgc_gene_hits_as_uhgg_strain_gene doesn't look for an XJIN_BENCHMARK
# inference first.
ruleorder: select_dominant_xjin_sfacts_strains_for_benchmark > alias_spgc_gene_hits_as_uhgg_strain_gene > select_strain_gene_hits


# NOTE: (2023-06-13) This renaming brings spgc inferences into alignment with spanda and panphlan so
# that the accuracy assessment rules below can be generic between all three.
rule alias_spgc_xjin_inferences_to_match_other_tools:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/r.proc.gene{pangenome_params}.spgc-fit.strain_gene.tsv",
    input:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene{pangenome_params}.spgc-fit.strain_gene.tsv",
    shell:
        alias_recipe


localrules:
    alias_spgc_xjin_inferences_to_match_other_tools,


# NOTE: (2023-06-20) UHGG accuracy gets its own rule, separate from the
# other units, because it's assigned based on tiling depth with a particular
# centroidA, centroidB, mapping strategy, etc.
# Also note that UHGG as an assessment unit isn't great. I could try switching
# to "best hit" UHGG assignment, but I think the accuracy would look
# misleadingly low...
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


# # NOTE: This hard-codes the "xjin_" prefix and finds the strain with the most
# # of these samples.
# # TODO: (2023-06-13) Decide if I should be selecting the panphlan/spanda strain
# # based on something like this.
# rule match_xjin_strains:
#     output:
#         samples="data/group/xjin_hmp2/{stem}.spgc_ss-{ss_params}.xjin_sample_count.tsv",
#     input:
#         samples="data/group/xjin_hmp2/{stem}.spgc_ss-{ss_params}.strain_samples.tsv",
#     params:
#         pattern="^xjin_"
#     shell:
#         """
#         awk '$2~/{params.pattern}/ {{print $0}}' {input} \
#                 | cut -f1 \
#                 | sort \
#                 | uniq -c \
#                 | awk -v OFS="\t" '{{print $2,$1}}' \
#             > {output}
#         """


# NOTE: This rule is required because multiple reference genomes may be present for a given species.
rule compile_reference_genome_accuracy_info:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.{unit}-xjin_strain_summary.tsv",
    input:
        script="scripts/compile_reference_genome_accuracy_info2.py",
        reference_genome_accuracy=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stem}.{genome}.{w.unit}-reconstruction_accuracy.tsv"
            for genome in species_genomes(w.species)
        ],
    params:
        args=lambda w: [
            f"{genome}=data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stem}.{genome}.{w.unit}-reconstruction_accuracy.tsv"
            for genome in species_genomes(w.species)
        ],
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {params.args} > {output}"




rule xjin_compare_tool_accuracy_single_species_single_unit:
    output:
        touch(
            "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.{unit}-accuracy.ALL_TOOLS.flag"
        ),
    input:
        spgc="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.spgc-fit.{unit}-xjin_strain_summary.tsv",
        spanda="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.spanda-s2.{unit}-xjin_strain_summary.tsv",
        panphlan="data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.panphlan.{unit}-xjin_strain_summary.tsv",
    shell:
        "echo {input} > {output}"


localrules:
    xjin_compare_tool_accuracy_single_species_single_unit,


rule xjin_accuracy_grid:
    output:
        touch("data/group/XJIN_BENCHMARK/{stem}.ACCURACY_BENCHMARK_GRID.flag"),
    input:
        lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{species}/{w.stem}.{tool}.{unit}-xjin_strain_summary.tsv"
            for species, tool, unit in product(
                config["genome"].species_id.unique(),
                [
                    "spgc-fit",
                "spanda-s2",
                "panphlan",
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


