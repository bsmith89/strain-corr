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




# NOTE: This rule takes the super long filename and turns it into a much shorter one for benchmarking
rule alias_spgc_gene_hits_for_benchmarking:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gene{pangenome_params}.spgc-fit.uhgg-strain_gene.tsv",
    input:
        source=lambda w: (
            "data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gtpro.{sfacts_stem}.gene{pangenome_params}.spgc_{spgc_stem}.uhgg-strain_gene.tsv".format(
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
            "data/group/xjin_ucfmt_hmp2/species/sp-{w.species}/{w.stem}.gtpro.{sfacts_stem}.gene{w.pangenome_params}.spgc_specgene-ref-t25-p95_ss-all_t-30_thresh-corr0-depth{w.thresh}.uhgg-strain_gene.tsv".format(
                w=w,
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, "xjin_ucfmt_hmp2")
                ],
            )
        ),
    shell:
        alias_recipe


# FIXME: Consider moving this to a different snakefile.
# FIXME: Use the hamming distance calculator script rather than repeating myself.
rule predict_sfacts_strain_gene_content_by_nearest_neighbor_matching:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.gene{centroid}.nnmatched-m{min_diss}.uhgg-strain_gene.tsv",
    input:
        script="scripts/predict_gene_content_by_nearest_neighbor.py",
        spgc_geno="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb.mgtp.nc",
        ref_gene="data/species/sp-{species}/gene{centroid}_{dbv}.reference_copy_number.nc",
        ref_meta="ref/uhgg_genomes_all_4644.tsv",
    params:
        min_geno_diss=lambda w: int(w.min_diss) / 1000,
        min_completeness=0.9,
        max_contamination=0.02,  # FIXME: Make sure I'm explainingt his procedure correctly.
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
        spgc_geno="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.mgtp.nc",
        strain_sample="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-xjin-all.strain_samples.tsv",
        species_depth="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stemA}.gene{pangen_params}.spgc{spgc_params}.species_depth.tsv",
        mgen_list="meta/XJIN_BENCHMARK/mgen.tsv",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.strain_geno} {wildcards.strain} {input.spgc_geno} {input.species_depth} {input.strain_sample} {input.mgen_list} {output}"


# NOTE: Which is better place to alias-in strain matching requirements? alias_species_specific_sfacts_comm or here?
rule alias_spgc_strain_match_for_benchmarking:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gtpro.sfacts-fit.gene{pangenome_params}.spgc_{spgc_stem}.{strain}.geno_matching_stats.tsv",
    input:
        source=lambda w: (
            "data/group/XJIN_BENCHMARK/species/sp-{species}/{stem}.gtpro.{sfacts_stem}.gene{pangenome_params}.spgc_{spgc_stem}.{strain}.geno_matching_stats.tsv".format(
                species=w.species,
                stem=w.stem,
                pangenome_params=w.pangenome_params,
                strain=w.strain,
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, "xjin_ucfmt_hmp2")
                ],
                spgc_stem=w.spgc_stem,
            ),
        ),
    shell:
        alias_recipe


# NOTE: (2023-10-18) By order the aliasing rule before the ambigious rule, I'm (counterintuitively)
# causing the aliasing to happen later in the pipeline.
ruleorder: alias_spgc_strain_match_for_benchmarking > match_strains_to_genomes_based_on_genotype


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
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.{stemB}.{strain}.uhggtop-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best-c{centroidB}.uhggtop-strain_gene.tsv",  # FIXME: convert to a strain_gene.tsv


rule xjin_accuracy_grid:
    output:
        touch("data/group/XJIN_BENCHMARK/{stem}.ACCURACY_BENCHMARK_GRID.flag"),
    input:
        bench=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{species}/{w.stem}.{tool}.{genome}.{unit}-reconstruction_accuracy.tsv"
            for species in config["species_group"]["xjin_ucfmt_hmp2"]
            for genome, tool, unit in product(
                species_group_genomes(species, "XJIN_BENCHMARK"),
                [
                    "spgc-fit",
                "spgc-depth200",
                "spanda-s4",
                "panphlan",
            ],
            [
                "uhggtop",
                "eggnog",
            ],
        )
        if species != "TODO"
        ],
    shell:
        "cat {input} > {output}"


rule xjin_spgc_strain_match_grid:
    output:
        touch("data/group/XJIN_BENCHMARK/{stem}.STRAIN_MATCH_BENCHMARK_GRID.flag"),
    input:
        sfacts_match=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{species}/{w.stem}.{genome}.geno_matching_stats.tsv"
            for species in config["species_group"]["xjin_ucfmt_hmp2"]
            for genome in species_group_genomes(species, "XJIN_BENCHMARK")
        ],
    shell:
        "cat {input} > {output}"


localrules:
    xjin_accuracy_grid,
