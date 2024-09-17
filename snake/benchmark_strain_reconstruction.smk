# NOTE: This rule takes the super long filename and turns it into a much shorter one for benchmarking
rule alias_spgc_gene_hits_for_benchmarking:
    output:
        "data/group/xjin/species/sp-{species}/{stem}.gene{pangenome_params}.spgc-fit.uhgg-strain_gene.tsv",
    input:
        source=lambda w: (
            "data/group/xjin/species/sp-{species}/{stem}.gtpro.{sfacts_stem}.gene{pangenome_params}.spgc_{spgc_stem}.uhgg-strain_gene.tsv".format(
                species=w.species,
                stem=w.stem,
                pangenome_params=w.pangenome_params,
                spgc_stem=get_spgc_stem(config, w.species, "xjin"),
                sfacts_stem=get_sfacts_stem(config, w.species, "xjin"),
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_spgc_gene_hits_for_benchmarking,


rule alias_spgc_depth_only_gene_hits_for_benchmarking:
    output:
        "data/group/xjin/species/sp-{species}/{stem}.gene{pangenome_params}.spgc-depth{thresh}.uhgg-strain_gene.tsv",
    input:
        # WARNING: This hard-coding of the spgc params might be problematic, even though most of them don't matter...
        source=lambda w: (
            "data/group/xjin/species/sp-{w.species}/{w.stem}.gtpro.{sfacts_stem}.gene{w.pangenome_params}.spgc_specgene-ref-filt-p95_ss-all_t-10_thresh-corr0-depth{w.thresh}.uhgg-strain_gene.tsv".format(
                w=w,
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, "xjin")
                ],
            )
        ),
    shell:
        alias_recipe


use rule compute_reference_and_spgc_pairwise_genotype_masked_hamming_distance as compute_benchmark_and_spgc_pairwise_genotype_masked_hamming_distance with:
    output:
        "data/group/xjin/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.genome_pdist-mask{thresh}-pseudo{pseudo}.pkl",
    input:
        script="scripts/calculate_ref_and_spgc_pairwise_genotype_masked_hamming_distance.py",
        spgc_agg_mgtp="data/group/xjin/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.mgtp.nc",
        ref_geno="data/group/xjin/species/sp-{species}/strain_genomes.gtpro.mgtp.nc",  # TODO: Confirm this is built correctly.


rule match_strains_to_genomes_based_on_genotype:
    output:
        "data/group/xjin/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.geno_matching_stats.tsv",
    input:
        script="scripts/match_strains_to_genomes.py",
        cdist="data/group/xjin/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.genome_pdist-mask10-pseudo10.pkl",
        spgc_agg_mgtp="data/group/xjin/species/sp-{species}/{stemA}.gtpro.{sfacts_params}.spgc_ss-all.mgtp.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.cdist} {input.spgc_agg_mgtp} {output}"


# NOTE: (2023-06-20) UHGG accuracy gets its own rule, separate from the
# other units, because it's assigned based on tiling depth with a particular
# centroidA, centroidB, mapping strategy, etc.
# Also note that UHGG tiles as an assessment unit isn't great. I could try
# switching to "best hit" UHGG assignment, but I think the accuracy would look
# misleadingly low for all of the tools...
# NOTE (2024-04-04): This rule is not used directly, but it's the parent to two other rules that
# are actually used.
# TODO (2024-04-04): Consider replacing this rule with one of its children.
rule assess_infered_strain_accuracy_uhgg_tiles:
    output:
        "data/group/xjin/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.uhggtiles-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/xjin/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/group/xjin/species/sp-{species}/genome/{strain}.tiles-l100-o99.gene{pangenome_params}.gene_matching-t30.uhggtiles-strain_gene.tsv",
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
        "data/group/xjin/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.{unit}-reconstruction_accuracy.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/xjin/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{unit}-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.{unit}-strain_gene.tsv",


use rule assess_infered_strain_accuracy_uhgg_tiles as assess_infered_strain_accuracy_uhgg_best_hit with:
    output:
        "data/group/xjin/species/sp-{species}/{stemA}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.{stemB}.{strain}.uhggtop-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/xjin/species/sp-{species}/{stemA}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best-c{centroidB}.uhggtop-strain_gene.tsv",  # FIXME: convert to a strain_gene.tsv


rule collect_xjin_benchmark_accuracy_grid_for_one_species:
    output:
        touch(
            "data/group/xjin/species/sp-{species}/{stem}.ACCURACY_BENCHMARK_GRID_ONE_SPECIES.flag"
        ),
    input:
        bench=lambda w: [
            f"data/group/xjin/species/sp-{w.species}/{w.stem}.{tool}.{genome}.{unit}-reconstruction_accuracy.tsv"
            for genome, tool, unit in product(
                species_group_genomes(w.species, "xjin"),
                [
                    "spgc-fit",
                    "spgc-depth200",
                    "spanda-s6",
                    "panphlan",
                ],
                [
                    "uhggtop",
                    "eggnog",
                    "cog",
                ],
            )
        ],
    shell:
        "echo {input} > {output}"


localrules:
    collect_xjin_benchmark_accuracy_grid_for_one_species,


rule collect_xjin_benchmark_accuracy_grid_for_all_species:
    output:
        touch("data/group/xjin/{stem}.ACCURACY_BENCHMARK_GRID.flag"),
    input:
        bench=lambda w: [
            f"data/group/xjin/species/sp-{species}/{w.stem}.ACCURACY_BENCHMARK_GRID_ONE_SPECIES.flag"
            for species in config["species_group"]["xjin"]
        ],
    shell:
        "echo {input} > {output}"


localrules:
    collect_xjin_benchmark_accuracy_grid_for_all_species,


rule collect_xjin_benchmark_spgc_strain_match:
    output:
        touch("data/group/xjin/{stem}.STRAIN_MATCH_BENCHMARK_GRID.flag"),
    input:
        sfacts_match=lambda w: [
            "data/group/xjin/species/sp-{species}/{w.stem}.geno_matching_stats.tsv".format(
                    species=config["genome"].loc[genome].species_id, w=w
                )
                for genome in config["genome_group"]["xjin"]
            if config["genome"].loc[genome].species_id != "UNKNOWN"
        ],
    shell:
        "echo {input} > {output}"


localrules:
    collect_xjin_benchmark_spgc_strain_match,


# FIXME: The input is pretty messy. Might need a convenience function to collect this list of species.
rule collect_xjin_benchmark_strain_meta:
    output:
        touch(
            "data/group/xjin/r.proc.{pang_stem}.spgc-fit.STRAIN_META_BENCHMARK_GRID.flag"
        ),
    input:
        sfacts_match=lambda w: [
            "data/group/xjin/species/sp-{species}/r.proc.gtpro.sfacts-fit.{w.pang_stem}.spgc-fit.strain_meta-s95-d100-a0-pos100-std25.tsv".format(
                    species=config["genome"].loc[genome].species_id, w=w
                )
                for genome in config["genome_group"]["xjin"]
            if config["genome"].loc[genome].species_id != "UNKNOWN"
        ],
    shell:
        "echo {input} > {output}"


localrules:
    collect_xjin_benchmark_strain_meta,


# FIXME: The input is pretty messy. Might need a convenience function to collect this list of species.
# DEPRECATED: Use .SELECT_SPECIES.flag
rule collect_xjin_benchmark_species_depth:
    output:
        touch("data/group/xjin/r.proc.{pang_stem}.SPECIES_DEPTH_BENCHMARK_GRID.flag"),
    input:
        sfacts_match=lambda w: [
            "data/group/xjin/species/sp-{species}/r.proc.{w.pang_stem}.spgc_specgene-ref-filt-p95.species_depth.tsv".format(
                    species=config["genome"].loc[genome].species_id, w=w
                )
                for genome in config["genome_group"]["xjin"]
            if config["genome"].loc[genome].species_id != "UNKNOWN"
        ],
    shell:
        "echo {input} > {output}"


localrules:
    collect_xjin_benchmark_species_depth,


# TODO: Make this rule work
rule collect_xjin_benchmark_grid_files:
    output:
        touch(
            "data/group/xjin/r.proc.gtpro.sfacts-fit.{gene_stem}.spgc-fit.BENCHMARK_GRID.flag"
        ),
    input:
        "data/group/xjin/r.proc.gtpro.sfacts-fit.STRAIN_MATCH_BENCHMARK_GRID.flag",
        "data/group/xjin/r.proc.{gene_stem}.SPECIES_DEPTH_BENCHMARK_GRID.flag",
        "data/group/xjin/r.proc.{gene_stem}.ACCURACY_BENCHMARK_GRID.flag",
        "data/group/xjin/r.proc.{gene_stem}.spgc-fit.STRAIN_META_BENCHMARK_GRID.flag",


localrules:
    collect_xjin_benchmark_grid_files,
