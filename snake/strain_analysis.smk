rule download_cog_definitions:
    output:
        "ref/cog-20.meta.tsv",
    params:
        url="https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab",
    shell:
        curl_recipe


rule download_cog_categories:
    output:
        "ref/cog-20.categories.tsv",
    params:
        url="https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab",
    shell:
        curl_recipe


# FIXME: This rule should be implemented in SPGC itself.
rule compile_spgc_strain_metadata:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_meta.tsv",
    input:
        script="scripts/compile_spgc_results_metadata.py",
        # TODO: Take genotype QC and metadata from a purpose-built script.
        agg_mgtp="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.mgtp.nc",
        spgc="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.nc",
    params:
        ambig_thresh=0.1,
    conda:
        "conda/toolz5.yaml"
    shell:
        "{input.script} {input.spgc} {input.agg_mgtp} {params.ambig_thresh} {output}"


rule count_ref_geno_positions:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gtpro.geno.npositions.tsv",
    input:
        script="scripts/count_ref_geno_positions.py",
        geno="data/species/sp-{species}/midasdb_{dbv}.gtpro.mgtp.nc",
    params:
        ambig_thresh=0.01,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.geno} {params.ambig_thresh} {output}"


rule collect_filtering_metadata:
    output:
        "data/group/{group}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/{group}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        # spgc="data/group/{group}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.nc",
        # sample_to_strain="data/group/{group}/{stem}.spgc_ss-{ss}.strain_samples.tsv",
    params:
        min_species_genes_frac=95 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=0 / 1000,
        min_geno_positions=100,
        max_log_gene_depth_ratio_std=25 / 100,
    shell:
        "{input.script} {input.meta} {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {params.max_log_gene_depth_ratio_std} {params.min_geno_positions} {output}"


# FIXME: This is generic to either v15 or v20, but points at v15 only for
# npositions data because it should be the same as v20 and I don't want to re-run it.
rule collect_metadata_for_uhgg_ref_strains_new:
    output:
        meta="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta.tsv",
    input:
        script="scripts/extract_metadata_midasdb_v15.py",
        meta="ref/midasdb_uhgg_{dbv}/metadata/2023-11-11-genomes-all_metadata.tsv",
        genome_to_species="ref/midasdb_uhgg_{dbv}/genomes.tsv",
        pos="data/species/sp-{species}/midasdb_v15.gtpro.geno.npositions.tsv",
        genes="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
        # FIXME: Rename the above.
    shell:
        "{input.script} {input.meta} {input.genome_to_species} {input.pos} {input.genes} {wildcards.species} {output}"


rule collect_filtering_metadata_for_uhgg_ref_strains_new:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta-complete90-contam5-pos{pos}.tsv",
    input:
        script="scripts/filter_ref_strains_v15.py",
        meta="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}.strain_meta.tsv",
    params:
        min_completeness=90 / 100,
        max_contamination=5 / 100,
        min_positions=lambda w: int(w.pos),
    shell:
        "{input.script} {input.meta} {params.min_completeness} {params.max_contamination} {params.min_positions} {output}"


# TODO: Replace this rule with a genotype concatenation rule to achieve the
# same effect (harnessing the generic rule below.)
rule compute_reference_and_spgc_pairwise_genotype_masked_hamming_distance:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask{thresh}-pseudo{pseudo}.pkl",
    input:
        script="scripts/calculate_ref_and_spgc_pairwise_genotype_masked_hamming_distance.py",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb_v15.gtpro.mgtp.nc",  # FIXME: Hard-coded to avoid having to re-run for v20.
    params:
        ambiguity_threshold=lambda w: int(w.thresh) / 100,
        pseudo=lambda w: int(w.pseudo) / 10,
    conda:
        "conda/sfacts.yaml"
    resources:
        walltime_hr=12,
    shell:
        "{input.script} {input.spgc_agg_mgtp} {input.ref_geno} {params.ambiguity_threshold} {params.pseudo} {output}"


rule compute_pairwise_genotype_masked_hamming_distance:
    output:
        "{stem}.geno_pdist-mask{thresh}-pseudo{pseudo}.pkl",
    input:
        script="scripts/calculate_pairwise_genotype_masked_hamming_distance.py",
        mgtp="{stem}.mgtp.nc",
    params:
        ambiguity_threshold=lambda w: int(w.thresh) / 100,
        pseudo=lambda w: int(w.pseudo) / 10,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.mgtp} {params.ambiguity_threshold} {params.pseudo} {output}"


# SPGC and Ref genome dereplication clustering (based on genotype dissimilarity)
# NOTE: Strains are not filtered before co-clustering.
rule cocluster_reference_and_spgc_genomes_based_on_genotype_dissimilarity:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-{diss_params}.coclust-{thresh}.tsv",
    input:
        script="scripts/cluster_from_cdmat_pkl.py",
        pickle="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-{diss_params}.pkl",
    params:
        thresh=lambda w: int(w.thresh) / 1000,
    conda:
        "conda/toolz2.yaml"
    shell:
        "{input.script} {input.pickle} {params.thresh} {output}"


rule combine_all_ref_and_spgc_metadata_and_clustering:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta_spgc_and_ref.tsv",
    input:
        script="scripts/combine_all_ref_and_spgc_metadata_and_clustering.py",
        spgc="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
        ref="data/species/sp-{species}/midasdb_{dbv}.gene{centroidB}.strain_meta-complete90-contam5-pos0.tsv",
        clust="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.coclust-10.tsv",
    shell:
        "{input.script} {input.spgc} {input.ref} {input.clust} {output}"


rule compare_spgc_and_ref_dissimilarities:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_diss_spgc_and_ref.tsv",
    input:
        script="scripts/compare_spgc_and_ref_dissimilarities.py",
        meta="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta_spgc_and_ref.tsv",
        geno_pdist="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.pkl",
        uhgg_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
        eggnog_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.eggnog-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
    shell:
        "{input.script} {input.meta} {input.geno_pdist} {input.uhgg_pdist} {input.eggnog_pdist} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_ref_genomes:
    output:
        "{stem}/midasdb.gene{centroid}_{dbv}.{unit}-strain_gene.prevalence.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stem}/midasdb.gene{centroid}_{dbv}.{unit}-strain_gene.tsv",
        filt="{stem}/midasdb_{dbv}.gene{centroid}.strain_meta-complete90-contam5-pos0.tsv",
    params:
        pseudo=0,
    shell:
        "{input.script} {input.gene} {input.filt} {params.pseudo} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_spgc_genomes:
    output:
        "{stemA}.spgc{stemB}.{unit}-strain_gene.prevalence.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stemA}.spgc{stemB}.{unit}-strain_gene.tsv",
        filt="{stemA}.spgc{stemB}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    params:
        pseudo=0,
    shell:
        "{input.script} {input.gene} {input.filt} {params.pseudo} {output}"


rule count_pangenome_fractions_across_genomes:
    output:
        "{stem}.{unit}-strain_gene.prevalence_class_fraction.tsv",
    input:
        script="scripts/count_pangenome_fractions_across_genomes.py",
        prevalence="{stem}.{unit}-strain_gene.prevalence.tsv",
        gene="{stem}.{unit}-strain_gene.tsv",
    params:
        core_vs_shell=0.9,
        shell_vs_cloud=0.15,
    shell:
        "{input.script} {input.prevalence} {input.gene} {params.core_vs_shell} {params.shell_vs_cloud} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
# SPGC and Ref genome pdist based on filtered gene content dissimilarity (jaccard).
# TODO: Drop unecessary parts of the filenames. Only stems that are important are species and centroidB
rule compute_reference_and_spgc_pairwise_filtered_gene_content_dissimilarities:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.genefilt-ratio{ratio}_{metric}_pdist.pkl",
    input:
        script="scripts/compute_reference_and_spgc_pairwise_filtered_gene_content_dissimilarities.py",
        spgc_gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.tsv",
        ref_gene="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.{unit}-strain_gene.tsv",
        ref_filt="data/species/sp-{species}/midasdb_{dbv}.gene75.strain_meta-complete90-contam5-pos0.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    params:
        # NOTE: 0 is a special "ratio" where no filtering will be performed.
        maximum_prevalence_ratio=lambda w: int(w.ratio),
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.spgc_gene} {input.ref_gene} {input.spgc_filt} {input.ref_filt} {params.maximum_prevalence_ratio} {wildcards.metric} {output}"


# TODO: Add this dissimilarity to the analysis metadata table.
rule compute_reference_and_spgc_pairwise_batch_corrected_gene_content_dissimilarities:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
    input:
        script="scripts/compute_reference_and_spgc_pairwise_batch_corrected_gene_content_dissimilarities.py",
        spgc_gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.tsv",
        ref_gene="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.{unit}-strain_gene.tsv",  # FIXME: Rename this file.
        ref_filt="data/species/sp-{species}/midasdb_{dbv}.gene75.strain_meta-complete90-contam5-pos0.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    conda:
        "conda/toolz4.yaml"
    resources:
        walltime_hr=12,
    shell:
        "{input.script} {input.spgc_gene} {input.ref_gene} {input.spgc_filt} {input.ref_filt} cosine {output}"


rule cluster_genes_based_on_cooccurence_in_spgc_strains:
    output:
        "data/group/{stem}.uhgg-strain_gene.gene_clust-t10.tsv",
    input:
        script="scripts/cluster_genes_based_on_cooccurence.py",
        gene="data/group/{stem}.uhgg-strain_gene.tsv",
        filt="data/group/{stem}.strain_meta-s95-d100-a0-pos100-std25.tsv",
    params:
        thresh=10 / 100,
    conda:
        "conda/toolz2.yaml"
    shell:
        "{input.script} {input.gene} {input.filt} {params.thresh} {output}"


rule calculate_morans_i_for_ref_strains:
    output:
        "data/species/sp-{species}/midasdb.gene75_{dbv}.uhgg-strain_gene.morans_i.tsv",
    input:
        script="scripts/calculate_morans_i.py",
        gene="data/species/sp-{species}/midasdb.gene75_{dbv}.uhgg-strain_gene.tsv",
        filt="data/species/sp-{species}/midasdb.gene75_{dbv}.strain_meta-complete90-contam5-pos0.tsv",
        pdist="data/species/sp-{species}/midasdb_{dbv}.gtpro.geno_pdist-mask10-pseudo10.pkl",
    shell:
        "{input.script} {input.gene} {input.filt} {input.pdist} {output}"


rule calculate_morans_i_for_spgc_strains:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.morans_i.tsv",
    input:
        script="scripts/calculate_morans_i.py",
        gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s95-d100-a0-pos100-std25.tsv",
        pdist="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.pkl",
    shell:
        "{input.script} {input.gene} {input.filt} {input.pdist} {output}"


rule construct_gene_x_cog_category_matrix_v15:
    output:
        cog_matrix="data/species/sp-{species}/midasdb_v15.gene{centroid}_x_cog_category.tsv",
    input:
        script="scripts/construct_gene_x_cog_category_matrix.py",
        emapper="ref/midasdb_uhgg_v15/pangenomes/{species}/eggnog.tsv",
        clustering="ref/midasdb_uhgg_v15/pangenomes/{species}/gene_info.txt",
        cog_category="ref/cog-20.meta.tsv",
    params:
        agg=lambda w: {
            99: "centroid_99",
            95: "centroid_95",
            90: "centroid_90",
            85: "centroid_85",
            80: "centroid_80",
            75: "centroid_75",
        }[int(w.centroid)],
    shell:
        "{input.script} {input.emapper} {input.clustering} {input.cog_category} {params.agg} {output.cog_matrix}"


rule construct_gene_x_cog_category_matrix_v20:
    output:
        cog_matrix="data/species/sp-{species}/midasdb_v20.gene{centroid}_x_cog_category.tsv",
    input:
        script="scripts/construct_gene_x_cog_category_matrix.py",
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
        clustering="ref/midasdb_uhgg_v20/pangenomes/{species}/gene_info.txt",
        cog_category="ref/cog-20.meta.tsv",
    params:
        agg=lambda w: {
            99: "centroid_99",
            95: "centroid_95",
            90: "centroid_90",
            85: "centroid_85",
            80: "centroid_80",
            75: "centroid_75",
        }[int(w.centroid)],
    shell:
        "{input.script} {input.emapper} {input.clustering} {input.cog_category} {params.agg} {output.cog_matrix}"


rule compile_gene_metadata_v15:
    output:
        meta="data/species/sp-{species}/midasdb_v15.gene{centroid}_meta.tsv",
        cog_matrix="data/species/sp-{species}/midasdb_v15.gene{centroid}_x_cog_category.tsv",
    input:
        script="scripts/compile_midasdb_gene_metadata.py",
        # FIXME: Add agg parameter to the script to take votes on cog_category.
        emapper="ref/midasdb_uhgg_v15/pangenomes/{species}/eggnog.tsv",
        clustering="ref/midasdb_uhgg_v15/pangenomes/{species}/gene_info.txt",
        nlength="ref/midasdb_uhgg_v15/pangenomes/{species}/genes.len",
        cog_category="ref/cog-20.meta.tsv",
    wildcard_constraints:
        centroid="99|95|90|85|80|75",
    # params:
    # agg=lambda w: {
    #     99: "centroid_99",
    #     95: "centroid_95",
    #     90: "centroid_90",
    #     85: "centroid_85",
    #     80: "centroid_80",
    #     75: "centroid_75",
    # }[int(w.centroid)],
    shell:
        "{input.script} {input.emapper} {input.clustering} {input.cog_category} {input.nlength} {output.meta} {output.cog_matrix}"

rule compile_gene_metadata_v20:
    output:
        meta="data/species/sp-{species}/midasdb_v20.gene{centroid}_meta.tsv",
        cog_matrix="data/species/sp-{species}/midasdb_v20.gene{centroid}_x_cog_category.tsv",
    input:
        script="scripts/compile_midasdb_gene_metadata_v20.py",
        # FIXME: Add agg parameter to the script to take votes on cog_category.
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
        genes_info="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
        cog_category="ref/cog-20.meta.tsv",
    wildcard_constraints:
        centroid="99|95|90|85|80|75",
    shell:
        "{input.script} {input.emapper} {input.genes_info} {input.cog_category} {output.meta} {output.cog_matrix}"


rule perform_ibd_association_test_with_hmp2_strains:
    output:
        "data/group/hmp2/species/sp-{species}/{stem}.gtpro.sfacts-fit.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc-fit.{unit}-strain_gene.hmp2_mwas-f{frac_thresh}-n{num_thresh}.tsv",
    input:
        script="scripts/hmp2_mwas.py",
        comm="data/group/hmp2/species/sp-{species}/{stem}.gtpro.sfacts-fit.comm.tsv",
        gene="data/group/hmp2/species/sp-{species}/{stem}.gtpro.sfacts-fit.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc-fit.{unit}-strain_gene.tsv",
        filt="data/group/hmp2/species/sp-{species}/{stem}.gtpro.sfacts-fit.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc-fit.strain_meta-s95-d100-a0-pos100-std25.tsv",
        mgen="meta/hmp2/mgen.tsv",
        preparation="meta/hmp2/preparation.tsv",
        stool="meta/hmp2/stool.tsv",
        subject="meta/hmp2/subject.tsv",
    params:
        frac_thresh=lambda w: int(w.frac_thresh) / 100,
        num_thresh=lambda w: int(w.num_thresh),
    shell:
        "{input.script} {input.comm} {input.gene} {input.filt} {input.mgen} {input.preparation} {input.stool} {input.subject} {params.frac_thresh} {params.num_thresh} {output}"


# NOTE: This rule takes the super long filename and turns it into a much shorter one for, e.g., notebooks.
rule alias_sfacts_outputs:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.{stemB}",
    input:
        source=lambda w: (
            "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.{w.stemB}".format(
                w=w,
                sfacts_params=get_sfacts_stem(config, w.species, w.group),
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_sfacts_outputs,


# NOTE: This ruleorder section is a place to clear up ambiguity about whether the pipeline should be run on the full sfacts spec (yes)
# or on previously aliased files (no).
# When I get errors like this:
#
# AmbiguousRuleException:
# Rules export_sfacts_comm and alias_sfacts_outputs are ambiguous for the file data/group/hmp2/species/sp-102506/r.proc.gtpro.sfacts-fit.comm.tsv.
# Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
#
# I can just add that rule to the end of the list.
ruleorder: alias_sfacts_outputs > export_sfacts_comm > identify_strain_samples > aggregate_strain_metagenotype > compute_reference_and_spgc_pairwise_genotype_masked_hamming_distance > match_strains_to_genomes_based_on_genotype


# That forces the aliasing to the end of the computation.


# NOTE: This rule takes the super long filename and turns it into a much shorter one for, e.g., notebooks.
rule alias_spgc_analysis_outputs:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.gene99_{dbv}-{btv}-agg75.spgc-fit.{stemB}",
    input:
        source=lambda w: (
            "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.gene99_{w.dbv}-{w.btv}-agg75.spgc_{spgc_params}.{w.stemB}".format(
                w=w,
                spgc_params=get_spgc_stem(config, w.species, w.group),
                sfacts_params=get_sfacts_stem(config, w.species, w.group),
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_spgc_analysis_outputs,


# NOTE: This ruleorder section is a place to clear up ambiguity about whether the pipeline should be run on the full sfacts and spgc spec (yes)
# or on previously aliased files (no).
# When I get errors like this:
#
# > AmbiguousRuleException: Rules count_pangenome_fractions_across_genomes and
# alias_spgc_analysis_outputs are ambiguous for the file
# data/group/xjin_ucfmt_hmp2/species/sp-100022/r.proc.gtpro.sfacts-fit.gene99_new-v22-agg75.spgc-fit.uhgg-
# strain_gene.prevalence_class_fraction-hmp2.tsv.
#
# I can just add that rule to the end of the list.
ruleorder: alias_spgc_analysis_outputs > calculate_gene_prevalence_in_spgc_genomes > cluster_genes_based_on_cooccurence_in_spgc_strains > aggregate_uhgg_strain_gene_by_annotation > count_pangenome_fractions_across_genomes > match_strains_to_genomes_based_on_genotype


rule collect_analysis_files:
    output:
        "data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",
    input:
        # Gene content
        gene="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.tsv",
        eggnog_strain_gene="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.eggnog-strain_gene.tsv",
        comm="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.comm.tsv",
        spgc="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.nc",
        # Genotype
        geno="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.spgc_ss-all.mgtp.nc",
        # Strain analysis
        meta="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.strain_meta_spgc_and_ref.tsv",
        diss_stats="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.strain_diss_spgc_and_ref.tsv",  # Distance to nearest ref
        geno_diss_pdmat="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.spgc_ss-all.geno_uhgg-v20_pdist-mask10-pseudo10.pkl",
        uhgg_diss_pdmat="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
        eggnog_diss_pdmat="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.eggnog-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
        strain_samples="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.spgc_ss-all.strain_samples.tsv",
        # Gene analysis
        gene_clust="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.gene_clust-t10.tsv",
        # morans_i="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.morans_i.tsv",
        cog_cat_matrix="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_cog_category.tsv",  # Whats the difference between with and without *.emapper.gene75.*?
        cog_meta="ref/cog-20.meta.tsv",
        prevalence="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.prevalence.tsv",
        eggnog_prevalence="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.eggnog-strain_gene.prevalence.tsv",
        prev_class_frac="data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.uhgg-strain_gene.prevalence_class_fraction.tsv",
        species_depth="data/group/{group}/species/sp-{species}/r.proc.gene99_v20-v23-agg75.spgc_specgene-ref-filt-p95.species_depth.tsv",  # Does not require SPGC.
        # Gene metadata
        ref_eggnog_prevalence="data/species/sp-{species}/midasdb.gene75_v20.eggnog-strain_gene.prevalence.tsv",
        eggnog="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_eggnog.tsv",
        cog="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_cog.tsv",
        ko="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_ko.tsv",
        go="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_go.tsv",
        kegg_module="data/species/sp-{species}/midasdb_v20.emapper.gene75_x_kegg_module.tsv",
        amr="data/species/sp-{species}/midasdb_v20.gene75_x_amr.tsv",
        plasmid="data/species/sp-{species}/midasdb_v20.gene75_x_genomad_plasmid.tsv",
        phage="data/species/sp-{species}/midasdb_v20.gene75_x_genomad_virus.tsv",
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
        gene_meta="data/species/sp-{species}/midasdb_v20.gene75_meta.tsv",
    shell:
        "echo {input} > {output}"


localrules:
    collect_analysis_files,


rule collect_spgc_manuscript_analysis_files:
    input:
        "data/group/xjin/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.BENCHMARK_GRID.flag",
        "data/group/hmp2/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag.SELECT_SPECIES.flag",
        "data/group/ucfmt/species/sp-102506/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",


localrules:
    collect_spgc_manuscript_analysis_files,


rule collect_ucfmt_analysis_files:
    input:
        "data/group/ucfmt/species/sp-102506/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",
        "data/group/ucfmt/species/sp-101337/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",
        "data/group/ucfmt/species/sp-102478/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",
        "data/group/ucfmt/species/sp-101433/r.proc.gtpro.sfacts-fit.gene99_v20-v23-agg75.spgc-fit.all_analysis_files.flag",


localrules:
    collect_ucfmt_analysis_files,
