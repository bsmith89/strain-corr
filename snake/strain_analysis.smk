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
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.mgtp.nc",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        # NOTE: species_gene.list and species_depth.tsv are both extracted from the OUTPUTs of SPGC so that I don't need to pattern-match all the filename parameters.
        # See "extract_*" rules in the strain_reconstruction file.
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.species_gene.list",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.species_depth.tsv",
        strain_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.uhgg-strain_gene.tsv",
    params:
        ambig_thresh=0.1,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.species_gene} {input.species_depth} {input.spgc_agg_mgtp} {input.strain_samples} {input.strain_gene} {params.ambig_thresh} {output}"


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
        "data/group/{group}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/{group}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        sample_to_strain="data/group/{group}/{stem}.spgc_ss-{ss}.strain_samples.tsv",
    params:
        min_species_genes_frac=90 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=1 / 1000,
        min_geno_positions=100,
    shell:
        "{input.script} {input.meta} {input.sample_to_strain} DROP_THIS_ARG_TODO {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {params.min_geno_positions} {output}"


# TODO: Move this into a new snakefile with all of the reference database
# work.
rule ref_gene_copy_number_to_presence_table:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
    input:
        script="scripts/gene_copy_number_nc_to_strain_gene_tsv.py",
        # TODO: Change "copies" input file naming to `midasdb.gene{centroid}_{dbv}.reference_copy_number.nc`.
        copies="data/species/sp-{species}/gene{centroid}_{dbv}.reference_copy_number.nc",
    shell:
        "{input.script} {input.copies} {output}"


# TODO: Make this rule work for MIDAS2DB v1.0
# NOTE: Legacy of "_new" identifier for MIDAS2DB v1.0
# TODO (2023-09-15): This rule should probably be moved to
# strain_reconstruction, since it's an input to reference-based species-gene
# calling.
rule collect_filtering_metadata_for_uhgg_ref_strains_v1_as_new:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_new.strain_meta-complete90-contam5-pos100.tsv",
    input:
        script="scripts/filter_ref_strains_v10.py",  # TODO
        meta="ref/uhgg_genomes_all_4644.tsv",
        pos="data/species/sp-{species}/midasdb_new.gtpro.geno.npositions.tsv",
        # FIXME: Change 
        genes="data/species/sp-{species}/midasdb.gene{centroid}_new.uhgg-strain_gene.tsv",
    params:
        min_completeness=90 / 100,
        # NOTE (2023-09-19): Some species (e.g.) 102386 have refs with only >2% contam (duplication?);
        # TODO: Increase this threshold to 5%.
        max_contamination=5 / 100,
        min_positions=100,
    shell:
        "{input.script} {input.meta} {input.pos} {input.genes} {wildcards.species} {params.min_completeness} {params.max_contamination} {params.min_positions} {output}"


# NOTE: This rule is a stop-gap as I migrate to MIDAS2DB v2
# TODO (2023-09-15): This rule should probably be moved to
# strain_reconstruction, since it's an input to reference-based species-gene
# calling.
rule collect_filtering_metadata_for_uhgg_ref_strains_v15:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_v15.strain_meta-complete90-contam5-pos100.tsv",
    input:
        script="scripts/filter_ref_strains_v15.py",
        meta="ref/midasdb_uhgg_v15/metadata/2023-11-11-genomes-all_metadata.tsv",
        genome_to_species="ref/midasdb_uhgg_v15/genomes.tsv",
        pos="data/species/sp-{species}/midasdb_v15.gtpro.geno.npositions.tsv",
        genes="data/species/sp-{species}/midasdb.gene{centroid}_v15.uhgg-strain_gene.tsv",
    params:
        min_completeness=90 / 100,
        max_contamination=5 / 100,
        min_positions=100,
    shell:
        "{input.script} {input.meta} {input.genome_to_species} {input.pos} {input.genes} {wildcards.species} {params.min_completeness} {params.max_contamination} {params.min_positions} {output}"


# TODO: Replace this rule with a genotype concatenation rule to achieve the
# same effect (harnessing the generic rule below.)
rule compute_reference_and_spgc_pairwise_genotype_masked_hamming_distance:
    output:
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask{thresh}-pseudo{pseudo}.pkl",
    input:
        script="scripts/calculate_ref_and_spgc_pairwise_genotype_masked_hamming_distance.py",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb_{dbv}.gtpro.mgtp.nc",
    params:
        ambiguity_threshold=lambda w: int(w.thresh) / 100,
        pseudo=lambda w: int(w.pseudo) / 10,
    conda:
        "conda/sfacts.yaml"
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


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
# SPGC and Ref genome dereplication clustering (based on genotype dissimilarity)
rule cocluster_reference_and_spgc_genomes_based_on_genotype_dissimilarity:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-{diss_params}.coclust-{thresh}.tsv",
    input:
        script="scripts/cluster_from_cdmat_pkl.py",
        pickle="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-{diss_params}.pkl",
        # NOTE: I decided I don't want to filter strains before co-clustering.
        # # NOTE: The filtering recipe is only currently implemented for xjin_ucfmt_hmp2.
        # spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
        # ref_filt="data/species/sp-{species}/midasdb.strain_meta-complete90-contam5-pos0.tsv",
    params:
        thresh=lambda w: int(w.thresh) / 1000,
    conda:
        "conda/toolz2.yaml"
    shell:
        "{input.script} {input.pickle} {params.thresh} {output}"


rule combine_all_ref_and_spgc_metadata_and_clustering:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta_for_analysis.tsv",
    input:
        script="scripts/combine_all_ref_and_spgc_metadata_and_clustering.py",
        spgc="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
        ref="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.strain_meta-complete90-contam5-pos100.tsv",
        geno_pdist="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.pkl",
        # NOTE: genefilt-ratio0 means there is no filtering on gene prevalence ratio
        # higher values of R would filter genes with R-times greater or less prevalence in SPGC genomes compared to reference genomes.
        # As a result, these distances are not "batch corrected" with *-ratio0*.
        # uhgg_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.genefilt-ratio0_cosine_pdist.pkl",
        uhgg_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
        eggnog_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.eggnog-strain_gene.gene_batch_corrected_cosine_pdist.pkl",
        # eggnog_pdist="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.eggnog-strain_gene.genefilt-ratio0_cosine_pdist.pkl",
        clust="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.coclust-10.tsv",
    shell:
        "{input.script} {input.spgc} {input.ref} {input.geno_pdist} {input.uhgg_pdist} {input.eggnog_pdist} {input.clust} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_ref_genomes:
    output:
        "{stem}/midasdb.gene75_{dbv}.{unit}-strain_gene.prevalence.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stem}/midasdb.gene75_{dbv}.{unit}-strain_gene.tsv",
        filt="{stem}/midasdb.gene75_{dbv}.strain_meta-complete90-contam5-pos100.tsv",
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
        filt="{stemA}.spgc{stemB}.strain_meta-s90-d100-a1-pos100.tsv",
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
        ref_filt="data/species/sp-{species}/midasdb.gene75_{dbv}.strain_meta-complete90-contam5-pos100.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
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
        ref_gene="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.{unit}-strain_gene.tsv",
        ref_filt="data/species/sp-{species}/midasdb.gene75_{dbv}.strain_meta-complete90-contam5-pos100.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
    conda:
        "conda/toolz4.yaml"
    shell:
        "{input.script} {input.spgc_gene} {input.ref_gene} {input.spgc_filt} {input.ref_filt} cosine {output}"


rule cluster_genes_based_on_cooccurence_in_spgc_strains:
    output:
        "data/group/{stem}.uhgg-strain_gene.gene_clust-t10.tsv",
    input:
        script="scripts/cluster_genes_based_on_cooccurence.py",
        gene="data/group/{stem}.uhgg-strain_gene.tsv",
        filt="data/group/{stem}.strain_meta-s90-d100-a1-pos100.tsv",
    params:
        thresh=10 / 100,
    conda:
        "conda/toolz2.yaml"
    shell:
        "{input.script} {input.gene} {input.filt} {params.thresh} {output}"


rule cluster_genes_based_on_cooccurence_in_ref_strains:
    output:
        "data/species/{stem}.uhgg-strain_gene.gene_clust-t10.tsv",
    input:
        script="scripts/cluster_genes_based_on_cooccurence.py",
        gene="data/species/{stem}.uhgg-strain_gene.tsv",
        filt="data/species/{stem}.strain_meta-complete90-contam5-pos100.tsv",
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
        filt="data/species/sp-{species}/midasdb.gene75_{dbv}.strain_meta-complete90-contam5-pos100.tsv",
        pdist="data/species/sp-{species}/midasdb_{dbv}.gtpro.geno_pdist-mask10-pseudo10.pkl",
    shell:
        "{input.script} {input.gene} {input.filt} {input.pdist} {output}"


rule calculate_morans_i_for_spgc_strains:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.morans_i.tsv",
    input:
        script="scripts/calculate_morans_i.py",
        gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
        pdist="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_uhgg-{dbv}_pdist-mask10-pseudo10.pkl",
    shell:
        "{input.script} {input.gene} {input.filt} {input.pdist} {output}"


rule compile_gene_metadata:
    output:
        meta="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}_meta.tsv",
        cog_matrix="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}_x_cog_category.tsv",
    input:
        script="scripts/compile_midasdb_gene_metadata.py",
        # FIXME: Add agg parameter to the script to take votes on cog_category.
        emapper="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/eggnog.tsv",
        clustering="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
        nlength="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/genes.len",
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
        "{input.script} {input.emapper} {input.clustering} {params.agg} {input.cog_category} {input.nlength} {output.meta} {output.cog_matrix}"


rule perform_ibd_association_test_with_hmp2_strains:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.hmp2_mwas-f{frac_thresh}-n{num_thresh}.tsv",
    input:
        script="scripts/hmp2_mwas.py",
        comm="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.comm.tsv",
        gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_{dbv}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-s90-d100-a1-pos100.tsv",
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
rule alias_spgc_analysis_outputs:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.gene99_{dbv}-v22-agg75.spgc-fit.{stemB}",
    input:
        source=lambda w: (
            "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.gene99_{w.dbv}-v22-agg75.spgc_{spgc_params}.{w.stemB}".format(
                w=w,
                spgc_params=config["species_group_to_spgc_stem"][(w.species, w.group)],
                sfacts_params=config["species_group_to_sfacts_stem"][
                    (w.species, w.group)
                ],
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
ruleorder: alias_spgc_analysis_outputs > calculate_gene_prevalence_in_spgc_genomes > cluster_genes_based_on_cooccurence_in_spgc_strains > aggregate_uhgg_strain_gene_by_annotation


rule collect_analysis_files:
    output: "data/group/{group}/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v15-v22-agg75.spgc-fit.uhgg-strain_gene.all_analysis_files.flag"
    input:
        gene="data/group/hmp2/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v15-v22-agg75.spgc-fit.uhgg-strain_gene.tsv",
        meta="data/group/hmp2/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v15-v22-agg75.spgc-fit.strain_meta_for_analysis.tsv",
        gene_clust="data/group/hmp2/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v15-v22-agg75.spgc-fit.uhgg-strain_gene.gene_clust-t10.tsv",
        morans_i="data/group/hmp2/species/sp-{species}/r.proc.gtpro.sfacts-fit.gene99_v15-v22-agg75.spgc-fit.uhgg-strain_gene.morans_i.tsv",
        cog_cat_matrix="data/species/sp-{species}/midasdb_v15.gene75_x_cog_category.tsv",
    shell:
        "echo {input} > {output}"
