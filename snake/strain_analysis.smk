rule download_cog_definitions:
    output: "ref/cog-20.meta.tsv"
    params:
        url="https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab"
    shell:
        curl_recipe

rule download_cog_categories:
    output: "ref/cog-20.categories.tsv"
    params:
        url="https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab"
    shell:
        curl_recipe


rule compile_spgc_strain_metadata:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_meta.tsv",
    input:
        script="scripts/compile_spgc_results_metadata.py",
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_gene.list",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_depth.tsv",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.mgtp.nc",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        strain_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.uhgg-strain_gene.tsv",
    params:
        ambig_thresh=0.1
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.species_gene} {input.species_depth} {input.spgc_agg_mgtp} {input.strain_samples} {input.strain_gene} {params.ambig_thresh} {output}"


rule count_ref_geno_positions:
    output:
        "data/species/sp-{species}/midasdb.geno.npositions.tsv",
    input:
        script="scripts/count_ref_geno_positions.py",
        geno="data/species/sp-{species}/midasdb.mgtp.nc",
    params:
        ambig_thresh=0.01,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.geno} {params.ambig_thresh} {output}"


rule collect_filtering_metadata_for_hmp2_spgc_strains:
    output:
        "data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        sample_to_strain="data/group/xjin_ucfmt_hmp2/{stem}.spgc_ss-{ss}.strain_samples.tsv",
        mgen="meta/hmp2/mgen.tsv",
    params:
        min_species_genes_frac=90 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=1 / 1000,
        min_geno_positions=100,
    shell:
        "{input.script} {input.meta} {input.sample_to_strain} <(cut -f1 {input.mgen}) {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {params.min_geno_positions} {output}"


rule ref_gene_copy_number_to_presence_table:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_new.uhgg-strain_gene.tsv",
    input:
        script="scripts/gene_copy_number_nc_to_strain_gene_tsv.py",
        # TODO: Change "copies" input file naming to `midasdb.gene{centroid}_new.reference_copy_number.nc`.
        copies="data/species/sp-{species}/gene{centroid}_new.reference_copy_number.nc",
    shell:
        "{input.script} {input.copies} {output}"


# TODO (2023-09-15): This rule should probably be moved to
# strain_reconstruction, since it's an input to reference-based species-gene
# calling.
rule collect_filtering_metadata_for_uhgg_ref_strains:
    output:
        "data/species/sp-{species}/midasdb.gene{centroid}_new.strain_meta-complete90-contam5-pos100.tsv",
    input:
        script="scripts/filter_ref_strains.py",  # TODO
        meta="ref/uhgg_genomes_all_4644.tsv",
        pos="data/species/sp-{species}/midasdb.geno.npositions.tsv",
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


rule compute_reference_and_spgc_pairwise_genotype_masked_hamming_distance:
    output:
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_pdist-mask{thresh}-pseudo{pseudo}.pkl",
    input:
        script="scripts/calculate_ref_and_spgc_pairwise_genotype_masked_hamming_distance.py",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb.mgtp.nc",
    params:
        ambiguity_threshold=lambda w: int(w.thresh) / 100,
        pseudo=lambda w: int(w.pseudo) / 10,
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.spgc_agg_mgtp} {input.ref_geno} {params.ambiguity_threshold} {params.pseudo} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
# SPGC and Ref genome dereplication clustering (based on genotype dissimilarity)
rule cocluster_reference_and_spgc_genomes_based_on_genotype_dissimilarity:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_pdist-{diss_params}.coclust-10.tsv",
    input:
        script="scripts/cluster_from_cdmat_pkl.py",
        pickle="data/group/{group}/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_pdist-{diss_params}.pkl",
        # NOTE: I decided I don't want to filter strains before co-clustering.
        # # NOTE: The filtering recipe is only currently implemented for xjin_ucfmt_hmp2.
        # spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
        # ref_filt="data/species/sp-{species}/midasdb.strain_meta-complete90-contam5-pos0.tsv",
    params:
        thresh=10 / 1000,
    conda:
        "conda/toolz2.yaml"
    shell:
        "{input.script} {input.pickle} {params.thresh} {output}"


rule combine_all_ref_and_spgc_metadata_and_clustering:
    output:
        "data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.strain_meta_for_analysis.tsv",
    input:
        script="scripts/combine_all_ref_and_spgc_metadata_and_clustering.py",
        spgc="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
        ref="data/species/sp-{species}/midasdb.gene{centroidB}_new.strain_meta-complete90-contam5-pos100.tsv",
        geno_pdist="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_pdist-mask10-pseudo10.pkl",
        # TODO: Consider other prevalence thresholds besides 5x? 0x means no filtering.
        gene_raw_pdist="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.genefilt-ratio0_jaccard_pdist.pkl",
        gene_filt_pdist="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.genefilt-ratio10_jaccard_pdist.pkl",
        clust="data/group/xjin_ucfmt_hmp2/species/sp-{species}/{stem}.spgc_ss-{ss}.geno_pdist-mask10-pseudo10.coclust-10.tsv",
    shell:
        "{input.script} {input.spgc} {input.ref} {input.geno_pdist} {input.gene_raw_pdist} {input.gene_filt_pdist} {input.clust} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_ref_genomes:
    output:
        "{stem}/midasdb.gene75_new.{unit}-strain_gene.prevalence.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stem}/midasdb.gene75_new.{unit}-strain_gene.tsv",
        filt="{stem}/midasdb.gene75_new.strain_meta-complete90-contam5-pos100.tsv",
    params:
        pseudo=0,
    shell:
        "{input.script} {input.gene} {input.filt} {params.pseudo} {output}"


# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
rule calculate_gene_prevalence_in_spgc_genomes:
    output:
        "{stemA}.spgc{stemB}.{unit}-strain_gene.prevalence-hmp2.tsv",
    input:
        script="scripts/strain_gene_to_prevalence.py",
        gene="{stemA}.spgc{stemB}.{unit}-strain_gene.tsv",
        filt="{stemA}.spgc{stemB}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
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


use rule count_pangenome_fractions_across_genomes as count_pangenome_fractions_across_hmp2_genomes with:
    output:
        "{stemA}.spgc{stemB}.{unit}-strain_gene.prevalence_class_fraction-hmp2.tsv",
    input:
        script="scripts/count_pangenome_fractions_across_genomes.py",
        prevalence="{stemA}.spgc{stemB}.{unit}-strain_gene.prevalence-hmp2.tsv",
        gene="{stemA}.spgc{stemB}.{unit}-strain_gene.tsv",



# NOTE: Split from `compile_spgc_to_ref_strain_report_new`:
# SPGC and Ref genome pdist based on filtered gene content dissimilarity (jaccard).
# TODO: Drop unecessary parts of the filenames. Only stems that are important are species and centroidB
rule compute_reference_and_spgc_pairwise_filtered_gene_content_dissimilarities:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.genefilt-ratio{ratio}_{metric}_pdist.pkl",
    input:
        script="scripts/compute_reference_and_spgc_pairwise_filtered_gene_content_dissimilarities.py",
        spgc_gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.{unit}-strain_gene.tsv",
        ref_gene="data/species/sp-{species}/midasdb.gene{centroidB}_new.{unit}-strain_gene.tsv",
        ref_filt="data/species/sp-{species}/midasdb.gene75_new.strain_meta-complete90-contam5-pos100.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv"
    params:
        # NOTE: 0 is a special "ratio" where no filtering will be performed.
        maximum_prevalence_ratio=lambda w: int(w.ratio),
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.spgc_gene} {input.ref_gene} {input.spgc_filt} {input.ref_filt} {params.maximum_prevalence_ratio} {wildcards.metric} {output}"


rule cluster_genes_based_on_cooccurence_in_spgc_strains:
    output:
        "data/group/{stem}.uhgg-strain_gene.gene_clust-t10.tsv",
    input:
        script="scripts/cluster_genes_based_on_cooccurence.py",
        gene="data/group/{stem}.uhgg-strain_gene.tsv",
        filt="data/group/{stem}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
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


rule calculate_morans_i_for_both_ref_and_spgc_strains:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.morans_i.tsv",
    input:
        script="scripts/calculate_morans_i.py",
        ref_gene="data/species/sp-{species}/midasdb.gene75_new.uhgg-strain_gene.tsv",
        ref_filt="data/species/sp-{species}/midasdb.gene75_new.strain_meta-complete90-contam5-pos100.tsv",
        spgc_gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        spgc_filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
        pdist="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_pdist-mask10-pseudo10.pkl",
    shell:
        "{input.script} {input.spgc_gene} {input.spgc_filt} {input.ref_gene} {input.ref_filt} {input.pdist} {output}"


rule compile_gene_metadata:
    output:
        meta="data/species/sp-{species}/midasdb.gene_meta.tsv",
        cog_matrix="data/species/sp-{species}/midasdb.gene_x_cog_category.tsv",
        # cog="data/species/sp-{species}/pangenome_new.centroids.emapper.gene_x_eggnog.tsv"  # NOTE: This file is already created and maps to eggnog IDs.
    input:
        script="scripts/compile_midasdb_gene_metadata.py",
        emapper="data/species/sp-{species}/pangenome_new.centroids.emapper.d/proteins.emapper.annotations",
        nlength="ref/midasdb_uhgg_new/pangenomes/{species}/cluster_info.txt",
        cog_category="ref/cog-20.meta.tsv",
    shell:
        "{input.script} {input.emapper} {input.cog_category} {input.nlength} {output.meta} {output.cog_matrix}"


rule perform_ibd_association_test_with_hmp2_strains:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.hmp2_mwas-f{frac_thresh}-n{num_thresh}.tsv",
    input:
        script="scripts/hmp2_mwas.py",
        comm="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.comm.tsv",
        gene="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        filt="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-hmp2-s90-d100-a1-pos100.tsv",
        mgen="meta/hmp2/mgen.tsv",
        preparation="meta/hmp2/preparation.tsv",
        stool="meta/hmp2/stool.tsv",
        subject="meta/hmp2/subject.tsv",
    params:
        frac_thresh=lambda w: int(w.frac_thresh) / 100,
        num_thresh=lambda w: int(w.num_thresh),
    shell:
        "{input.script} {input.comm} {input.gene} {input.filt} {input.mgen} {input.preparation} {input.stool} {input.subject} {params.frac_thresh} {params.num_thresh} {output}"


rule collect_filtering_metadata_for_ucfmt_spgc_strains:
    output:
        "data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta-ucfmt-s90-d100-a1-pos100.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        sample_to_strain="data/group/xjin_ucfmt_hmp2/{stem}.spgc_ss-{ss}.strain_samples.tsv",
        mgen="meta/ucfmt/mgen.tsv",
    params:
        min_species_genes_frac=90 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=1 / 1000,
        min_geno_positions=100,
    shell:
        "{input.script} {input.meta} {input.sample_to_strain} <(cut -f1 {input.mgen}) {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {params.min_geno_positions} {output}"

rule select_dominant_ucfmt_donor_strain_genes:
    output:
        "data/group/xjin_ucfmt_hmp2/species/sp-{species}/r.proc.gtpro.{sfacts}.gene{stemB}.{spgc}.uhgg-strain_gene-ucfmt.tsv",
    input:
        script="scripts/select_dominant_ucfmt_donor_strain_genes.py",
        comm="data/group/xjin_ucfmt_hmp2/species/sp-{species}/r.proc.gtpro.{sfacts}.comm.tsv",
        spgc_filt="data/group/xjin_ucfmt_hmp2/species/sp-{species}/r.proc.gtpro.{sfacts}.gene{stemB}.{spgc}.strain_meta-ucfmt-s90-d100-a1-pos100.tsv",
        mgen_meta="meta/ucfmt/mgen.tsv",
        sample_meta="meta/ucfmt/sample.tsv",
        subject_meta="meta/ucfmt/subject.tsv",
        gene="data/group/xjin_ucfmt_hmp2/species/sp-{species}/r.proc.gtpro.{sfacts}.gene{stemB}.{spgc}.uhgg-strain_gene.tsv",
    params:
        detection_thresh=0.2,
    shell:
        "{input.script} {input.comm} {input.spgc_filt} {params.detection_thresh} {input.mgen_meta} {input.sample_meta} {input.subject_meta} {input.gene} {output}"

# FIXME: This notebook is too intense and slowing down my development. Split out components
# into individual scripts for
# - Count of core/shell/cloud genes (based on both SPGC and Ref definitions) in all genomes spgc/ref.
# - Compiled gene annotation metadata
# - MWAS
# - Moran's I (centered) (SPGC vs. Ref)
# - Gene co-clustering (SPGC vs. Ref)
# - Gene-content dissimilarity (filtered by mismatched prevalence in refs/spgc)
rule compile_spgc_to_ref_strain_report_new:
    output:
        gene_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.gene_stats.tsv",
        gene_meta="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.gene_meta.tsv",
        spgc_strain_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_strain_stats.tsv",
        ref_strain_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.ref_strain_stats.tsv",
        html="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_ref_comparison.html",
    input:
        nb="nb/prototype_compare_spgc_and_refs.ipynb",
        species_taxonomy="ref/gtpro/species_taxonomy_ext.tsv",
        sample_to_spgc="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.strain_samples.tsv",
        sfacts_fit="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.world.nc",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.mgtp.nc",
        # ref_geno="data/species/sp-{species}/gtpro_ref.mgtp.nc",  # Original
        ref_geno="data/species/sp-{species}/midasdb.mgtp.nc",  # Re-calculated
        spgc_meta="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        ref_gene_copy_number_uhgg="data/species/sp-{species}/gene{centroidB}_new.reference_copy_number.nc",
        spgc_gene_uhgg="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.uhgg-strain_gene.tsv",
        spgc_gene_uhgg_depth="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_depth_ratio.tsv",
        spgc_gene_uhgg_corr="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_correlation.tsv",
        uhgg_x_eggnog="data/species/sp-{species}/pangenome_new.centroids.emapper.gene_x_eggnog.tsv",
        uhgg_x_top_eggnog="data/species/sp-{species}/pangenome_new.centroids.emapper.gene_x_top_eggnog.tsv",
        uhgg_gene_length="ref/midasdb_uhgg_new/pangenomes/{species}/cluster_info.txt",
        geno_pdist="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_pdist-mask10-pseudo10.pkl",
        gene_pdist="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.gene_filt5_jaccard_pdist.pkl",
        gene_annotations="data/species/sp-{species}/pangenome_new.centroids.emapper.d/proteins.emapper.annotations",
        uhgg_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}_new-{pang}-agg{centroidB}.depth2.nc",
        species_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}.species_depth.tsv",
        mgen="meta/hmp2/mgen.tsv",
        preparation="meta/hmp2/preparation.tsv",
        stool="meta/hmp2/stool.tsv",
        subject="meta/hmp2/subject.tsv",
    params:
        input_path_args=lambda w, input: [
            f"-p {k}_inpath {v}" for k, v in input.items()
        ],
        output_path_args=lambda w, output: [
            f"-p {k}_outpath {v}" for k, v in output.items()
        ],
        extra_args=lambda w: [
            f"-r species_id {w.species} -r show_unimportant_figures ''"
        ],
    log:
        nb="nb/papermill/data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_ref_comparison.ipynb",
    conda:
        "conda/papermill.yaml"
    threads: 1
    resources:
        mem_mb=10_000,
        pmem=10_000,
        walltime_hr=12,
        papermill_flags="--no-progress-bar",
    shell:
        """
        echo "Running {log.nb}"
        papermill {input.nb} {log.nb} {params.extra_args} {resources.papermill_flags} {params.input_path_args} {params.output_path_args}
        echo "Rendering {output.html}"
        jupyter nbconvert --to html --embed-images {log.nb} --stdout > {output.html}
        """

# NOTE: This rule takes the super long filename and turns it into a much shorter one for benchmarking
rule alias_spgc_analysis_outputs:
    output:
        # "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.gene{w.gene_params}.spgc-fit.{stemB}",
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.sfacts-fit.gene99_new-v22-agg75.spgc-fit.{stemB}",
    input:
        # "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.gene{w.gene_params}.{spgc_params}.{w.stemB}".format(
        source=lambda w: (
            "data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{sfacts_params}.gene99_new-v22-agg75.spgc_{spgc_params}.{w.stemB}".format(
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
ruleorder: alias_spgc_analysis_outputs > calculate_gene_prevalence_in_spgc_genomes > count_pangenome_fractions_across_hmp2_genomes > cluster_genes_based_on_cooccurence_in_spgc_strains > select_dominant_ucfmt_donor_strain_genes > aggregate_uhgg_strain_gene_by_annotation
