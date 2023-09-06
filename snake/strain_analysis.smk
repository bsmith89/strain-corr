rule compile_spgc_strain_metadata:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_meta.tsv",
    input:
        script="scripts/compile_spgc_results_metadata.py",
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_gene.list",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_depth.tsv",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_mgtp.nc",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        strain_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_gene.tsv",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.species_gene} {input.species_depth} {input.spgc_agg_mgtp} {input.strain_samples} {input.strain_gene} {output}"


rule filter_hmp2_spgc_strains:
    output:
        "data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_filt-s90-d100-a1.tsv",
    input:
        script="scripts/filter_spgc_strains.py",
        meta="data/group/xjin_ucfmt_hmp2/{stem}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        sample_to_strain="data/group/xjin_ucfmt_hmp2/{stem}.spgc_ss-{ss}.strain_samples.tsv",
        mgen="meta/hmp2/mgen.tsv",
    params:
        min_species_genes_frac=90 / 100,
        min_total_depth=100 / 100,
        gene_count_outlier_alpha=1 / 1000,
    shell:
        "{input.script} {input.meta} {input.sample_to_strain} <(cut -f1 {input.mgen}) {params.min_species_genes_frac} {params.min_total_depth} {params.gene_count_outlier_alpha} {output}"


rule compute_reference_and_spgc_pairwise_genotype_dissimilarities:
    output:
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.geno_pdist.pkl",
    input:
        script="scripts/calculate_ref_and_spgc_pairwise_genotype_dissimilarity.py",
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.strain_mgtp.nc",
        ref_geno="data/species/sp-{species}/midasdb.geno.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.spgc_agg_mgtp} {input.ref_geno} {output}"


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
        spgc_agg_mgtp="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.strain_mgtp.nc",
        # ref_geno="data/species/sp-{species}/gtpro_ref.mgtp.nc",  # Original
        ref_geno="data/species/sp-{species}/midasdb.geno.nc",  # Re-calculated
        spgc_meta="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        ref_gene_copy_number_uhgg="data/species/sp-{species}/gene{centroidB}_new.reference_copy_number.nc",
        spgc_gene_uhgg="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_gene.tsv",
        spgc_gene_uhgg_depth="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_depth_ratio.tsv",
        spgc_gene_uhgg_corr="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}_new-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_correlation.tsv",
        uhgg_x_eggnog="data/species/sp-{species}/pangenome_new.centroids.emapper.gene_x_eggnog.tsv",
        uhgg_x_top_eggnog="data/species/sp-{species}/pangenome_new.centroids.emapper.gene_x_top_eggnog.tsv",
        uhgg_gene_length="ref/midasdb_uhgg_new/pangenomes/{species}/cluster_info.txt",
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
