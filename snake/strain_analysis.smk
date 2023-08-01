rule compile_strain_spgc_metadata:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_meta.tsv",
    input:
        script="scripts/compile_spgc_results_metadata.py",
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_gene.list",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{gene_params}.spgc_specgene-{specgene_params}.species_depth.tsv",
        strain_fit="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        strain_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{gene_params}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{t}_thresh-{thresh_params}.strain_gene.tsv",
    conda:
        "conda/sfacts.yaml"
    shell:
        "{input.script} {input.species_gene} {input.species_depth} {input.strain_fit} {input.strain_samples} {input.strain_gene} {output}"


rule compile_spgc_to_ref_strain_report:
    output:
        gene_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.gene_stats.tsv",
        spgc_strain_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_strain_stats.tsv",
        ref_strain_stats="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.ref_strain_stats.tsv",
        html="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_ref_comparison.html",
    input:
        nb="nb/prototype_compare_spgc_and_refs.ipynb",
        species_taxonomy="ref/gtpro/species_taxonomy_ext.tsv",
        sample_to_spgc="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.strain_samples.tsv",
        sfacts_fit="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.world.nc",
        ref_geno="data/species/sp-{species}/gtpro_ref.mgtp.nc",
        spgc_meta="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        ref_gene_copy_number_uhgg="ref/midasdb_uhgg_pangenomes/{species}/gene{centroidB}.reference_copy_number.nc",
        spgc_gene_uhgg="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_gene.tsv",
        spgc_gene_uhgg_depth="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_depth_ratio.tsv",
        spgc_gene_uhgg_corr="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_correlation.tsv",
        uhgg_x_eggnog="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_eggnog.tsv",
        uhgg_x_top_eggnog="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_top_eggnog.tsv",
        uhgg_gene_length="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
        gene_annotations="data/species/sp-{species}/pangenome.centroids.emapper.d/proteins.emapper.annotations",
        uhgg_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}-{pang}-agg{centroidB}.depth2.nc",
        species_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}.species_depth.tsv",
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
        extra_args=lambda w: [f"-r species_id {w.species}"],
    log:
        nb="nb/papermill/data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.spgc_ref_comparison.ipynb",
    conda:
        "conda/papermill.yaml"
    threads: 1
    resources:
        mem_mb=30_000,
        pmem=30_000,
        walltime_hr=12,
    shell:
        """
        papermill {input.nb} {log.nb} {params.extra_args} {params.input_path_args} {params.output_path_args}
        jupyter nbconvert --to html --embed-images {log.nb} --stdout > {output.html}
        """


use rule compile_spgc_to_ref_strain_report as compile_ucfmt_donor_strain_report with:
    output:
        html="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.donor_comparison.html",
    input:
        nb="nb/prototype_compare_ucfmt_donor_strains.ipynb",
        species_taxonomy="ref/gtpro/species_taxonomy_ext.tsv",
        sample_to_spgc="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.spgc_ss-{ss}.strain_samples.tsv",
        sfacts_fit="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.world.nc",
        ref_geno="data/species/sp-{species}/gtpro_ref.mgtp.nc",
        spgc_meta="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_meta.tsv",
        ref_gene_copy_number_uhgg="ref/midasdb_uhgg_pangenomes/{species}/gene{centroidB}.reference_copy_number.nc",
        spgc_gene_uhgg="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.strain_gene.tsv",
        spgc_gene_uhgg_depth="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_depth_ratio.tsv",
        spgc_gene_uhgg_corr="data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}.strain_correlation.tsv",
        uhgg_x_eggnog="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_eggnog.tsv",
        uhgg_x_top_eggnog="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_top_eggnog.tsv",
        uhgg_gene_length="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
        gene_annotations="data/species/sp-{species}/pangenome.centroids.emapper.d/proteins.emapper.annotations",
        uhgg_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}-{pang}-agg{centroidB}.depth2.nc",
        species_depth="data/group/{group}/species/sp-{species}/r.proc.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}.species_depth.tsv",
        mgen="meta/ucfmt/mgen.tsv",
        sample="meta/ucfmt/sample.tsv",
        subject="meta/ucfmt/subject.tsv",
    log:
        nb="nb/papermill/data/group/{group}/species/sp-{species}/{stem}.gtpro.{fit}.gene{centroidA}-{pang}-agg{centroidB}.spgc_specgene-{specgene}_ss-{ss}_t-{t}_thresh-{thresh}.donor_comparison.ipynb",
