
rule select_species_core_genes_from_reference_new:
    output:
        species_gene="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.spgc_specgene-ref-t{trim_quantile}-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes.py",
        copy_number="data/species/sp-{species}/gene{centroid}_{dbv}.reference_copy_number.nc",
    params:
        trim_quantile=lambda w: float(w.trim_quantile) / 100,
        prevalence=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.copy_number} {params.trim_quantile} {params.prevalence} {output}"


# FIXME: Rename a bunch of files to match the */midasdb_v15.gene75.* format.
rule select_species_core_genes_from_reference_by_filtered_set_prevalence:
    output:
        species_gene="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.spgc_specgene-ref-filt-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes2.py",
        prevalence="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.prevalence.tsv",
    params:
        threshold=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.prevalence} {params.threshold} {output}"


# NOTE: (2023-12-01) This rule is the first step in implementing the
# new packaged SPGC pipeline.
rule export_gene_depth_table_from_netcdf:
    output:
        "{stem}.depth2.tsv.gz",
    input:
        script="scripts/export_gene_depth_table_from_netcdf.py",
        depth="{stem}.depth2.nc",
    shell:
        "{input.script} {input.depth} {output}"


# TODO: Drop the ss-all part here and in the partitioning rule.
rule run_spgc:
    output:
        "data/group/{group}/species/sp-{species}/{proc_stem}.gtpro.{sfacts_stem}.gene{centroidA}_{dbv}-{pang_stem}-agg{centroidB}.spgc_specgene-{specgene}_ss-all_t-10_thresh-corr{cthresh}-depth{dthresh}.nc",
    input:
        depth="data/group/{group}/species/sp-{species}/{proc_stem}.gene{centroidA}_{dbv}-{pang_stem}-agg{centroidB}.depth2.tsv.gz",
        partition="data/group/{group}/species/sp-{species}/{proc_stem}.gtpro.{sfacts_stem}.spgc_ss-all.strain_samples.tsv",
        species_genes="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.spgc_specgene-{specgene}.species_gene.list",
    params:
        trim_frac_species_genes=0.15,
        species_free_thresh=1e-4,
        depth_ratio_thresh=lambda w: int(w.dthresh) / 1000,
        corr_thresh=lambda w: int(w.cthresh) / 1000,
    conda:
        "conda/toolz4.yaml"
    shell:
        """
        spgc run --full-output \
                --trim-frac-species-genes {params.trim_frac_species_genes} \
             --species-free-thresh {params.species_free_thresh} \
             --depth-ratio-thresh {params.depth_ratio_thresh} \
             --correlation-thresh {params.corr_thresh} \
             {input.depth} \
             {input.species_genes} \
             {input.partition} \
             {output}
        """


rule extract_strain_gene_hits_from_spgc_netcdf:
    output:
        "{stemA}.spgc_{spgc_params}.uhgg-strain_gene.tsv",
    input:
        script="scripts/extract_strain_gene_hits_from_spgc_netcdf.py",
        ncdf="{stemA}.spgc_{spgc_params}.nc",
    shell:
        "{input.script} {input.ncdf} {output}"


rule calculate_species_depth_directly:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
    input:
        depth="data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_{dbv}-{btp}-agg{centroidB}.depth2.tsv.gz",
        species_genes="data/species/sp-{species}/midasdb.gene{centroidB}_{dbv}.spgc_specgene-{specgene_params}.species_gene.list",
    params:
        trim_frac_species_genes=0.15,
    conda:
        "conda/toolz4.yaml"
    shell:
        """
        spgc estimate_species_depth --trim-frac-species-genes {params.trim_frac_species_genes} {input.depth} {input.species_genes} {output}
        """


# FIXME: Skip the "mofiy_strain_samples_file_format rule by exporting this correctly in the first place.
# NOTE: In this new formulation, I include ALL strain-pure samples, including those
# below the previously considered minimum depth threshold.
# If I want to exclude these samples, I should consider dropping them during
# the correlation and/or depth-ratio calculation steps.
rule identify_strain_samples:
    output:
        "{stem}.spgc_ss-all.strain_samples.tsv",
    input:
        "{stem}.comm.tsv",
    params:
        frac_thresh=0.95,  # NOTE: This must be greater than 50%
    shell:
        """
        awk -v OFS='\t' -v thresh={params.frac_thresh} 'NR > 1 && $3 > thresh {{print $1,$2}}' {input} > {output}
        """


rule aggregate_strain_metagenotype:
    output:
        "{stem}.gtpro.{strain_fit_params}.spgc_ss-all.mgtp.nc",
    input:
        script="scripts/aggregate_strain_metagenotypes_across_strain_samples.py",
        mgtp="{stem}.gtpro.mgtp.nc",
        mapping="{stem}.gtpro.{strain_fit_params}.spgc_ss-all.strain_samples.tsv",
    conda:
        "conda/sfacts.yaml"
    resources:
        mem_mb=10_000,
    shell:
        "{input.script} {input.mapping} {input.mgtp} {output}"
