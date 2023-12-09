# NOTE: I think (2023-12-06) that this rule (and
# filter_midasdb_all_gene_annotations_by_centroid_new as well) may not actually
# be used for anything any more.  I think I can test this by adding a new input
# file that doesn't exist. Now, if the pipeline requires this file in any way,
# it'll fail.
rule combine_midasdb_all_gene_annotations_new:
    output:
        "data/species/sp-{species}/midasdb_uhgg_{dbv}.gene_annotations.tsv",
    input:
        genome=lambda w: [
            f"ref/midasdb_uhgg_{dbv}/gene_annotations/{w.species}/{genome}/{genome}.tsv.lz4"
            for genome in config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species]
        ],
        # If this rule is a pre-requisite for anything, the whole
        # pipeline will fail:
        fail='__FAIL__',
    params:
        genome_pattern="ref/midasdb_uhgg_{dbv}/gene_annotations/{species}/$genome/$genome.tsv.lz4",
        genome_list=lambda w: config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species],
    shell:
        """
        for genome in {params.genome_list}
        do
            echo -n . >&2
            lz4 -dc {params.genome_pattern} \
                    | awk '$1 != "locus_tag" && $2 != "gene"'
        done > {output}
        echo "" >&2
        """


# FIXME: (2023-12-06) Does it really make sense to filter just to the centroid name?
# Don't we need to aggregate in a more intelligent way?
rule filter_midasdb_all_gene_annotations_by_centroid_new:
    output:
        "data/species/sp-{species}/midasdb_uhgg_{dbv}.gene{centroid}_annotations.tsv",
    input:
        annot="data/species/sp-{species}/midasdb_uhgg_{dbv}.gene_annotations.tsv",
        centroids_list="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
    params:
        col=lambda w: {"99": 2, "95": 3, "90": 4, "85": 5, "80": 6, "75": 7}[w.centroid],
    shell:
        """
        grep -Ff \
            <( \
                lz4cat {input.centroids_list} \
                | cut -f{params.col} \
                | sed '1,1d' \
                | sort \
                | uniq \
                ) \
            {input.annot} \
            > {output}
        """


rule select_species_core_genes_from_reference_new:
    output:
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroid}_{dbv}.spgc_specgene-ref-t{trim_quantile}-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes.py",
        copy_number="data/species/sp-{species}/gene{centroid}_{dbv}.reference_copy_number.nc",
    params:
        trim_quantile=lambda w: float(w.trim_quantile) / 100,
        prevalence=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.copy_number} {params.trim_quantile} {params.prevalence} {output}"


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
        species_genes="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}_{dbv}.spgc_specgene-{specgene}.species_gene.list",
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


rule extract_species_depth_from_spgc_netcdf:
    output:
        "{stemA}.spgc_{spgc_params}.species_depth.tsv",
    input:
        script="scripts/extract_species_depth_from_spgc_netcdf.py",
        ncdf="{stemA}.spgc_{spgc_params}.nc",
    shell:
        "{input.script} {input.ncdf} {output}"

# NOTE: This is really just a stop-gap to make a downstream rule easier to
# parse. Pattern matching gets hard with all these pipeline parameters.
# See rule "compile_spgc_strain_metadata".
rule extract_species_gene_list_from_spgc_netcdf:
    output:
        "{stemA}.spgc_{spgc_params}.species_gene.list",
    input:
        script="scripts/extract_species_gene_list_from_spgc_netcdf.py",
        ncdf="{stemA}.spgc_{spgc_params}.nc",
    shell:
        "{input.script} {input.ncdf} {output}"

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
