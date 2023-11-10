

rule combine_midasdb_all_gene_annotations_new:
    output:
        "data/species/sp-{species}/midasdb_uhgg_new.gene_annotations.tsv",
    input:
        genome=lambda w: [
            f"ref/midasdb_uhgg_new/gene_annotations/{w.species}/{genome}/{genome}.tsv.lz4"
            for genome in config["midasdb_uhgg_new_species_genome"][w.species]
        ],
    params:
        genome_pattern="ref/midasdb_uhgg_new/gene_annotations/{species}/$genome/$genome.tsv.lz4",
        genome_list=lambda w: config["midasdb_uhgg_new_species_genome"][w.species],
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




rule filter_midasdb_all_gene_annotations_by_centroid_new:
    output:
        "data/species/sp-{species}/midasdb_uhgg_new.gene{centroid}_annotations.tsv",
    input:
        annot="data/species/sp-{species}/midasdb_uhgg_new.gene_annotations.tsv",
        centroids_list="ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt",
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
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroid}_new.spgc_specgene-ref-t{trim_quantile}-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes.py",
        copy_number="data/species/sp-{species}/gene{centroid}_new.reference_copy_number.nc",
    params:
        trim_quantile=lambda w: float(w.trim_quantile) / 100,
        prevalence=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.copy_number} {params.trim_quantile} {params.prevalence} {output}"


# TODO: Go through and figure out where the "-t{trim_quantile}-" infix shows up.
# NOTE: Unlike select_species_core_genes_from_reference_new, this approach doesn't filter
# down to single copy genes or anything like that. Just looks for high-prevalence genes.
rule select_species_core_genes_from_filtered_reference_new:
    output:
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroid}_new.spgc_specgene-ref2-p{prevalence}.species_gene.list",
    input:
        script="scripts/select_high_prevalence_species_genes2.py",
        prevalence="data/species/sp-{species}/midasdb.gene{centroid}_new.strain_gene.prevalence-pseudo0.tsv",
    params:
        thresh=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.prevalence} {params.thresh} {output}"


rule select_species_core_genes_de_novo:
    output:
        species_gene="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo-n{n_genes}.species_gene.list",
        species_corr="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo-n{n_genes}.species_correlation.tsv",
    input:
        script="scripts/select_highly_correlated_species_genes.py",
        species_depth="{stemA}/{stemB}.gtpro.species_depth.tsv",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
    params:
        n_marker_genes=lambda w: int(w.n_genes),
    shell:
        """
        {input.script} {input.species_depth} {wildcards.species} {input.gene_depth} {params.n_marker_genes} {output.species_gene} {output.species_corr}
        """


# TODO: Use strain-partitioning from a different step to avoid redundant code.
rule select_species_core_genes_de_novo_with_dereplication:
    output:
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_specgene-denovo2-t30-n{n_genes}.species_gene.list",
        species_corr="data/group/{group}/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_specgene-denovo2-t30-n{n_genes}.species_correlation.tsv",
    input:
        script="scripts/select_highly_correlated_species_genes_with_dereplication.py",
        species_depth="data/group/{group}/{stemA}.gtpro.species_depth.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{pangenome_params}.depth2.nc",
    params:
        diss_thresh=0.3,
        trnsf_root=2,
        n_marker_genes=lambda w: int(w.n_genes),
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {wildcards.species} \
                {params.diss_thresh} \
                {input.gene_depth} \
                {params.trnsf_root} \
                {params.n_marker_genes} \
                {output.species_gene} \
                {output.species_corr}
        """




# TODO: Use this in place of midasuhgg* everywhere.
rule alias_species_genes_from_reference_to_match_de_novo_paths_new:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_new-{bowtie_params}-agg{centroidB}.spgc_specgene-ref-{specgene_params}.species_gene.list",
    input:
        "data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}_new.spgc_specgene-ref-{specgene_params}.species_gene.list",
    wildcard_constraints:
        centroidA="99|95|90|85|80|75",
        centroidB="99|95|90|85|80|75",
    shell:
        alias_recipe


# TODO: Use this in place of midasuhgg* everywhere.
rule alias_species_genes_from_filtered_reference_to_match_de_novo_paths_new:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.gene{centroidA}_new-{bowtie_params}-agg{centroidB}.spgc_specgene-ref2-{specgene_params}.species_gene.list",
    input:
        "data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}_new.spgc_specgene-ref2-{specgene_params}.species_gene.list",
    wildcard_constraints:
        centroidA="99|95|90|85|80|75",
        centroidB="99|95|90|85|80|75",
    shell:
        alias_recipe


localrules:
    alias_species_genes_from_reference_to_match_de_novo_paths_new,


rule calculate_species_depth_from_core_genes:
    output:
        species_depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene{specgene_params}.species_depth.tsv",
    input:
        script="scripts/calculate_species_depth_from_core_genes.py",
        species_gene="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene{specgene_params}.species_gene.list",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
    params:
        trim_frac=0.15,  # FIXME: Add this as a parameter to spgc-string.
    shell:
        """
        {input.script} {input.species_gene} {input.gene_depth} {params.trim_frac} {output.species_depth}
        """


rule combine_species_depth_from_select_species:
    output:
        "data/group/{group}/r.{proc}.gene{pangenome_params}.spgc_specgene{specgene_params}.species_depth.tsv"
    input:
        species_depth=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/r.{w.proc}.gene{w.pangenome_params}.spgc_specgene{w.specgene_params}.species_depth.tsv"
            for species in config["species_group"][w.group]
        ],
    params:
        species_list=lambda w: config["species_group"][w.group],
        species_pattern="data/group/{group}/species/sp-$species/r.{proc}.gene{pangenome_params}.spgc_specgene{specgene_params}.species_depth.tsv",
    shell:
        """
        for species in {params.species_list}
        do
            echo -n . >&2
            file={params.species_pattern}
            awk -v OFS='\t' -v species=$species '{{print $1,species,$2}}' $file
        done > {output}
        echo "" >&2
        """



# NOTE: Hub-rule; comment out to reduce DAG building time.
# NOTE: Because I use a species-specific species_depth.tsv,
# in many of these files, samples with 0-depth (or maybe samples with
# 0-depth for ALL species genes), are not even listed.
# Therefore, the final list of "species-free samples" does not include
# these at all.
# TODO: Consider if this is a problem.
rule identify_species_free_samples:
    output:
        "{stemA}/species/sp-{species}/{stemB}.spgc_specgene{specgene_params}.species_free_samples.list",
    input:
        script="scripts/identify_species_free_samples.py",
        species_depth="{stemA}/species/sp-{species}/{stemB}.spgc_specgene{specgene_params}.species_depth.tsv",
    params:
        thresh=0.0001,
    shell:
        """
        {input.script} \
                {params.thresh} \
                {input.species_depth} \
                {output}
        """


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
        frac_thresh=0.95,
    shell:
        """
        awk -v OFS='\t' -v thresh={params.frac_thresh} 'NR == 1 || $3 > thresh {{print $2,$1}}' {input} > {output}
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


# rule partition_strain_samples_method1:
#     output:
#         nospecies="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.species_free_samples.list",
#         strain="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_samples.tsv",
#     input:
#         script="scripts/partition_strain_samples.py",
#         species_depth="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.species_depth.tsv",
#         strain_frac="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.comm.tsv",
#     params:
#         frac_thresh=0.95,
#         absent_thresh=0.0001,
#         present_thresh=0.5,
#     shell:
#         """
#         {input.script} \
#                 {input.species_depth} \
#                 {input.strain_frac} \
#                 {params.frac_thresh} \
#                 {params.absent_thresh} \
#                 {params.present_thresh} \
#                 {output.nospecies} \
#                 {output.strain}
#         """
#
#
# use rule partition_strain_samples_method1 as partition_strain_samples_method2 with:
#     output:
#         nospecies="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.spgc.species_free_samples.list",
#         strain="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.spgc.strain_samples.tsv",
#     input:
#         script="scripts/partition_strain_samples.py",
#         species_depth="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.species_depth.tsv",
#         strain_frac="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.comm.tsv",


rule calculate_strain_specific_correlation_and_depth_ratio_of_genes:
    output:
        corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_correlation.tsv",
        depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_depth_ratio.tsv",
    input:
        script="scripts/calculate_strain_partitioned_gene_stats.py",
        nospecies_samples="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_free_samples.list",
        strain_partitions="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
    params:
        trnsfm=lambda w: float(w.trnsfm) / 10,
    shell:
        """
        {input.script} \
                {input.nospecies_samples} \
                {input.strain_partitions} \
                {input.species_depth} \
                {input.gene_depth} \
                {params.trnsfm} \
                {output.corr} \
                {output.depth}
        """


# rule calculate_correlation_and_depth_quantiles_relative_to_species_genes:
#     output:
#         corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_corr_quantile.tsv",
#         depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_depth_quantile.tsv",
#     input:
#         script="scripts/calculate_strain_gene_scores.py",
#         # species_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.species_gene2-n500.list",
#         species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.species_gene-trim25-prev95.list",
#         strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_correlation.tsv",
#         strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
#     shell:
#         """
#         {input.script} \
#                 {input.species_gene} \
#                 {input.strain_corr} \
#                 {input.strain_depth} \
#                 {output.corr} \
#                 {output.depth}
#         """


# NOTE: Because the species_gene file is not specific to this group, the stem is different. I therefore cannot simply
# parameterize the filename by specgene parameters.
rule pick_strain_gene_thresholds_by_quantiles_clipped:
    output:
        "{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-corrq{corr}-depthq{depth}.strain_gene_threshold.tsv",
    input:
        script="scripts/pick_strain_gene_thresholds.py",
        species_gene="{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_gene.list",  # An alias of the pangenome and sfacts agnostic version.
        strain_corr="{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_correlation.tsv",
        strain_depth="{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_depth_ratio.tsv",
    wildcard_constraints:
        corr=integer_wc,
        depth=integer_wc,
    params:
        strain_corr_quantile=lambda w: float(w.corr) / 1000,
        strain_depth_quantile=lambda w: float(w.depth) / 1000,
        min_corr=0.4,
        max_corr=0.8,
        min_depth=0.1,
        max_depth=0.5,
    shell:
        """
        {input.script} \
                {input.species_gene} \
                {input.strain_corr} \
                {params.strain_corr_quantile} \
                {params.min_corr} \
                {params.max_corr} \
                {input.strain_depth} \
                {params.strain_depth_quantile} \
                {params.min_depth} \
                {params.max_depth} \
                {output}
        """


use rule pick_strain_gene_thresholds_by_quantiles_clipped as pick_strain_gene_thresholds_by_quantiles_clipped_fixed_depth with:
    output:
        "{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-corrq{corr}-depth{depth}.strain_gene_threshold.tsv",
    params:
        strain_corr_quantile=lambda w: float(w.corr) / 1000,
        strain_depth_quantile=0.5,  # NOTE: Dummy value. Min=max; therefore depth is fixed. Resulting depth_threshold_high is actually the median.
        min_depth=lambda w: float(w.depth) / 1000,
        max_depth=lambda w: float(w.depth) / 1000,
        min_corr=0.4,  # NOTE: There's a parallel parameter also set for the parent rule.
        max_corr=0.8,


use rule pick_strain_gene_thresholds_by_quantiles_clipped as pick_strain_gene_thresholds_by_fixed_depth_and_corr with:
    output:
        "{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-corr{corr}-depth{depth}.strain_gene_threshold.tsv",
    params:
        strain_corr_quantile=0.5,  # Dummy value. Min=max; therefore corr is fixed.
        strain_depth_quantile=0.5,  # NOTE: Dummy value. Min=max; therefore depth is fixed. Resulting depth_threshold_high is actually the median.
        min_depth=lambda w: float(w.depth) / 1000,
        max_depth=lambda w: float(w.depth) / 1000,
        min_corr=lambda w: float(w.corr) / 1000,
        max_corr=lambda w: float(w.corr) / 1000,


rule pick_strain_gene_thresholds_by_grid_search:
    output:
        "{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-alpha{alpha}.strain_gene_threshold.tsv",
    input:
        script="scripts/pick_strain_gene_thresholds_by_grid_search.py",
        species_gene="{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_gene.list",
        strain_corr="{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_correlation.tsv",
        strain_depth="{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_depth_ratio.tsv",
    params:
        resolution_corr=50,
        resolution_depth=50,
        alpha=lambda w: float(w.alpha) / 100,
    shell:
        """
        {input.script} \
                {input.species_gene} \
                {input.strain_corr} \
                {input.strain_depth} \
                {params.resolution_corr} \
                {params.resolution_depth} \
                {params.alpha} \
                {output}
        """


rule select_strain_gene_hits:
    output:
        "{stem}.spgc_{spgc_stem}_thresh-{thresh_params}.uhgg-strain_gene.tsv",
    input:
        script="scripts/select_strain_genes.py",
        thresholds="{stem}.spgc_{spgc_stem}_thresh-{thresh_params}.strain_gene_threshold.tsv",
        strain_corr="{stem}.spgc_{spgc_stem}.strain_correlation.tsv",
        strain_depth="{stem}.spgc_{spgc_stem}.strain_depth_ratio.tsv",
    shell:
        """
        {input.script} \
                {input.strain_corr} \
                {input.strain_depth} \
                {input.thresholds} \
                {output}
        """


rule aggregate_uhgg_strain_gene_by_annotation:
    output:
        "{stemA}/species/sp-{species}/{stemB}.{agg}-strain_gene.tsv",
        # "data/group/{group}/species/sp-{species}/{stem}.{agg}-strain_gene.tsv",
    wildcard_constraints:
        agg="eggnog|top_eggnog|cog|ko",
    input:
        script="scripts/aggregate_uhgg_strain_gene_by_annotation.py",
        uhgg="{stemA}/species/sp-{species}/{stemB}.uhgg-strain_gene.tsv",
        # uhgg="data/group/{group}/species/sp-{species}/{stem}.uhgg-strain_gene.tsv",
        agg="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_{agg}.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.uhgg} {input.agg} {wildcards.agg} {output}"




rule convert_midasdb_species_gene_list_to_reference_genome_table_new:
    output:
        "data/species/sp-{species}/gene{centroid}_new.reference_copy_number.nc",
    input:
        script="scripts/convert_gene_info_to_genome_table.py",
        genes="ref/midasdb_uhgg_new/pangenomes/{species}/gene_info.txt",
    shell:
        "{input.script} {input.genes} centroid_{wildcards.centroid} {output}"
