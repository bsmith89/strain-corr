rule combine_midasdb_all_gene_annotations:
    output:
        "ref/midasdb_uhgg_gene_annotations/sp-{species}.gene_annotations.tsv",
    input:
        "ref/midasdb_uhgg_gene_annotations/{species}",
    shell:
        """
        find {input} -name '*.tsv.lz4' \
                | xargs lz4cat \
                | awk '$1 != "locus_tag" && $2 != "gene"' \
            > {output}
        """


# TODO: The aggregated gene_annotations.tsv files could/should go into each
# individual species directory.
rule filter_midasdb_all_gene_annotations_by_centroid:
    output:
        "ref/midasdb_uhgg_gene_annotations/sp-{species}.gene{centroid}_annotations.tsv",
    input:
        annot="ref/midasdb_uhgg_gene_annotations/sp-{species}.gene_annotations.tsv",
        centroids_list="ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
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


rule select_species_core_genes_from_reference:
    output:
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroid}.species_gene-trim{trim_quantile}-prev{prevalence}.list",
    input:
        script="scripts/select_high_prevalence_species_genes.py",
        copy_number="ref/midasdb_uhgg_pangenomes/{species}/gene{centroid}.reference_copy_number.nc",
    params:
        trim_quantile=lambda w: float(w.trim_quantile) / 100,
        prevalence=lambda w: float(w.prevalence) / 100,
    shell:
        "{input.script} {input.copy_number} {params.trim_quantile} {params.prevalence} {output}"


# TODO: Consider doing this *after* aggregating by strain. This will prevent a
# very very abundant/prevalent strain from dominating and degrading the core
# gene selection.
rule select_species_core_genes_de_novo:
    output:
        species_gene="{stemA}/species/sp-{species}/{stemB}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_gene-n700.list",
        species_corr="{stemA}/species/sp-{species}/{stemB}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_correlation.tsv",
    input:
        script="scripts/select_highly_correlated_species_genes.py",
        species_depth="{stemA}/{stemB}.gtpro.species_depth.tsv",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    params:
        n_marker_genes=700,
    shell:
        """
        {input.script} {input.species_depth} {wildcards.species} {input.gene_depth} {params.n_marker_genes} {output.species_gene} {output.species_corr}
        """


rule calculate_species_depth_from_core_genes:
    output:
        species_depth="{stemA}/species/sp-{species}/{stemB}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",
    input:
        script="scripts/calculate_species_depth_from_core_genes.py",
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.species_gene-trim25-prev95.list",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    params:
        trim_frac=0.05,
    shell:
        """
        {input.script} {input.species_gene} {input.gene_depth} {params.trim_frac} {output.species_depth}
        """


rule partition_strain_samples:
    output:
        nospecies="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.gene{centroidA}-{params}-agg{centroidB}.spgc.species_free_samples.list",
        strain="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_samples.tsv",
    input:
        script="scripts/partition_strain_samples.py",
        species_depth="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",  # NOTE: Consider using GT-Pro depths here instead of species_depth, to avoid circular dependency in picking species genes. Previously: "{stemA}/{stemB}.gtpro.species_depth.tsv",
        strain_frac="data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.comm.tsv",
    params:
        frac_thresh=0.95,
        absent_thresh=0.0001,
        present_thresh=0.5,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {input.strain_frac} \
                {params.frac_thresh} \
                {params.absent_thresh} \
                {params.present_thresh} \
                {output.nospecies} \
                {output.strain}
        """


rule calculate_strain_specific_correlation_and_depth_ratio_of_genes:
    output:
        corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
        depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
    input:
        script="scripts/calculate_strain_partitioned_gene_stats.py",
        nospecies_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.species_free_samples.list",
        strain_partitions="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_samples.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    shell:
        """
        {input.script} \
                {input.nospecies_samples} \
                {input.strain_partitions} \
                {input.species_depth} \
                {input.gene_depth} \
                {output.corr} \
                {output.depth}
        """


rule calculate_correlation_and_depth_quantiles_relative_to_species_genes:
    output:
        corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_corr_quantile.tsv",
        depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_quantile.tsv",
    input:
        script="scripts/calculate_strain_gene_scores.py",
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.species_gene-trim25-prev95.list",  # FIXME: Are these the right way to define a species core gene?
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
    shell:
        """
        {input.script} \
                {input.species_gene} \
                {input.strain_corr} \
                {input.strain_depth} \
                {output.corr} \
                {output.depth}
        """


rule pick_strain_gene_thresholds:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc-corr{corr_quant}-depth{depth_quant}.strain_gene_threshold.tsv",
    input:
        script="scripts/pick_strain_gene_thresholds.py",
        species_gene="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.species_gene-trim25-prev95.list",  # FIXME: Are these the right way to define a species core gene? # species_gene="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_gene.list",
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
    params:
        strain_corr_quantile=lambda w: float(w.corr_quant) / 1000,
        strain_depth_quantile=lambda w: float(w.depth_quant) / 1000,
        min_corr=0.2,
        max_corr=0.9,
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


rule convert_midasdb_species_gene_list_to_reference_genome_table:
    output:
        "ref/midasdb_uhgg_pangenomes/{species}/gene{centroid}.reference_copy_number.nc",
    input:
        script="scripts/convert_gene_info_to_genome_table.py",
        genes="ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
    shell:
        "{input.script} {input.genes} centroid_{wildcards.centroid} {output}"


rule assess_infered_strain_accuracy:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.{strain}.gene_content_reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        gene_matching="data/species/sp-{species}/genome/{strain}.midas_uhgg_pangenome-blastn.gene_matching-c{centroidB}-t95.tsv",
        thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc-corr10-depth10.strain_gene_threshold.tsv",
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
    shell:
        """
        {input.script} \
                {input.gene_matching} \
                {input.strain_corr} \
                {input.strain_depth} \
                {input.thresholds} \
                {output}
        """


rule assess_strain_accuracy_for_species:
    output:
        "data/group/xjin_hmp2/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.gene_content_reconstruction_accuracy.ALL_STRAINS.flag",
    input:
        reference_strain_accuracy=lambda w: [
            f"data/group/xjin_hmp2/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.params}-agg{w.centroidB}.spgc.{strain}.gene_content_reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],


rule collect_files_for_strain_assessment:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.refit-{stemC}.gene{centroidA}-{params}-agg{centroidB}.spgc-corr{corr_quant}-depth{depth_quant}.strain_files.flag",
    input:
        sfacts="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        refit="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.refit-{stemC}.world.nc",
        strain_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
        strain_depth_ratio="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
        strain_fraction="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",
        gtpro_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.species_depth.tsv",
        species_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_correlation.tsv",
        species_gene_de_novo="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_gene-n700.list",
        species_gene_reference="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.species_gene-trim25-prev95.list",
        strain_thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc-corr{corr_quant}-depth{depth_quant}.strain_gene_threshold.tsv",
        strain_corr_quantile="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_corr_quantile.tsv",
        strain_depth_quantile="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_quantile.tsv",
        gene_annotations="ref/midasdb_uhgg_gene_annotations/sp-{species}.gene{centroidB}_annotations.tsv",
        depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
        reference_copy_number="ref/midasdb_uhgg_pangenomes/{species}/gene{centroidB}.reference_copy_number.nc",
        gtpro_reference_genotype="data/species/sp-{species}/gtpro_ref.mgtp.nc",
        strain_blastn_midas=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.midas_uhgg_pangenome-blastn.tsv"
            for strain in species_genomes(w.species)
        ],
        strain_blastn_self=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.{strain}-blastn.tsv"
            for strain in species_genomes(w.species)
        ],
        strain_blastn_ratio=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.midas_uhgg_pangenome-blastn.bitscore_ratio-c75.tsv"
            for strain in species_genomes(w.species)
        ],
        strain_gene_lengths=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.prodigal-single.cds.nlength.tsv"
            for strain in species_genomes(w.species)
        ],
        strain_genotype=lambda w: [
            f"data/species/sp-{w.species}/strain_genomes.gtpro.mgtp.nc"
        ],
        reference_strain_accuracy=lambda w: [
            f"data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.params}-agg{w.centroidB}.spgc.{strain}.gene_content_reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
        reference_strain_mapping="data/group/{group}/species/sp-{species}/ALL_STRAINS.tiles-l100-o99.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    params:
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        "echo {input} {params.cluster_info} | tee {output}"
