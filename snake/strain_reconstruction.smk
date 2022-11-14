rule combine_midasdb_all_gene_annotations:
    output:
        "ref/midasdb_uhgg.sp-{species}.gene_annotations.tsv",
    input:
        flag="data/species/sp-{species}/download_uhgg_gene_annotations_all_tsv.flag",
    shell:
        """
        find ref/midasdb_uhgg/gene_annotations/{wildcards.species}/ -name '*.tsv.lz4' \
                | xargs lz4cat \
                | awk '$1 != "locus_tag" && $2 != "gene"' \
            > {output}
        """


rule filter_midasdb_all_gene_annotations_by_centroid:
    output:
        "ref/midasdb_uhgg.sp-{species}.gene{centroid}_annotations.tsv",
    input:
        annot="ref/midasdb_uhgg.sp-{species}.gene_annotations.tsv",
        centroids_list="ref/midasdb_uhgg/pangenomes/{species}/gene_info.txt.lz4",
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


rule convert_genes_tally_to_cluster_depth:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.midas_gene.depth.nc",
    input:
        script="scripts/convert_genes_tally_to_depth.py",
        midasdir="data/group/{group}/species/sp-{species}/{stem}.midas_merge/genes",
        meta="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    params:
        inpath="data/group/{group}/species/sp-{species}/{stem}.midas_merge/genes/{species}/{species}.genes_reads.tsv.lz4",
        assumed_read_length=125,
    shell:
        """
        {input.script} \
                <(lz4cat {params.inpath} | tqdm --bytes) \
                {input.meta} \
                {params.assumed_read_length} \
                {output}
        """


rule aggregate_gene_depth_by_centroid:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.depth.nc",
    input:
        script="scripts/aggregate_gene_depth_by_centroid.py",
        depth="data/group/{group}/species/sp-{species}/{stem}.midas_gene.depth.nc",
        meta="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    params:
        aggregate_genes_by=lambda w: {
            "99": "centroid_99",
            "95": "centroid_95",
            "90": "centroid_90",
            "85": "centroid_85",
            "80": "centroid_80",
            "75": "centroid_75",
        }[w.centroid],
    shell:
        """
        {input.script} \
                {input.depth} \
                {input.meta} \
                {params.aggregate_genes_by} \
                {output}
        """


rule calculate_species_correlation_of_genes:
    output:
        corr="data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.species_correlation.tsv",
        sample_depth="data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.species_depth.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.species_depth_ratio.tsv",
        threshold="data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.species_corr_threshold.tsv",
    input:
        script="scripts/calculate_species_correlation_of_genes.py",
        species_depth="data/group/{group}/{stem}.gtpro.species_depth.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stem}.midas_gene{centroid}.depth.nc",
    params:
        transformation_root=0.33,
        n_marker_genes=1000,
        trim_frac=0.05,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {wildcards.species} \
                {input.gene_depth} \
                {params.transformation_root} \
                {params.n_marker_genes} \
                {params.trim_frac} \
                {output.corr} \
                {output.sample_depth} \
                {output.gene_depth} \
                {output.threshold}
        """


rule calculate_strain_specific_correlation_of_genes:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_correlation.tsv",
    input:
        script="scripts/calculate_strain_specific_correlation_of_genes.py",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.species_depth.tsv",
        strain_frac="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
    params:
        strain_frac_thresh=0.95,
        species_depth_thresh_abs=0.0001,
        species_depth_thresh_pres=0.5,
        transformation_root=0.33,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {params.species_depth_thresh_abs} \
                {params.species_depth_thresh_pres} \
                {input.strain_frac} \
                {params.strain_frac_thresh} \
                {input.gene_depth} \
                {params.transformation_root} \
                {output}
        """


rule calculate_cross_species_strain_specific_correlation_of_genes:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_by_species_correlation.nc",
    input:
        script="scripts/calculate_cross_species_strain_specific_correlation_of_genes.py",
        species_depth="data/group/{group}/{stemA}.gtpro.species_depth.tsv",
        strain_frac="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
    params:
        strain_frac_thresh=0.95,
        species_depth_thresh_abs=0.0001,
        species_depth_thresh_pres=0.5,
        transformation_root=0.33,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {wildcards.species} \
                {params.species_depth_thresh_abs} \
                {params.species_depth_thresh_pres} \
                {input.strain_frac} \
                {params.strain_frac_thresh} \
                {input.gene_depth} \
                {params.transformation_root} \
                {output}
        """


rule calculate_strain_specific_gene_depth_ratio:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_depth_ratio.tsv",
    input:
        script="scripts/calculate_strain_specific_gene_depth_ratio.py",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.species_depth.tsv",
        strain_frac="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
    params:
        strain_frac_thresh=0.95,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {input.strain_frac} \
                {params.strain_frac_thresh} \
                {input.gene_depth} \
                {output}
        """


rule pick_strain_gene_thresholds:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_gene_threshold.tsv",
    input:
        script="scripts/pick_strain_gene_thresholds.py",
        species_corr="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.species_correlation.tsv",
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_depth_ratio.tsv",
    params:
        strain_corr_quantile_strict=0.1,
        strain_corr_quantile_moderate=0.05,
        strain_corr_quantile_lenient=0.01,
        strain_depth_quantile=0.05,
        species_corr_threshold=0.98,
    shell:
        """
        {input.script} \
                {input.species_corr} \
                {params.species_corr_threshold} \
                {input.strain_corr} \
                {params.strain_corr_quantile_strict} \
                {params.strain_corr_quantile_moderate} \
                {params.strain_corr_quantile_lenient} \
                {input.strain_depth} \
                {params.strain_depth_quantile} \
                {output}
        """


rule convert_midasdb_species_gene_list_to_reference_genome_table:
    output:
        "data/species/sp-{species}/midas_gene{centroid}.reference_copy_number.nc",
    input:
        script="scripts/convert_gene_info_to_genome_table.py",
        genes="ref/midasdb_uhgg/pangenomes/{species}/gene_info.txt.lz4",
    shell:
        "{input.script} {input.genes} centroid_{wildcards.centroid} {output}"


rule collect_files_for_strain_assessment:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_files.flag",
    input:
        sfacts="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        strain_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_correlation.tsv",
        strain_depth_ratio="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_depth_ratio.tsv",
        strain_fraction="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        gene_centroids="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.species_depth.tsv",
        gtpro_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.species_depth.tsv",
        species_correlation="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.species_correlation.tsv",
        strain_thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_gene_threshold.tsv",
        gene_annotations="ref/midasdb_uhgg.sp-{species}.gene{centroid}_annotations.tsv",
        midas_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
        reference_copy_number="data/species/sp-{species}/midas_gene{centroid}.reference_copy_number.nc",
    shell:
        "echo {input} | tee {output}"
