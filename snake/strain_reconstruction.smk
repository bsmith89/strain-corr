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


rule convert_genes_tally_from_multi_to_cluster_depth:
    output:
        "{stemA}/species/sp-{species}/{stemB}.midas_gene.depth.nc",
    input:
        script="scripts/convert_genes_tally_to_depth.py",
        midasdir="{stemA}/species/sp-{species}/{stemB}.midas_merge_from_multi/genes",
        midasdb=ancient("ref/midasdb_uhgg"),
    params:
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
        inpath="{stemA}/species/sp-{species}/{stemB}.midas_merge_from_multi/genes/{species}/{species}.genes_reads.tsv.lz4",
        assumed_read_length=125,
    resources:
        mem_mb=64_000,
        pmem=64_000,
    shell:
        """
        {input.script} \
                <(lz4cat {params.inpath} | tqdm --bytes) \
                {params.cluster_info} \
                {params.assumed_read_length} \
                {output}
        """


rule aggregate_gene_depth_by_centroid:
    output:
        "{stemA}/species/sp-{species}/{stemB}.midas_gene{centroid}.depth.nc",
    input:
        script="scripts/aggregate_gene_depth_by_centroid.py",
        depth="{stemA}/species/sp-{species}/{stemB}.midas_gene.depth.nc",
        midasdb=ancient("ref/midasdb_uhgg"),
    params:
        aggregate_genes_by=lambda w: {
            "99": "centroid_99",
            "95": "centroid_95",
            "90": "centroid_90",
            "85": "centroid_85",
            "80": "centroid_80",
            "75": "centroid_75",
        }[w.centroid],
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        """
        {input.script} \
                {input.depth} \
                {params.cluster_info} \
                {params.aggregate_genes_by} \
                {output}
        """


rule calculate_species_correlation_of_genes:
    output:
        corr="{stemA}/species/sp-{species}/{stemB}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_correlation.tsv",
        sample_depth="{stemA}/species/sp-{species}/{stemB}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_depth.tsv",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_depth_ratio.tsv",
        threshold="{stemA}/species/sp-{species}/{stemB}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_corr_threshold.tsv",
    input:
        script="scripts/calculate_species_correlation_of_genes.py",
        species_depth="{stemA}/{stemB}.gtpro.species_depth.tsv",
        gene_depth="{stemA}/species/sp-{species}/{stemB}.midas_gene{centroid}.depth.nc",
    params:
        exponent=lambda w: float(w.exponent) / 100,
        n_marker_genes=1000,
        trim_frac=0.25,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {wildcards.species} \
                {input.gene_depth} \
                {params.exponent} \
                {params.n_marker_genes} \
                {params.trim_frac} \
                {output.corr} \
                {output.sample_depth} \
                {output.gene_depth} \
                {output.threshold}
        """


rule calculate_strain_specific_correlation_of_genes:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_correlation.tsv",
    input:
        script="scripts/calculate_strain_specific_correlation_of_genes.py",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_depth.tsv",
        strain_frac="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
    params:
        strain_frac_thresh=0.95,
        species_depth_thresh_abs=0.0001,
        species_depth_thresh_pres=0.5,
        exponent=lambda w: float(w.exponent) / 100,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {params.species_depth_thresh_abs} \
                {params.species_depth_thresh_pres} \
                {input.strain_frac} \
                {params.strain_frac_thresh} \
                {input.gene_depth} \
                {params.exponent} \
                {output}
        """




rule calculate_strain_specific_gene_depth_ratio:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_depth_ratio.tsv",
    input:
        script="scripts/calculate_strain_specific_gene_depth_ratio.py",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_depth.tsv",
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
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_gene_threshold.tsv",
    input:
        script="scripts/pick_strain_gene_thresholds.py",
        species_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_correlation.tsv",
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_depth_ratio.tsv",
    params:
        strain_corr_quantile_strict=0.075,
        strain_corr_quantile_moderate=0.025,
        strain_corr_quantile_lenient=0.01,
        strain_depth_quantile=0.05,
        n_marker_genes=1000,
    shell:
        """
        {input.script} \
                {input.species_corr} \
                {params.n_marker_genes} \
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
        "ref/midasdb_uhgg_pangenomes/{species}/midas_gene{centroid}.reference_copy_number.nc",
    input:
        script="scripts/convert_gene_info_to_genome_table.py",
        genes="ref/midasdb_uhgg_pangenomes/{species}/gene_info.txt.lz4",
    shell:
        "{input.script} {input.genes} centroid_{wildcards.centroid} {output}"


rule collect_files_for_strain_assessment:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_files.flag",
    input:
        sfacts="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        strain_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_correlation.tsv",
        strain_depth_ratio="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_depth_ratio.tsv",
        strain_fraction="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_depth.tsv",
        gtpro_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.species_depth.tsv",
        species_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.midas_gene{centroid}.spgc-e{exponent}.species_correlation.tsv",
        strain_thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.midas_gene{centroid}.spgc-e{exponent}.strain_gene_threshold.tsv",
        gene_annotations="ref/midasdb_uhgg_gene_annotations/sp-{species}.gene{centroid}_annotations.tsv",
        midas_depth="data/group/{group}/species/sp-{species}/{stemA}.midas_gene{centroid}.depth.nc",
        reference_copy_number="ref/midasdb_uhgg_pangenomes/{species}/midas_gene{centroid}.reference_copy_number.nc",
        midasdb=ancient("ref/midasdb_uhgg"),
        gtpro_reference_genotype="data/species/sp-{species}/gtpro_ref.mgtp.nc",
    params:
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        "echo {input} {params.cluster_info} | tee {output}"
