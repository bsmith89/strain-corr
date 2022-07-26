rule convert_genes_tally_to_cluster_depth:
    output:
        "data_temp/sp-{species}.{group}.a.{stem}.midas_gene.depth.nc",
    input:
        script="scripts/convert_genes_tally_to_depth.py",
        midasdir="data_temp/sp-{species}.{group}.a.{stem}.midas_merge2/genes",
        meta="ref_temp/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    params:
        inpath="data_temp/sp-{species}.{group}.a.{stem}.midas_merge2/genes/{species}/{species}.genes_reads.tsv.lz4",
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
        "data_temp/sp-{species}.{stem}.midas_gene{centroid}.depth.nc",
    input:
        script="scripts/aggregate_gene_depth_by_centroid.py",
        depth="data_temp/sp-{species}.{stem}.midas_gene.depth.nc",
        meta="ref_temp/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
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
        corr="data_temp/sp-{species}.{stem}.midas_gene{centroid}.species_correlation.tsv",
        sample_depth="data_temp/sp-{species}.{stem}.midas_gene{centroid}.species_depth.tsv",
        gene_depth="data_temp/sp-{species}.{stem}.midas_gene{centroid}.species_depth_ratio.tsv",
    input:
        script="scripts/calculate_species_correlation_of_genes.py",
        species_depth="data/{stem}.gtpro.species_depth.tsv",
        gene_depth="data_temp/sp-{species}.{stem}.midas_gene{centroid}.depth.nc",
    params:
        transformation_root=3,
        corr_thresh=0.95,
        trim_frac=0.25,
    shell:
        """
        {input.script} \
                {input.species_depth} \
                {wildcards.species} \
                {input.gene_depth} \
                {params.transformation_root} \
                {params.corr_thresh} \
                {params.trim_frac} \
                {output.corr} \
                {output.sample_depth} \
                {output.gene_depth}
        """


rule calculate_strain_specific_correlation_of_genes:
    output:
        "data_temp/sp-{species}.{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_correlation.tsv",
    input:
        script="scripts/calculate_strain_specific_correlation_of_genes.py",
        species_depth="data_temp/sp-{species}.{stemA}.midas_gene75.species_depth.tsv",  # NOTE: gene75-based estimates
        strain_frac="data_temp/sp-{species}.{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data_temp/sp-{species}.{stemA}.midas_gene{centroid}.depth.nc",
    params:
        strain_frac_thresh=0.95,
        species_depth_thresh_abs=0.01,
        species_depth_thresh_pres=2.0,
        transformation_root=3,
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


rule calculate_strain_specific_gene_depth_ratio:
    output:
        "data_temp/sp-{species}.{stemA}.gtpro.{stemB}.midas_gene{centroid}.strain_depth_ratio.tsv",
    input:
        script="scripts/calculate_strain_specific_gene_depth_ratio.py",
        species_depth="data_temp/sp-{species}.{stemA}.midas_gene75.species_depth.tsv",  # NOTE: gene75-based estimates
        strain_frac="data_temp/sp-{species}.{stemA}.gtpro.{stemB}.comm.tsv",
        gene_depth="data_temp/sp-{species}.{stemA}.midas_gene{centroid}.depth.nc",
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
