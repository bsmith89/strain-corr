# TODO: Write script.
# This will shuffle the order of lines in the samples map (reproducibly based on seed)
# and then take the top n_samples for each strain and write out a strain_samples
# file, but with fewer lines.
rule subsample_strain_samples_for_benchmarking:
    output:
        "{stem}.spgc-ss-seed{seed}-n{n_samples}.strain_samples.tsv",
    input:
        script="scripts/subsample_strain_samples_for_benchmarking.py",
        samples="{stem}.spgc-allsamples.strain_samples.tsv",
    params:
        seed=lambda w: int(w.seed),
        n_samples=lambda w: int(w.n_samples),
    shell:
        """
        {input.script} {input.seed} {input.n_samples} {input.samples} {output}
        """


# # TODO: Parameterize this rule (or a child of this rule) by the partition of
# # the strain_samples map.
# # NOTE: I could also parameterize the main rule, instead of splitting this one
# # out, but it will require handling filenames somehow.
# use rule calculate_strain_specific_correlation_and_depth_ratio_of_genes as calculate_strain_specific_correlation_and_depth_ratio_of_genes_with_subsampling with:
#     output:
#         corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_correlation.tsv",
#         depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{params}-agg{centroidB}.spgc.strain_depth_ratio.tsv",
#     input:
#         script="scripts/calculate_strain_partitioned_gene_stats.py",
#         nospecies_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_free_samples.list",
#         strain_partitions="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.strain_samples.tsv",
#         species_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.gene{centroidA}-{params}-agg{centroidB}.spgc.species_depth.tsv",
#         gene_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",


rule assess_infered_strain_accuracy:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_thresh-{thresh_params}.{strain}.gene_content_reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        # script="scripts/assess_gene_content_reconstruction_accuracy2.py",
        # gene_matching="data/species/sp-{species}/genome/{strain}.midas_uhgg_pangenome-blastn.gene_matching-c{centroidB}-t95.tsv",
        gene_matching="data/group/{group}/species/sp-{species}/genome/{strain}.tiles-l100-o99.gene{centroidA}-{bowtie_params}-agg{centroidB}.gene_matching-t30.tsv",
        thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_thresh-{thresh_params}.strain_gene_threshold.tsv",
        strain_corr="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}.strain_correlation.tsv",
        strain_depth="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}.strain_depth_ratio.tsv",
    shell:
        """
        {input.script} \
                {input.gene_matching} \
                {input.strain_corr} \
                {input.strain_depth} \
                {input.thresholds} \
                {output}
        """


# rule assess_strain_accuracy_for_species:
#     output:
#         "data/group/xjin_hmp2/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.gene_content_reconstruction_accuracy.ALL_STRAINS.flag",
#     input:
#         reference_strain_accuracy=lambda w: [
#             f"data/group/xjin_hmp2/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.params}-agg{w.centroidB}.spgc.{strain}.gene_content_reconstruction_accuracy.tsv"
#             for strain in species_genomes(w.species)
#         ],


rule collect_files_for_strain_assessment:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_thresh-corr{corr_thresh}-depth{depth_thresh}.strain_files.flag",
    input:
        sfacts="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        strain_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}.strain_correlation.tsv",
        strain_depth_ratio="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}.strain_depth_ratio.tsv",
        strain_fraction="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        # strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_samples.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
        gtpro_depth="data/group/{group}/{stemA}.gtpro.species_depth.tsv",
        species_correlation="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo2-t30-n500.species_correlation.tsv",
        species_gene_de_novo="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo-n500.species_gene.list",
        species_gene_de_novo2="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo2-t30-n500.species_gene.list",
        species_gene_reference="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.spgc_specgene-ref-t25-p95.species_gene.list",
        strain_thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_thresh-corr{corr_thresh}-depth{depth_thresh}.strain_gene_threshold.tsv",
        # strain_corr_quantile="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_corr_quantile.tsv",
        # strain_depth_quantile="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_depth_quantile.tsv",
        gene_annotations="ref/midasdb_uhgg_gene_annotations/sp-{species}.gene{centroidB}_annotations.tsv",
        depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
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
            f"data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.bowtie_params}-agg{w.centroidB}.spgc_specgene-{w.specgene_params}_ss-{w.ss_params}_thresh-corr{w.corr_thresh}-depth{w.depth_thresh}.{strain}.gene_content_reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
        reference_strain_accuracy_depth_only=lambda w: [
            f"data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.bowtie_params}-agg{w.centroidB}.spgc_specgene-{w.specgene_params}_ss-{w.ss_params}_thresh-corr0-depth{w.depth_thresh}.{strain}.gene_content_reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
        reference_strain_mapping="data/group/{group}/species/sp-{species}/ALL_STRAINS.tiles-l100-o99.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
    params:
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        "echo {input} {params.cluster_info} | tr ' ' '\n' > {output}"
