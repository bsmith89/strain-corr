rule select_xjin_samples:
    output:
        "{stem}.spgc_ss-xjin-all.strain_samples.tsv",
    input:
        "{stem}.spgc_ss-all.strain_samples.tsv",
    params:
        pattern="^xjin_",
    shell:
        """
        awk 'NR==1 || $2~/{params.pattern}/ {{print $0}}' {input} > {output}
        """


# This will shuffle the order of lines in the samples map (reproducibly based on seed)
# and then take the top n_samples for each strain and write out a strain_samples
# file, but with fewer lines.
rule subsample_xjin_samples_for_benchmarking:
    output:
        "{stem}.spgc_ss-xjin-seed{seed}-n{n_samples}.strain_samples.tsv",
    input:
        script="scripts/subsample_strain_samples_for_benchmarking.py",
        samples="{stem}.spgc_ss-xjin-all.strain_samples.tsv",
    params:
        seed=lambda w: int(w.seed),
        n_samples=lambda w: int(w.n_samples),
    shell:
        """
        {input.script} {params.seed} {params.n_samples} {input.samples} {output}
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


rule select_highest_depth_xjin_samples_for_benchmarking:
    output:
        "{stemA}.gtpro.{stemB}.spgc_ss-xjin-deepest-n{n_samples}.strain_samples.tsv",
    input:
        script="scripts/select_highest_depth_strain_samples_for_benchmarking.py",
        depth="{stemA}.gtpro.species_depth.tsv",
        samples="{stemA}.gtpro.{stemB}.spgc_ss-xjin-all.strain_samples.tsv",
    params:
        n_samples=lambda w: int(w.n_samples),
    shell:
        """
        {input.script} {params.n_samples} {input.depth} {input.samples} {output}
        """


rule alias_spgc_gene_hits_as_uhgg_strain_gene:
    output:
        "{stemA}.spgc{stemB}.uhgg-strain_gene.tsv",
    input:
        "{stemA}.spgc{stemB}.strain_gene.tsv",
    shell:
        alias_recipe


ruleorder: alias_spgc_gene_hits_as_uhgg_strain_gene > aggregate_uhgg_strain_gene_by_annotation


localrules:
    alias_spgc_gene_hits_as_uhgg_strain_gene,


rule aggregate_uhgg_strain_gene_by_annotation:
    output:
        "data/group/{group}/species/sp-{species}/{stem}.{agg}-strain_gene.tsv",
    wildcard_constraints:
        agg="eggnog|top_eggnog|cog|ko",
    input:
        uhgg="data/group/{group}/species/sp-{species}/{stem}.uhgg-strain_gene.tsv",
        agg="data/species/sp-{species}/pangenome.centroids.emapper.gene_x_{agg}.tsv",
    group:
        "assess_gene_inference_benchmark"
    run:
        gene_table = pd.read_table(input.uhgg, index_col="gene_id")
        agg = pd.read_table(input.agg, index_col="gene_id").squeeze()
        result = (
            gene_table.astype(bool).join(agg).groupby(wildcards.agg).any().astype(int)
        )
        result.to_csv(output[0], sep="\t")


# NOTE: (2023-06-13) This renaming brings spgc inferences into alignment with spanda and panphlan so
# that the accuracy assessment rules below can be generic between all three.
rule alias_spgc_xjin_hmp2_inferences_as_xjin_benchmark:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_{stemB}.uhgg-strain_gene.tsv",
    input:
        "data/group/xjin_hmp2/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_{stemB}.uhgg-strain_gene.tsv",
    shell:
        alias_recipe


# NOTE: (2023-06-13) This rule order should mean that
# alias_spgc_gene_hits_as_uhgg_strain_gene doesn't look for an XJIN_BENCHMARK
# inference first.
ruleorder: alias_spgc_xjin_hmp2_inferences_as_xjin_benchmark > alias_spgc_gene_hits_as_uhgg_strain_gene


# NOTE: (2023-06-20) UHGG accuracy gets its own rule, separate from the
# other units, because it's assigned based on tiling depth with a particular
# centroidA, centroidB, mapping strategy, etc.
# Also note that UHGG as an assessment unit isn't great. I could try switching
# to "best hit" UHGG assignment, but I think the accuracy would look
# misleadingly low...
# TODO: Probably worth considering using the BLAST-based strain-ORF annotation
# approach, though.
rule assess_infered_strain_accuracy_uhgg:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.uhgg-reconstruction_accuracy.tsv",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.uhgg-strain_gene.tsv",
        truth="data/group/xjin_hmp2/species/sp-{species}/genome/{strain}.tiles-l100-o99.gene{pangenome_params}.gene_matching-t30.uhgg-strain_gene.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        """
        {input.script} \
                {input.truth} \
                {input.infer} \
                {output}
        """


ruleorder: assess_infered_strain_accuracy_uhgg > assess_infered_strain_accuracy_other_unit


use rule assess_infered_strain_accuracy_uhgg as assess_infered_strain_accuracy_other_unit with:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{strain}.{unit}-reconstruction_accuracy.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko",
    input:
        script="scripts/assess_gene_content_reconstruction_accuracy.py",
        infer="data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.{stemB}.{unit}-strain_gene.tsv",
        truth="data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.{unit}-strain_gene.tsv",


# # NOTE: This hard-codes the "xjin_" prefix and finds the strain with the most
# # of these samples.
# # TODO: (2023-06-13) Decide if I should be selecting the panphlan/spanda strain
# # based on something like this.
# rule match_xjin_strains:
#     output:
#         samples="data/group/xjin_hmp2/{stem}.spgc_ss-{ss_params}.xjin_sample_count.tsv",
#     input:
#         samples="data/group/xjin_hmp2/{stem}.spgc_ss-{ss_params}.strain_samples.tsv",
#     params:
#         pattern="^xjin_"
#     shell:
#         """
#         awk '$2~/{params.pattern}/ {{print $0}}' {input} \
#                 | cut -f1 \
#                 | sort \
#                 | uniq -c \
#                 | awk -v OFS="\t" '{{print $2,$1}}' \
#             > {output}
#         """


rule compile_reference_genome_accuracy_info_for_spgc:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{pangenome_params}.spgc_specgene-{specgene_params}_ss-xjin-{ss_params}_t-{trnsfm}_thresh-{thresh_params}.{unit}-xjin_strain_summary.tsv",
    input:
        script="scripts/compile_reference_genome_accuracy_info.py",
        strain_samples="data/group/xjin_hmp2/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-xjin-{ss_params}.strain_samples.tsv",
        species_free_samples="data/group/xjin_hmp2/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_specgene-{specgene_params}.species_free_samples.list",
        species_depth="data/group/xjin_hmp2/species/sp-{species}/{stemA}.gene{pangenome_params}.spgc_specgene-{specgene_params}.species_depth.tsv",
        strain_thresh="data/group/xjin_hmp2/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{pangenome_params}.spgc_specgene-{specgene_params}_ss-xjin-{ss_params}_t-{trnsfm}_thresh-{thresh_params}.strain_gene_threshold.tsv",
        reference_genome_accuracy=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.pangenome_params}.spgc_specgene-{w.specgene_params}_ss-xjin-{w.ss_params}_t-{w.trnsfm}_thresh-{w.thresh_params}.{strain}.{w.unit}-reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
        strain_meta="data/group/xjin_hmp2/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{pangenome_params}.spgc_specgene-{specgene_params}_ss-xjin-{ss_params}_t-{trnsfm}_thresh-{thresh_params}.strain_meta.tsv",
    params:
        args=lambda w: [
            f"{strain}=data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.pangenome_params}.spgc_specgene-{w.specgene_params}_ss-xjin-{w.ss_params}_t-{w.trnsfm}_thresh-{w.thresh_params}.{strain}.{w.unit}-reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.strain_meta} {input.strain_samples} {input.species_free_samples} {input.species_depth} {input.strain_thresh} {params.args} > {output}"


# TODO: Combine all reference genomes into one table.
rule compile_reference_genome_accuracy_info_for_spanda:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.spanda{spanda_params}.{unit}-xjin_strain_summary.tsv",
    input:
        reference_genome_accuracy=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gene{w.pangenome_params}.spanda{w.spanda_params}.{strain}.{w.unit}-reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
    params:
        pattern=lambda w: f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gene{w.pangenome_params}.spanda{w.spanda_params}.{{strain}}.{w.unit}-reconstruction_accuracy.tsv",
        strain_list=lambda w: species_genomes(w.species),
    group:
        "assess_gene_inference_benchmark"
    run:
        result = []
        for genome_id in params.strain_list:
            path = params.pattern.format(strain=genome_id)
            data = pd.read_table(path).assign(genome_id=genome_id)
            highest_f1_strain = data.set_index("strain").f1.idxmax()
            data = data.assign(
                dominant_strain=highest_f1_strain
            )  # Don't try and pick the highest-abundance strain.
            result.append(data)
        result = pd.concat(result).assign(
            total_num_reference_genomes=len(params.strain_list)
        )
        result.to_csv(output[0], sep="\t", index=False)


# TODO: Combine all reference genomes into one table.
rule compile_reference_genome_accuracy_info_for_panphlan:
    output:
        "data/group/XJIN_BENCHMARK/species/sp-{species}/{stemA}.gene{pangenome_params}.panphlan.{unit}-xjin_strain_summary.tsv",
    input:
        reference_genome_accuracy=lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gene{w.pangenome_params}.panphlan.{strain}.{w.unit}-reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
    params:
        pattern=lambda w: f"data/group/XJIN_BENCHMARK/species/sp-{w.species}/{w.stemA}.gene{w.pangenome_params}.panphlan.{{strain}}.{w.unit}-reconstruction_accuracy.tsv",
        strain_list=lambda w: species_genomes(w.species),
    group:
        "assess_gene_inference_benchmark"
    run:
        result = []
        for genome_id in params.strain_list:
            path = params.pattern.format(strain=genome_id)
            data = pd.read_table(path).assign(genome_id=genome_id)
            highest_f1_strain = data.set_index("strain").f1.idxmax()
            data = data.assign(
                dominant_strain=highest_f1_strain
            )  # Don't try and pick the highest-abundance strain.
            result.append(data)
        result = pd.concat(result).assign(
            total_num_reference_genomes=len(params.strain_list)
        )
        result.to_csv(output[0], sep="\t", index=False)


rule collect_files_for_strain_assessment:
    output:
        "data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-corr{corr_thresh}-depth{depth_thresh}.strain_files.flag",
    input:
        sfacts="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.world.nc",
        strain_correlation="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_correlation.tsv",
        strain_depth_ratio="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}.strain_depth_ratio.tsv",
        strain_fraction="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.comm.tsv",
        strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.spgc_ss-{ss_params}.strain_samples.tsv",
        # strain_samples="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc.strain_samples.tsv",
        species_depth="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_depth.tsv",
        gtpro_depth="data/group/{group}/{stemA}.gtpro.species_depth.tsv",
        species_correlation="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo2-t30-n500.species_correlation.tsv",
        species_gene="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}.species_gene.list",
        species_gene_de_novo="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo-n500.species_gene.list",
        species_gene_de_novo2="data/group/{group}/species/sp-{species}/{stemA}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-denovo2-t30-n500.species_gene.list",
        species_gene_reference="data/species/sp-{species}/midasuhgg.pangenome.gene{centroidB}.spgc_specgene-ref-t25-p95.species_gene.list",
        strain_thresholds="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-{ss_params}_t-{trnsfm}_thresh-corr{corr_thresh}-depth{depth_thresh}.strain_gene_threshold.tsv",
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
            f"data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.bowtie_params}-agg{w.centroidB}.spgc_specgene-{w.specgene_params}_ss-{w.ss_params}_t-{w.trnsfm}_thresh-corr{w.corr_thresh}-depth{w.depth_thresh}.{strain}.gene_content_reconstruction_accuracy.tsv"
            for strain in species_genomes(w.species)
        ],
        # reference_strain_accuracy_depth_only=lambda w: [
        #     f"data/group/{w.group}/species/sp-{w.species}/{w.stemA}.gtpro.{w.stemB}.gene{w.centroidA}-{w.bowtie_params}-agg{w.centroidB}.spgc_specgene-{w.specgene_params}_ss-{w.ss_params}_t-{w.trnsfm}_thresh-corr0-depth{w.depth_thresh}.{strain}.gene_content_reconstruction_accuracy.tsv"
        #     for strain in species_genomes(w.species)
        # ],
        # xjin_strain_summary="data/group/{group}/species/sp-{species}/{stemA}.gtpro.{stemB}.gene{centroidA}-{bowtie_params}-agg{centroidB}.spgc_specgene-{specgene_params}_ss-xjin-{ss_params}_t-{trnsfm}_thresh-corr{corr_thresh}-depth{depth_thresh}.xjin_strain_summary.tsv"
        reference_strain_mapping="data/group/{group}/species/sp-{species}/ALL_STRAINS.tiles-l100-o99.gene{centroidA}-{bowtie_params}-agg{centroidB}.depth2.nc",
    params:
        cluster_info="ref/midasdb_uhgg/pangenomes/{species}/cluster_info.txt",
    shell:
        "echo {input} {params.cluster_info} | tr ' ' '\n' > {output}"


rule xjin_benchmarking_grid_single_species_single_unit:
    output:
        touch(
            "data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.{unit}-accuracy.xjin_benchmark_grid.flag"
        ),
    input:
        spgc0=lambda w: [
            (
                f"data/group/XJIN_BENCHMARK/{w.stemA}.gtpro.{w.stemB}.gene{w.pangenome_params}"
                f".spgc_specgene-{specgene_params}"
                f"_ss-xjin-deepest-n{max_samples_param}"
                f"_t-{trnsfm_params}"
                f"_thresh-corr{corr_thresh_param}-depth{depth_thresh_param}"
                f".{w.unit}-xjin_strain_summary.tsv"
            )
            for specgene_params, max_samples_param, trnsfm_params, depth_thresh_param, corr_thresh_param, in product(
                ["ref-t25-p95"],
                [1, 10, 20, 999],
                [10, 30],
                [250],
                [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550],
            )
        ],
        # # NOTE: (2023-06-13) This was previously for comparing reference and
        # # denovo species-gene picking.
        # spgc1=lambda w: [
        #     f"data/group/XJIN_BENCHMARK/{w.stemA}.gtpro.{w.stemB}.gene{w.pangenome_params}.spgc_specgene-{specgene_params}_ss-xjin-{ss_params}_t-{trnsfm_params}_thresh-{thresh_params}.{w.unit}-xjin_strain_summary.tsv"
        #     for specgene_params, ss_params, trnsfm_params, thresh_params in product(
        #         ["denovo2-t30-n800"],
        #         ["all"],
        #         [30],
        #         ["corr200-depth250"],
        #     )
        # ],
        spanda="data/group/XJIN_BENCHMARK/{stemA}.gene{pangenome_params}.spanda-s2.{unit}-xjin_strain_summary.tsv",
        panphlan="data/group/XJIN_BENCHMARK/{stemA}.gene{pangenome_params}.panphlan.{unit}-xjin_strain_summary.tsv",


rule xjin_benchmarking_grid_single_species:
    output:
        touch(
            "data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.all_accuracy_xjin_benchmark_grid.flag"
        ),
    input:
        uhgg="data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.uhgg-accuracy.xjin_benchmark_grid.flag",
        ko="data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.ko-accuracy.xjin_benchmark_grid.flag",
        cog="data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.cog-accuracy.xjin_benchmark_grid.flag",
        top_eggnog="data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.top_eggnog-accuracy.xjin_benchmark_grid.flag",
        eggnog="data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.eggnog-accuracy.xjin_benchmark_grid.flag",


rule xjin_benchmarking_grid_all_species:
    output:
        touch(
            "data/group/XJIN_BENCHMARK/{stemA}.gtpro.{stemB}.gene{pangenome_params}.all_accuracy_xjin_benchmark_grid.ALL_XJIN_SPECIES.flag"
        ),
    input:
        lambda w: [
            f"data/group/XJIN_BENCHMARK/species/sp-{species}/{w.stemA}.gtpro.{w.stemB}.gene{w.pangenome_params}.all_accuracy_xjin_benchmark_grid.flag"
            for species in config["genome"].species_id.unique()
            if species != "TODO"
        ],
