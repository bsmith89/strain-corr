# NOTE: Emapper can take a long time to run.
# Testing can be done on 100035 because it has very few genes in its pangenome.
rule eggnog_mapper_translated_orfs:
    output:
        "{stem}.emapper.d/proteins.emapper.annotations",
        "{stem}.emapper.d/proteins.emapper.hits",
        "{stem}.emapper.d/proteins.emapper.seed_orthologs",
    input:
        fasta="{stem}.tran.fa",
        db="ref/eggnog_mapper_db",
    params:
        outdir="{stem}.emapper.d",
        tax_scope="auto",
        sensmode="more-sensitive",
        mapper="diamond",
    conda:
        "conda/eggnog.yaml"
    threads: 24
    resources:
        walltime_hr=240,
        mem_mb=20_000,
        pmem=20_000 // 24,
    shell:
        """
        tmpdir=$(mktemp -d)
        export EGGNOG_DATA_DIR={input.db}
        rm -rf {params.outdir}.temp {params.outdir}
        mkdir -p {params.outdir}.temp
        emapper.py \
                -m {params.mapper} \
                -i {input.fasta} \
                --itype proteins \
                --sensmode {params.sensmode} \
                --go_evidence all \
                --dbmem \
                --tax_scope {params.tax_scope} \
                --temp_dir $tmpdir \
                --override \
                --cpu {threads} \
                --output_dir {params.outdir}.temp \
                --output 'proteins'
        mv {params.outdir}.temp {params.outdir}
        # TODO  # Test on 100035 because it has very few genes
        """


# note: this rule takes the entire eggnog output and assigns raw annotations
# to the features. in a later step, i'll aggregate these annotations
# to the higher-level-centroid by voting or something.
# note: the */eggnog.tsv file is annotatins on c99s, not unique features.
# todo: rename all these scripts to "..._emapper_output_to_gene99_x...".
rule parse_midasdb_emapper_annotations_to_gene99_x_unit_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.emapper.gene99_x_{unit}.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_{unit}.py",
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
    shell:
        "{input.script} {input.emapper} {output}"


# Specialized script for gene_x_cog_category.
# Required because the COG categories assigned by emapper are not the same as those found
# in ref/cog-20.meta.tsv
# We therefore combine the two, and give every gene all categories in either file.
rule parse_midasdb_emapper_annotations_to_gene99_x_cog_category_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.emapper.gene99_x_cog_category.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_cog_category_v20.py",
        emapper="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/eggnog.tsv",
        cog_category="ref/cog-20.meta.tsv",
    shell:
        "{input.script} {input.emapper} {input.cog_category} {output}"


ruleorder: parse_midasdb_emapper_annotations_to_gene99_x_cog_category_v20 > parse_midasdb_emapper_annotations_to_gene99_x_unit_v20


# NOTE: Genomad annotations are done at the unique gene level.
rule parse_genomad_annotations_to_gene_x_accession_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.gene_x_genomad_{annot}.tsv",
    input:
        script="scripts/parse_genomad_annotations_to_gene_x_accession_v20.py",
        annot="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/genomad_{annot}.tsv",
    shell:
        "{input.script} {input.annot} {output}"


# NOTE: resfinder annotations are done at the unique gene level.
rule parse_resfinder_annotations_to_gene_x_accession_v20:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gene_x_amr.tsv",
    input:
        script="scripts/parse_resfinder_annotations_to_gene_x_accession.py",
        annot="ref/midasdb_uhgg_v20/pangenomes/{species}/annotation/resfinder.tsv",
    shell:
        "{input.script} {input.annot} {output}"


# Take annotations at the c99 level and combine them into centroidNN
# annotations using voting at the unique feature (gene_id) level.
rule aggregate_gene99_annotations_to_higher_centroid_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.emapper.gene{centroid}_x_{unit}.tsv",
    input:
        script="scripts/aggregate_gene99_annotations_to_higher_centroid.py",
        annot="data/species/sp-{species}/midasdb_v20.emapper.gene99_x_{unit}.tsv",
        clust="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    wildcard_constraints:
        centroid="95|90|85|80|75",  # Note no 99; this would potential result in a circular dependency.
    params:
        agg=lambda w: {
            95: "centroid_95",
            90: "centroid_90",
            85: "centroid_85",
            80: "centroid_80",
            75: "centroid_75",
        }[int(w.centroid)],
    resources:
        walltime_hr=10,
    shell:
        "{input.script} {input.annot} {input.clust} {params.agg} {output}"


# Take annotations at the gene level and combine them into centroidNN
# annotations using voting at the unique feature (gene_id) level.
rule aggregate_gene_annotations_to_higher_centroid_v20:
    output:
        "data/species/sp-{species}/midasdb_v20.gene{centroid}_x_{unit}.tsv",
    input:
        script="scripts/aggregate_gene_annotations_to_higher_centroid.py",
        annot="data/species/sp-{species}/midasdb_v20.gene_x_{unit}.tsv",
        clust="ref/midasdb_uhgg_v20/pangenomes/{species}/genes_info.tsv",
    params:
        agg=lambda w: {
            99: "centroid_99",
            95: "centroid_95",
            90: "centroid_90",
            85: "centroid_85",
            80: "centroid_80",
            75: "centroid_75",
        }[int(w.centroid)],
    resources:
        walltime_hr=10,
    shell:
        "{input.script} {input.annot} {input.clust} {params.agg} {output}"


rule parse_strain_emapper_annotations_to_gene_x_unit:
    output:
        "data/species/sp-{species}/genome/{genome}.prodigal-single.cds.emapper.gene_x_{unit}.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_{unit}.py",
        emapper="data/species/sp-{species}/genome/{genome}.prodigal-single.cds.emapper.d/proteins.emapper.annotations",
    shell:
        "{input.script} {input.emapper} {output}"


rule aggregate_strain_emapper_output_by_unit:
    output:
        "data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.{unit}-strain_gene.tsv",
    input:
        script="scripts/aggregate_emapper_output_by_unit.py",
        data="data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.gene_x_{unit}.tsv",
    shell:
        "{input.script} {input.data} {wildcards.unit} {wildcards.strain} {output}"


rule aggregate_midasdb_reference_gene_by_emapper_annotation:
    output:
        "data/species/sp-{species}/{stemB}.gene{centroid}_{dbv}.{unit}-strain_gene.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko|amr",
    input:
        script="scripts/aggregate_uhgg_strain_gene_by_annotation.py",
        uhgg="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
        mapping="data/species/sp-{species}/midasdb_{dbv}.emapper.gene{centroid}_x_{unit}.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.uhgg} {input.mapping} {wildcards.unit} {output}"


rule aggregate_midasdb_reference_gene_by_amr_annotation:
    output:
        "data/species/sp-{species}/{stemB}.gene{centroid}_{dbv}.amr-strain_gene.tsv",
    input:
        script="scripts/aggregate_uhgg_strain_gene_by_annotation.py",
        uhgg="data/species/sp-{species}/midasdb.gene{centroid}_{dbv}.uhgg-strain_gene.tsv",
        mapping="data/species/sp-{species}/midasdb_{dbv}.gene{centroid}_x_amr.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.uhgg} {input.mapping} accession_no {output}"


rule dbCAN_annotate_translated_orfs:
    output:
        dir=directory("{stem}.dbcan.d"),
    input:
        fasta="{stem}.tran.fa",
        db="ref/dbcan",
    conda:
        "conda/dbcan.yaml"
    threads: 4
    resources:
        walltime_hr=24,
        mem_mb=20_000,
        pmem=20_000 // 4,
    shell:
        """
        run_dbcan \
                {input.fasta} \
                protein \
                --db_dir {input.db} \
                --tools hmmer diamond \
                --tf_cpu {threads} --stp_cpu {threads} --dia_cpu {threads} --hmm_cpu {threads} --dbcan_thread {threads} \
                --out_dir {output.dir}
        """


rule aggregate_uhgg_strain_gene_by_emapper_annotation:
    output:
        "data/{stemA}/species/sp-{species}/{stemB}.gene{centroidA}_{dbv}-{pang_params}-agg{centroidB}.{stemC}.{unit}-strain_gene.tsv",
    wildcard_constraints:
        unit="eggnog|top_eggnog|cog|ko",
        centroidB="95|90|85|90|75",
    input:
        script="scripts/aggregate_uhgg_strain_gene_by_annotation.py",
        uhgg="data/{stemA}/species/sp-{species}/{stemB}.gene{centroidA}_{dbv}-{pang_params}-agg{centroidB}.{stemC}.uhgg-strain_gene.tsv",
        mapping="data/species/sp-{species}/midasdb_{dbv}.emapper.gene{centroidB}_x_{unit}.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.uhgg} {input.mapping} {wildcards.unit} {output}"


rule aggregate_uhgg_strain_gene_by_amr_annotation:
    output:
        "data/{stemA}/species/sp-{species}/{stemB}.gene{centroidA}_{dbv}-{pang_params}-agg{centroidB}.{stemC}.amr-strain_gene.tsv",
    wildcard_constraints:
        centroidB="95|90|85|90|75",
    input:
        script="scripts/aggregate_uhgg_strain_gene_by_annotation.py",
        uhgg="data/{stemA}/species/sp-{species}/{stemB}.gene{centroidA}_{dbv}-{pang_params}-agg{centroidB}.{stemC}.uhgg-strain_gene.tsv",
        mapping="data/species/sp-{species}/midasdb_{dbv}.gene{centroidB}_x_amr.tsv",
    group:
        "assess_gene_inference_benchmark"
    shell:
        "{input.script} {input.uhgg} {input.mapping} accession_no {output}"
