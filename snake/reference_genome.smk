rule link_reference_genome:
    output:
        "data/species/sp-{species}/genome/{genome}.fn",
    input:
        lambda w: config["genome"].loc[w.genome].genome_path,
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe

rule link_midasdb_reference_genome:
    output:
        "data/species/sp-{species}/genome/midasdb_{dbv}/{genome}.fn",
    input: "ref/midasdb_uhgg_{dbv}/mags/{species}/{genome}.fa"
    wildcard_constraints:
        genome=noperiod_wc,
    shell:
        alias_recipe


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


rule parse_emapper_output_to_gene_x_unit:
    output:
        "{stem}.emapper.gene_x_{unit}.tsv",
    input:
        script="scripts/parse_emapper_output_to_gene_x_{unit}.py",
        emapper="{stem}.emapper.d/proteins.emapper.annotations",
    shell:
        "{input.script} {input.emapper} {output}"


rule aggregate_strain_emapper_output_by_unit:
    output:
        "data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.{agg}-strain_gene.tsv",
    input:
        "data/species/sp-{species}/genome/{strain}.prodigal-single.cds.emapper.gene_x_{agg}.tsv",
    run:
        strain_gene_x_agg = pd.read_table(input[0])
        result = (
            (strain_gene_x_agg.groupby(wildcards.agg).apply(len) > 0)
            .astype(int)
            .to_frame(name=wildcards.strain)
        ).rename_axis(index="gene_id")
        result.to_csv(output[0], sep="\t")


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


rule normalize_genome_sequence:
    output:
        "{stem}.norm.fn",
    input:
        "{stem}.fn",
    shell:
        "sed '/^>/!s/[^ACGT]/N/g' {input} > {output}"


rule tile_reference_genome:
    output:
        "{stem}.tiles-l{length}-o{overlap}.fn",
    input:
        script="scripts/tile_fasta.py",
        fn="{stem}.fn",
    wildcard_constraints:
        genome=noperiod_wc,
        length=integer_wc,
        overlap=integer_wc,
    params:
        length=lambda w: int(w.length),
        overlap=lambda w: int(w.overlap),
    shell:
        "{input.script} {params.length} {params.overlap} {input.fn} > {output}"


rule genome_fasta_to_fastq:
    """
    Convert a FASTA formatted file into FASTQ.
    \
    Input/output patterns are limited to files found in */genome/* in order to
    prevent circular dependencies.
    \
    """
    output:
        "{stemA}/genome/{stemB}.fq.gz",
    input:
        "{stemA}/genome/{stemB}.fn",
    conda:
        "conda/seqtk.yaml"
    shell:
        "seqtk seq -F '#' {input} | gzip -c > {output}"


rule combine_strain_genome_gtpro_data_loadable:
    output:
        "data/group/{group}/species/sp-{species}/strain_genomes.gtpro.tsv.bz2",
    input:
        strain=lambda w: [
            f"data/species/sp-{w.species}/genome/{strain}.norm.tiles-l500-o31.gtpro_parse.tsv.bz2"
            for strain in species_group_genomes(w.species, w.group)
        ],
    params:
        strain_list=lambda w: species_group_genomes(w.species, w.group),
        pattern="data/species/sp-{species}/genome/$strain.norm.tiles-l500-o31.gtpro_parse.tsv.bz2",
    shell:
        """
        for strain in {params.strain_list}
        do
            path={params.pattern}
            ( \
                bzip2 -dc $path \
                | awk -v OFS='\t' -v strain=$strain -v species={wildcards.species} '$1==species {{print strain,$0}}' \
                | bzip2 -zc \
            )
        done > {output}

        """


rule combine_midasdb_reference_genome_gtpro_data_loadable:
    output:
        "data/species/sp-{species}/midasdb_{dbv}.gtpro.tsv.bz2",
    input:
        genome=lambda w: [
            f"data/species/sp-{w.species}/genome/midasdb_{w.dbv}/{genome}.norm.tiles-l500-o31.gtpro_parse.tsv.bz2"
            for genome in config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species]
        ],
    params:
        genome_list=lambda w: config[f"midasdb_uhgg_{w.dbv}_species_genome"][w.species],
        pattern="data/species/sp-{species}/genome/midasdb_{dbv}/$genome.norm.tiles-l500-o31.gtpro_parse.tsv.bz2",
    shell:
        """
        for genome in {params.genome_list}
        do
            path={params.pattern}
            ( \
                bzip2 -dc $path \
                | awk -v OFS='\t' -v strain=$genome -v species={wildcards.species} '$1==species {{print strain,$0}}' \
                | bzip2 -zc \
            )
        done > {output}

        """


rule alias_midas_uhgg_pangenome_cds_new:
    output:
        "data/species/sp-{species}/pangenome_{dbv}.centroids.fn",
    input:
        fasta="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/centroids.ffn",
    shell:
        alias_recipe


localrules:
    alias_midas_uhgg_pangenome_cds_new,




rule blastn_genome_against_midasdb_uhgg_new:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.tsv",
    input:
        query="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
        subject="data/species/sp-{species}/pangenome_{dbv}.centroids.fn",
    threads: 1
    shell:
        """
        blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
        """


rule assign_matching_genes_based_on_best_blastn_hit:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best.tsv",
    input:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.tsv",
    shell:
        """
        sort -k1,1 -k12,12rn {input} | sort -k1,1 -u | cut -f1,2 > {output}
        """


rule aggreggate_top_blastn_hits_by_midasdb_centroid:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best-c{centroid}.uhggtop-strain_gene.tsv",
    input:
        script="scripts/identify_strain_genes_from_top_blastn_hits.py",
        hit="{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-blastn.gene_matching-best.tsv",
        gene_clust="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/cluster_info.txt",
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
        {input.script} {input.gene_clust} {params.aggregate_genes_by} {input.hit} {output}
        """


rule blastn_genome_against_genome:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.{stemB}-blastn.tsv",
    input:
        query="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
        subject="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.fn",
    threads: 1
    shell:
        """
        blastn -query {input.query} -subject {input.subject} -max_target_seqs 100000 -num_threads {threads} -outfmt 6 > {output}
        """


use rule diamond_search_fa as blastp_midasdb_uhgg_against_genome with:
    output:
        "{stemA}/species/sp-{species}/genome/midas_uhgg_pangenome.{stemB}-blastp.tsv",
    input:
        fasta="data/species/sp-{species}/pangenome.centroids.tran.fa",
        db="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.dmnd",


use rule diamond_search_fa as blastp_genome_against_midasdb_uhgg with:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome-blastp.tsv",
    input:
        fasta="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.fa",
        db="data/species/sp-{species}/pangenome.centroids.tran.dmnd",
    params:
        db=lambda w, input: parse_diamond_db_from_path(input.db),
        extra_diamond_blastp_args="--ultra-sensitive --max-hsps 10000",


use rule diamond_search_fa as reciprocal_blastp_genome with:
    output:
        "{stemA}/species/sp-{species}/genome/{stemB}.{stemB}-blastp.tsv",
    input:
        fasta="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.fa",
        db="{stemA}/species/sp-{species}/genome/{stemB}.prodigal-single.cds.tran.dmnd",


rule calculate_bitscore_ratio_of_orfs_and_pangenome_genes_new:
    output:
        "data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-{blastn_or_p}.bitscore_ratio-c{centroid}.tsv",
    input:
        script="scripts/calculate_bitscore_ratio_for_gene_matching.py",
        # NOTE: orf_x_orf is used to get the maximum possible bitscore for each ORF.
        orf_x_orf="data/species/sp-{species}/genome/{stemB}.{stemB}-{blastn_or_p}.tsv",
        orf_x_midas="data/species/sp-{species}/genome/{stemB}.midas_uhgg_pangenome_{dbv}-{blastn_or_p}.tsv",
        gene_clust="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/cluster_info.txt",
    wildcard_constraints:
        blastn_or_p="blastn|blastp",
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
                {input.orf_x_orf} \
                {input.orf_x_midas} \
                {input.gene_clust} \
                {params.aggregate_genes_by} \
                {output}
        """


rule assign_matching_genes_based_on_bitscore_ratio:
    output:
        "data/species/sp-{species}/genome/{stemB}.gene_matching-c{centroid}-t{thresh}.tsv",
    input:
        ratio="data/species/sp-{species}/genome/{stemB}.bitscore_ratio-c{centroid}.tsv",
    params:
        thresh=lambda w: int(w.thresh) / 100,
    shell:
        """
        awk -v OFS='\t' -v thresh='{params.thresh}' '(NR > 1) && ($3 > thresh) {{print $1,$2}}' {input} > {output}
        """


use rule run_bowtie_multispecies_pangenome_v22_new as run_bowtie_multispecies_pangenome_on_reference_genome_tiles_v22_new with:
    output:
        "data/hash/{hash}/species/sp-{species}/genome/{genome}.tiles-{tile_params}.pangenomes{centroid}_{dbv}-v22.{bam_or_cram}",
    input:
        db="data/hash/{hash}/pangenomes{centroid}_{dbv}.bt2.d/centroids.bt2db",
        r1="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",
        r2="data/species/sp-{species}/genome/{genome}.tiles-{tile_params}.fq.gz",
    threads: 1


# NOTE: Hub-rule paired with it's parent.
use rule load_one_species_pangenome2_depth_into_netcdf_new as load_one_species_pangenome2_tile_depth_into_netcdf_new with:
    output:
        "data/group/{group}/species/sp-{species}/ALL_STRAINS.{stem}.gene{centroidA}_{dbv}-{bowtie_params}-agg{centroidB}.depth2.nc",
    input:
        script="scripts/merge_pangenomes_depth.py",
        samples=lambda w: [
            "data/hash/{_hash}/species/sp-{w.species}/genome/{genome}.{w.stem}.pangenomes{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                genome=genome,
                _hash=config["species_group_to_hash"][w.group],
            )
            for genome in species_group_genomes(w.species, w.group)
        ],
        gene_info="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/gene_info.txt",
        gene_length="ref/midasdb_uhgg_{dbv}/pangenomes/{species}/genes.len",
    params:
        args=lambda w: [
            "{genome}=data/hash/{_hash}/species/sp-{w.species}/genome/{genome}.{w.stem}.pangenomes{w.centroidA}_{w.dbv}-{w.bowtie_params}.gene_mapping_tally.tsv.lz4".format(
                w=w,
                genome=genome,
                _hash=config["species_group_to_hash"][w.group],
            )
            for genome in species_group_genomes(w.species, w.group)
        ],
        centroidB_col=lambda w: f"centroid_{w.centroidB}",


# NOTE: So-as to re-use the "gene matching" instrument designed for
# blastn-based accuracy assessment, here we "label" each gene with arbitrary,
# ascending integers. (In other words, genes are no longer mapped to ORFs, so we
# use ascending integers as "dummy" ORF IDs, in order to keep the file format the same.)
rule assign_matching_genes_based_on_tile_depth:
    output:
        "{stemA}/species/sp-{species}/genome/{strain}.tiles-{tile_params}.gene{centroidA}-{params}-agg{centroidB}.gene_matching-t{thresh}.uhggtiles-strain_gene.tsv",
    input:
        script="scripts/identify_strain_genes_from_tiling_depth.py",
        depth="{stemA}/species/sp-{species}/ALL_STRAINS.tiles-{tile_params}.gene{centroidA}-{params}-agg{centroidB}.depth2.nc",
    params:
        thresh=lambda w: int(w.thresh),
    shell:
        """
        {input.script} {input.depth} {wildcards.strain} {params.thresh} {output}
        """
