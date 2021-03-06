# Same as sfacts, but with an environment defined by conda/sfacts.yaml
# rather than hard-coded into the container.
# This causes problems with GPU sfacts.
rule start_ipython_sfacts:
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        ipython
        """


rule start_shell_sfacts:
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        bash
        """


# NOTE: Comment out this rule to speed up DAG-building time
rule load_metagenotype_from_merged_gtpro:
    output:
        "{stem}.gtpro.mgen.nc",
    input:
        "{stem}.gtpro.tsv.bz2",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts load --gtpro-metagenotype {input} {output}
        """


rule compile_species_variation_from_vcf:
    output:
        "data/sp-{species}.gtpro_ref.mgen.nc",
    input:
        script="scripts/vcf_to_sfacts.py",
        gtpro_snp_dict="ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
        vcf="raw/gtpro_refs/variation_in_species/{species}/core_snps.vcf.gz",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.gtpro_snp_dict} {input.vcf} {wildcards.species} {output}
        """


rule concatenate_gtpro_refs:
    output:
        "data/sp-{species}.{group}.a.{stem}.plus_refs.mgen.nc",
    input:
        script="scripts/stack_mgen_with_refs.py",
        mgen="data/sp-{species}.{group}.a.{stem}.mgen.nc",
        ref="data/sp-{species}.gtpro_ref.mgen.nc",
    params:
        multi=100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.mgen} {input.ref} {params.multi} {output}
        """


rule filter_metagenotype:
    output:
        "{stem}.filt-poly{poly}-cvrg{cvrg}.mgen.nc",
    input:
        "{stem}.mgen.nc",
    wildcard_constraints:
        poly="[0-9]+",
        cvrg="[0-9]+",
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts filter_mgen \
                --min-minor-allele-freq {params.poly} \
                --min-horizontal-cvrg {params.cvrg} \
                {input} {output}
        """


rule subset_metagenotype:
    output:
        "{stem}.ss-g{num_positions}-block{block_number}-seed{seed}.mgen.nc",
    input:
        "{stem}.mgen.nc",
    params:
        seed=lambda w: int(w.seed),
        num_positions=lambda w: int(w.num_positions),
        block_number=lambda w: int(w.block_number),
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts sample_mgen \
                --random-seed {params.seed} \
                --num-positions {params.num_positions} \
                --block-number {params.block_number} \
                {input} \
                {output}
        """


rule sfacts_nmf_approximation:
    output:
        "{stem}.approx-nmf-s{strain_exponent}-seed{seed}.world.nc",
    input:
        "{stem}.mgen.nc",
    params:
        seed=lambda w: int(w.seed),
        strain_exponent=lambda w: float(w.strain_exponent) / 100,
        alpha_genotype=0.1,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts nmf_init \
                --verbose \
                --random-seed {params.seed} \
                --strain-sample-exponent {params.strain_exponent} \
                --alpha-genotype {params.alpha_genotype} \
                {input} \
                {output}
        """


rule sfacts_clust_approximation:
    output:
        "{stem}.approx-clust-thresh{thresh}-s{strain_exponent}-seed{seed}.world.nc",
    input:
        "{stem}.mgen.nc",
    params:
        seed=lambda w: int(w.seed),
        thresh=lambda w: float(w.thresh) / 100,
        strain_exponent=lambda w: float(w.strain_exponent) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts clust_init \
                --verbose \
                --random-seed {params.seed} \
                --strain-sample-exponent {params.strain_exponent} \
                --thresh {params.thresh} \
                {input} \
                {output}
        """


rule fit_sfacts:
    output:
        fit="{stem}.fit-sfacts{strategy}-s{strain_exponent}-seed{seed}.world.nc",
        hist="{stem}.fit-sfacts{strategy}-s{strain_exponent}-seed{seed}.loss_history",
    input:
        mgen="{stem}.mgen.nc",
        strategy="meta/sfacts/strategy{strategy}.args",
    wildcard_constraints:
        strain_exponent="[0-9]+",
        nposition="[0-9]+",
        seed="[0-9]+",
    params:
        strain_exponent=lambda w: float(w.strain_exponent) / 100,
        seed=lambda w: int(w.seed),
    resources:
        walltime_hr=2,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts fit \
                @{input.strategy} \
                --verbose --device {resources.device} \
                --random-seed {params.seed} \
                --strain-sample-exponent {params.strain_exponent} \
                --history-outpath {output.hist} \
                -- {input.mgen} {output.fit}
        """


rule collapse_similar_strains:
    output:
        "{stem}.collapse-{diss}.world.nc",
    input:
        "{stem}.world.nc",
    params:
        diss=lambda w: float(w.diss) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        sfacts cleanup_fit --dissimilarity {params.diss} --discretized {input} {output}
        """


rule assign_plurality_strain:
    output:
        "{stem}.top_strain_only.world.nc",
    input:
        script="scripts/assign_plurality_strain.py",
        data="{stem}.world.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        {input.script} {input.data} {output}
        """


rule cleanup_fit:
    output:
        "{stem}.clean-diss{diss}-abund{abund}-entr{entr}.world.nc",
    wildcard_constraints:
        diss=noperiod_wc,
        abund=noperiod_wc,
        entr=noperiod_wc,
    input:
        "{stem}.world.nc",
    params:
        diss=lambda w: float(w.diss) / 100,
        abund=lambda w: float(w.abund) / 100,
        entr=lambda w: float(w.entr) / 10,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        sfacts cleanup_fit --dissimilarity {params.diss} --abundance {params.abund} --entropy {params.entr} {input} {output}
        """


rule export_sfacts_comm:
    output:
        "{stem}.comm.tsv",
    input:
        "{stem}.world.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        sfacts dump --community {output} {input}
        """


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group.
rule calculate_all_strain_depths:
    output:
        "data_temp/{group}.a.{proc_stem}.gtpro.{fit_stem}.strain_depth.tsv",
    input:
        script="scripts/merge_all_strains_depth.py",
        species="data/{group}.a.{proc_stem}.gtpro.species_depth.tsv",
        strains=lambda w: [
            f"data_temp/sp-{species}.{w.group}.a.{w.proc_stem}.gtpro.{w.fit_stem}.comm.tsv"
            for species in checkpoint_select_species_with_greater_max_coverage_gtpro(
                group=w.group,
                stem=w.proc_stem,
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    params:
        args=lambda w: [
            f"{species}=data_temp/sp-{species}.{w.group}.a.{w.proc_stem}.gtpro.{w.fit_stem}.comm.tsv"
            for species in checkpoint_select_species_with_greater_max_coverage_gtpro(
                group=w.group,
                stem=w.proc_stem,
                cvrg_thresh=0.2,
                num_samples=2,
                require_in_species_group=True,
            )
        ],
    shell:
        """
        {input.script} {input.species} {output} {params.args}
        """
