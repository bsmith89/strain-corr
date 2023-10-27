# Same as sfacts, but with an environment defined by conda/sfacts.yaml
# rather than hard-coded into the container.
# This causes problems with GPU sfacts.
rule compile_species_variation_from_vcf:
    output:
        "data/species/sp-{species}/gtpro_ref.mgtp.nc",
    input:
        script="scripts/vcf_to_sfacts.py",
        gtpro_snp_dict="ref/gtpro/variants_main.covered.hq.snp_dict.tsv",
        vcf="raw/gtpro_refs/variation_in_species/{species}/core_snps.vcf.gz",
    conda:
        "conda/sfacts.yaml"
    resources:
        mem_mb=10_000,
    shell:
        """
        {input.script} {input.gtpro_snp_dict} {input.vcf} {wildcards.species} {output}
        """


rule concatenate_gtpro_refs:
    output:
        "{stemA}/species/sp-{species}/{stemB}.plus_refs.mgtp.nc",
    input:
        script="scripts/stack_mgen_with_refs.py",
        mgen="{stemA}/species/sp-{species}/{stemB}.mgtp.nc",
        ref="data/species/sp-{species}/gtpro_ref.mgtp.nc",
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
        "{stem}.filt-poly{poly}-cvrg{cvrg}.mgtp.nc",
    input:
        "{stem}.mgtp.nc",
    wildcard_constraints:
        poly="[0-9]+",
        cvrg="[0-9]+",
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
    conda:
        "conda/sfacts.yaml"
    resources:
        mem_mb=12_000,
    shell:
        """
        python3 -m sfacts filter_mgen \
                --min-minor-allele-freq {params.poly} \
                --min-horizontal-cvrg {params.cvrg} \
                {input} {output}
        """


rule subset_metagenotype:
    output:
        "{stem}.ss-g{num_positions}-block{block_number}-seed{seed}.mgtp.nc",
    input:
        "{stem}.mgtp.nc",
    wildcard_constraints:
        num_positions=integer_wc,
        block_number=integer_wc,
        seed=integer_wc,
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


rule resample_metagenotype:
    output:
        "{stem}.rs-alpha{alpha}-g{num_positions}-seed{seed}.mgtp.nc",
    input:
        "{stem}.mgtp.nc",
    params:
        seed=lambda w: int(w.seed),
        num_positions=lambda w: int(w.num_positions),
        entropy_weight_alpha=lambda w: float(w.alpha),
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts sample_mgen \
                --random-seed {params.seed} \
                --with-replacement \
                --num-positions {params.num_positions} \
                --entropy-weighted-alpha {params.entropy_weight_alpha} \
                {input} \
                {output}
        """


rule sfacts_nmf_approximation:
    output:
        "{stem}.approx-nmf-s{strain_exponent}-seed{seed}.world.nc",
    input:
        "{stem}.mgtp.nc",
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


# FIXME: Seed parameterized filename is not really a parameter.
rule sfacts_clust_initialization:
    output:
        "{stem}.approx-clust-thresh{thresh}-s{strain_exponent}-seed0.world.nc",
    input:
        "{stem}.mgtp.nc",
    params:
        thresh=lambda w: float(w.thresh) / 100,
        strain_exponent=lambda w: float(w.strain_exponent) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts clust_init \
                --verbose \
                --strain-sample-exponent {params.strain_exponent} \
                --thresh {params.thresh} \
                {input} \
                {output}
        """


rule sfacts_clust_approximation:
    output:
        "{stem}.approx-clust2-thresh{thresh}-s{strain_exponent}.world.nc",
    input:
        "{stem}.mgtp.nc",
    params:
        thresh=lambda w: float(w.thresh) / 100,
        strain_exponent=lambda w: float(w.strain_exponent) / 100,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts clust_init \
                --verbose \
                --strain-sample-exponent {params.strain_exponent} \
                --thresh {params.thresh} \
                --frac 1.0 \
                --pseudo 0.0 \
                {input} \
                {output}
        """


rule calculate_metagenotype_pdist:
    output:
        "{stem}.mgtp.pdist.nc",
    input:
        mgen="{stem}.mgtp.nc",
    conda:
        "conda/sfacts.yaml"
    resources:
        walltime_hr=5,
    shell:
        "sfacts mgen_diss --verbose {input.mgen} {output}"


rule fit_sfacts:
    output:
        fit="{stem}.fit-sfacts{strategy}-s{strain_exponent}-seed{seed}.world.nc",
        hist="{stem}.fit-sfacts{strategy}-s{strain_exponent}-seed{seed}.loss_history",
    input:
        mgen="{stem}.mgtp.nc",
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
        gpu=1,
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


rule cleanup_sfacts_fit:
    output:
        "{stem}.clean-m{m}-e{e}.world.nc",
    input:
        "{stem}.world.nc",
    params:
        metagenotype_error_thresh=lambda w: int(w.m) / 100,
        entropy_error_thresh=lambda w: int(w.e) / 100,
        monte_carlo_draws=10,
        random_seed=0,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        sfacts cleanup_fit2 \
                --metagenotype-error {params.metagenotype_error_thresh} \
                --entropy-error {params.entropy_error_thresh} \
                --monte-carlo-draws {params.monte_carlo_draws} \
                --random-seed {params.random_seed} \
                {input} {output}
        """


rule refit_genotypes_sfacts:
    output:
        fit="data/{stemA}.fit-{stemB}.refit-sfacts{strategy}-seed{seed}.world.nc",
    input:
        fit="data/{stemA}.fit-{stemB}.world.nc",
        mgen="data/{stemA}.mgtp.nc",
        strategy="meta/sfacts/strategy{strategy}.args",
    wildcard_constraints:
        seed="[0-9]+",
    params:
        seed=lambda w: int(w.seed),
    resources:
        walltime_hr=2,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
        gpu=1,
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        python3 -m sfacts fit_geno \
                @{input.strategy} \
                --debug --device {resources.device} \
                --random-seed {params.seed} \
                -- {input.fit} {input.mgen} {output.fit}
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


rule alias_canonical_fit:
    output:
        "data/group/{group}/species/sp-{species}/r.{proc}.gtpro.sfacts-fit.world.nc",
    input:
        source=lambda w: (
            "data/group/{group}/species/sp-{species}/r.{proc}.gtpro.{sfacts_stem}.world.nc".format(
                group=w.group,
                species=w.species,
                proc=w.proc,
                sfacts_stem=config["species_group_to_sfacts_stem"][
                    (w.species, w.group)
                ],
            )
        ),
    shell:
        alias_recipe


localrules:
    alias_canonical_fit,


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


rule sfacts_metagenotypet_to_tsv:
    output:
        "{stem}.mgtp.tsv",
    input:
        "{stem}.mgtp.nc",
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        sfacts dump --metagenotype {output} {input}
        """


# NOTE: Hub-rule: Comment out this rule to reduce DAG-building time
# once it has been run for the focal group.
rule calculate_all_strain_depths:
    output:
        "data/group/{group}/r.{proc}.gtpro.{stem}.strain_depth.tsv",
    input:
        script="scripts/merge_all_strains_depth.py",
        species="data/group/{group}/r.{proc}.gtpro.species_depth.tsv",
        strains=lambda w: [
            f"data/group/{w.group}/species/sp-{species}/r.{w.proc}.gtpro.{w.stem}.comm.tsv"
            for species in config["species_group"][w.group]
        ],
    params:
        args=lambda w: [
            f"{species}=data/group/{w.group}/species/sp-{species}/r.{w.proc}.gtpro.{w.stem}.comm.tsv"
            for species in config["species_group"][w.group]
        ],
    shell:
        """
        {input.script} {input.species} {output} {params.args}
        """
