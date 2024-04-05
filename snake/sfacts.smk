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


rule calculate_all_strain_depths:  # Hub-rule
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
