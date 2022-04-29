# rule start_ipython_sfacts:
#     params:
#         sfacts_dev_path=config["software-dev-path"]["sfacts"],
#     container:
#         config["container"]["sfacts"]
#     shell:
#         """
#         export PYTHONPATH="{params.sfacts_dev_path}"
#         ipython
#         """


# Same as sfacts, but with an environment defined by conda/sfacts.yaml
# rather than hard-coded into the container.
# This causes problems with GPU sfacts.
rule start_ipython_sfacts:
    params:
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    container:
        config["container"]["mambaforge"]
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        ipython
        """


# rule start_jupyter_sfacts:
#     threads: MAX_THREADS
#     params:
#         sfacts_dev_path=config["software-dev-path"]["sfacts"],
#         port=config["jupyter_port"],
#     container:
#         config["container"]["sfacts"]
#     shell:
#         """
#         export PYTHONPATH="{params.sfacts_dev_path}"
#         jupyter lab --port={params.port}
#         """


rule start_jupyter_sfacts:
    threads: MAX_THREADS
    params:
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
        port=config["jupyter_port"],
    container:
        config["container"]["mambaforge"]
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        jupyter lab --port={params.port}
        """


rule load_metagenotype_from_merged_gtpro:
    output:
        "{stem}.gtpro.mgen.nc",
    input:
        "{stem}.gtpro_combine.tsv.bz2",
    params:
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    container:
        config["container"]["sfacts"]
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        python3 -m sfacts load --gtpro-metagenotype {input} {output}
        """


rule filter_metagenotype:
    output:
        "data/{stem}.filt-poly{poly}-cvrg{cvrg}.mgen.nc",
    input:
        "data/{stem}.mgen.nc",
    wildcard_constraints:
        poly="[0-9]+",
        cvrg="[0-9]+",
    params:
        poly=lambda w: float(w.poly) / 100,
        cvrg=lambda w: float(w.cvrg) / 100,
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    container:
        config["container"]["sfacts"]
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        python3 -m sfacts filter_mgen --min-minor-allele-freq {params.poly} --min-horizontal-cvrg {params.cvrg} {input} {output}
        """





rule fit_sfacts_strategy1:
    output:
        "{stem}.fit-sfacts1-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    input:
        "{stem}.mgen.nc",
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
        nposition="[0-9]+",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        rho_hyper=1.0,
        gamma_hyper=1e-10,
        pi_hyper=0.3,
        seed=lambda w: int(w.seed),
        model_name="model2",
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    container:
        config["container"]["mambaforge"]
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --no-nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} pi_hyper={params.pi_hyper} rho_hyper={params.rho_hyper} \
                --anneal-hyperparameters pi_hyper=1.0 \
                --anneal-wait 3000 --anneal-steps 6000 \
                -- {input} {output}
        """

rule fit_sfacts_strategy2:
    output:
        "{stem}.fit-sfacts2-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    input:
        "{stem}.mgen.nc",
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
        nposition="[0-9]+",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        rho_hyper=0.9,
        gamma_hyper=1e-10,
        pi_hyper=0.4,
        seed=lambda w: int(w.seed),
        model_name="model2",
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    container:
        config["container"]["mambaforge"]
    conda:
        "conda/sfacts.yaml"
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --no-nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                -- {input} {output}
        """
rule fit_sfacts_strategy_old:
    output:
        fit="{stem}.fit-sfacts-s{nstrain}-g{nposition}-seed{seed}.world.nc",
    wildcard_constraints:
        nstrain="[0-9]+",
        seed="[0-9]+",
        nposition="[0-9]+",
    input:
        data="{stem}.mgen.nc",
    params:
        nstrain=lambda w: int(w.nstrain),
        nposition=lambda w: int(w.nposition),
        precision=32,
        gamma_hyper=1e-6,
        rho_hyper=0.001,
        pi_hyper=0.3,
        seed=lambda w: int(w.seed),
        model_name="new_model_1",
        lag1=50,
        lag2=100,
        lr=0.05,
        min_learning_rate=1e-6,
        sfacts_dev_path=config["software-dev-path"]["sfacts"],
    resources:
        walltime_hr=36,
        pmem=5_000,
        mem_mb=5_000,
        device={0: "cpu", 1: "cuda"}[config["USE_CUDA"]],
        gpu_mem_mb={0: 0, 1: 5_000}[config["USE_CUDA"]],
    container:
        config["container"]["sfacts"]
    shell:
        """
        export PYTHONPATH="{params.sfacts_dev_path}"
        python3 -m sfacts fit -m {params.model_name}  \
                --verbose --device {resources.device} \
                --precision {params.precision} \
                --random-seed {params.seed} \
                --num-strains {params.nstrain} --num-positions {params.nposition} \
                --no-nmf-init \
                --hyperparameters gamma_hyper={params.gamma_hyper} \
                --hyperparameters pi_hyper={params.pi_hyper} \
                --hyperparameters rho_hyper={params.rho_hyper} \
                --optimizer-learning-rate {params.lr} \
                --min-optimizer-learning-rate {params.min_learning_rate} \
                --max-iter 1_000_000 --lag1 {params.lag1} --lag2 {params.lag2} \
                {input.data} \
                {output.fit}
        """