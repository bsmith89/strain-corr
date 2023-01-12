rule patch_sfinder:
    output:
        touch("build/sfinder_patched.flag"),
    input:
        patch="include/StrainFinder.patch",
        repo="include/StrainFinder",
    shell:
        """
        cd {input.repo}
        git apply ../StrainFinder.patch
        """


use rule start_jupyter as start_jupyter_sfinder with:
    conda:
        "conda/sfinder.yaml"


use rule start_ipython as start_ipython_sfinder with:
    conda:
        "conda/sfinder.yaml"


use rule start_shell as start_shell_sfinder with:
    conda:
        "conda/sfinder.yaml"


rule metagenotype_tsv_to_sfinder_aln:
    output:
        cpickle="{stem}.sfinder.aln.cpickle",
        indexes="{stem}.sfinder.aln.indexes.txt",
    input:
        script="scripts/metagenotype_to_sfinder_alignment.py",
        data="{stem}.mgtp.tsv",
    conda:
        "conda/sfinder.yaml"
    shell:
        """
        {input.script} {input.data} {output}
        """


localrules:
    metagenotype_tsv_to_sfinder_aln,


rule fit_sfinder:
    output:
        "{stem}.fit-sfinder-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
        max_runtime_s=72000,
    resources:
        walltime_sec=72000,
    shell:
        """
        rm -rf {output}
        include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 1 --dtol 1 --ntol 2 --max_time {resources.walltime_sec} --n_keep 5 --converge \
                --em_out {output}
        """


rule fit_sfinder_exhaustive:
    output:
        "{stem}.fit-sfinder2-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
        max_runtime_s=72000,
    resources:
        walltime_sec=72000,
    shell:
        """
        rm -rf {output}
        include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --exhaustive \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 1 --dtol 1 --ntol 2 --max_time {resources.walltime_sec} --n_keep 5 --converge \
                --em_out {output}
        """


rule fit_sfinder_global:
    output:
        "{stem}.fit-sfinder3-s{nstrain}-seed{seed}.em.cpickle",
    input:
        "{stem}.sfinder.aln.cpickle",
    conda:
        "conda/sfinder.yaml"
    params:
        nstrain=lambda w: int(w.nstrain),
        seed=lambda w: int(w.seed),
    shell:
        """
        rm -rf {output}
        include/StrainFinder/StrainFinder.py \
                --force_update --merge_out --msg \
                --aln {input} \
                --seed {params.seed} \
                -N {params.nstrain} \
                --max_reps 10 --dtol 1 --ntol 2 --max_time 7200 --n_keep 3 --converge \
                --em_out {output}
        """


rule parse_sfinder_cpickle:
    output:
        comm="{stem}.fit-sfinder{params}.comm.tsv",
        geno="{stem}.fit-sfinder{params}.geno.tsv",
    input:
        script="scripts/strainfinder_result_to_flatfiles.py",
        cpickle="{stem}.fit-sfinder{params}.em.cpickle",
        indexes="{stem}.sfinder.aln.indexes.txt",
    conda:
        "conda/sfinder.yaml"
    shell:
        """
        {input.script} {input.cpickle} {input.indexes} {output.comm} {output.geno}
        """


localrules:
    parse_sfinder_cpickle,
