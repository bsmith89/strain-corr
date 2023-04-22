use rule start_shell as start_shell_panphlan with:
    conda:
        "conda/panphlan.yaml"


rule download_panphlan_reference:
    output:
        directory("raw/ref/panphlan/{panphlan_species}"),
    params:
        outdir="raw/ref/panphlan",
    conda:
        "conda/panphlan.yaml"
    shell:
        """
        panphlan_download_pangenome.py -i {wildcards.panphlan_species} -o {params.outdir}
        """


rule link_panphlan_reference:
    output:
        directory("ref/panphlan/{species}"),
    input:
        dir=lambda w: ancient(
            "raw/ref/panphlan/" + config["species_to_panphlan"][w.species]
        ),
    shell:
        alias_recipe


rule run_panphlan_map:
    output:
        "data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_map.csv",
    input:
        pangenome="ref/panphlan/{species}",
        r1="data/reads/{mgen}/r1.{stem}.fq.gz",
        r2="data/reads/{mgen}/r2.{stem}.fq.gz",
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
    conda:
        "conda/panphlan.yaml"
    threads: 8
    resources:
        walltime_hr=12,
    shell:
        """
        fastq=$(mktemp)
        echo Unzipping {input.r1} and {input.r2} to $fastq
        gzip -dc {input.r1} {input.r2} > $fastq
        echo Running panphlan_map on {input.r1} and {input.r2} as $fastq
        panphlan_map.py --nproc {threads} \
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --indexes {input.pangenome}/{params.panphlan_species} \
                -i $fastq \
            -o {output}
        rm $fastq
        """


rule run_panphlan_profile:
    output:
        hit="data/group/{group}/species/sp-{species}/{stem}.panphlan_hit.tsv",
        cvrg="data/group/{group}/species/sp-{species}/{stem}.panphlan_cvrg.tsv",
    input:
        samples=lambda w: [
            f"data/species/sp-{w.species}/reads/{mgen}/{w.stem}.panphlan_map.csv"
            for mgen in config["mgen_group"][w.group]
        ],
        pangenome="ref/panphlan/{species}",
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
        sample_list=lambda w: config["mgen_group"][w.group],
        sample_pattern="data/species/sp-{species}/reads/$sample/{stem}.panphlan_map.csv",
    conda:
        "conda/panphlan.yaml"
    shell:
        """
        tmpdir=$(mktemp -d)
        for sample in {params.sample_list}
        do
            ln -s $(realpath {params.sample_pattern}) $tmpdir/$sample
        done
        ls $tmpdir
        panphlan_profiling.py -i $tmpdir \
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --o_matrix {output.hit} \
                --o_covmat {output.cvrg}
        """
                # --func_annot {input.pangenome}/panphlan_{params.panphlan_species}_annot.tsv \
                # --field 8 \
