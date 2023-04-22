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


rule run_panphlan:
    output:
        tally="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_map.csv",
        hit="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_hit.tsv",
        depth="data/species/sp-{species}/reads/{mgen}/r.{stem}.panphlan_depth.tsv",
    input:
        pangenome="ref/panphlan/{species}",
        r1="data/reads/{mgen}/r1.{stem}.fq.gz",
        r2="data/reads/{mgen}/r2.{stem}.fq.gz",
    params:
        panphlan_species=lambda w: config["species_to_panphlan"][w.species],
    conda:
        "conda/panphlan.yaml"
    threads: 2
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
            -o {output.tally}
        rm $fastq
        panphlan_profiling.py -i <(echo {output.tally}) \
                -p {input.pangenome}/{params.panphlan_species}_pangenome.tsv \
                --o_matrix {output.hit} \
                --o_covmat {output.depth}
        """
