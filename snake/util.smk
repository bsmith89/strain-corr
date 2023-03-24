# Install Jupyter Kernels


rule install_jupyter_kernel_default:
    params:
        name="default",
    conda:
        "conda/toolz.yaml"
    shell:
        """
        python -m ipykernel install --user --name={params.name} --env PATH $PATH
        """


rule start_jupyter:
    params:
        port=config["jupyter_port"],
    shell:
        "jupyter lab --port={params.port} --notebook-dir nb/"


rule start_ipython:
    shell:
        "ipython"


use rule start_ipython as start_ipython_toolz with:
    conda:
        "conda/toolz.yaml"


rule start_shell:
    shell:
        "bash"


use rule start_shell as start_shell_toolz with:
    conda:
        "conda/toolz.yaml"


use rule start_shell as start_shell_seqtk with:
    conda:
        "conda/seqtk.yaml"


use rule start_shell as start_shell_hsblastn with:
    conda:
        "conda/hsblastn.yaml"


rule visualize_rulegraph:
    output:
        "data/rulegraph.dot",
    input:
        "Snakefile",
    shell:
        dd(
            """
        snakemake --rulegraph all > {output}
        """
        )


rule generate_report:
    output:
        "fig/report.html",
    input:
        "Snakefile",
    shell:
        dd(
            """
        snakemake --forceall --report {output} all
        """
        )


rule dot_to_pdf:
    output:
        "fig/{stem}.pdf",
    input:
        "data/{stem}.dot",
    shell:
        dd(
            """
        dot -Tpdf < {input} > {output}
        """
        )


rule processed_notebook_to_html:
    output:
        "build/{stem}.ipynb.html",
    input:
        "build/{stem}.ipynb",
    shell:
        dd(
            """
        jupyter nbconvert -t html {input} {output}
        """
        )


rule serve_directory:
    params:
        port=config["server_port"],
    shell:
        """
        python3 -m http.server {params.port}
        """


rule query_db:
    output:
        "data/{db}.select_{query}.tsv",
    input:
        db="data/{db}.db",
        query="scripts/query/{query}.sql",
    shell:
        dd(
            """
        sqlite3 -header -separator '\t' {input.db} < {input.query} > {output}
        """
        )


rule config_debug:
    output:
        "config_debug.{config_key}",
    params:
        meta=lambda w: nested_dictlookup(config, *w.config_key.split(".")),
    shell:
        """
        echo "{wildcards.config_key}"
        echo "{params.meta}"
        false  # Recipe never succeeds.
        """
