# Install Jupyter Kernels


rule install_jupyter_kernel_native:
    params:
        name="native",
    shell:
        """
        python -m ipykernel install --user --name={params.name} --env PATH $PATH
        """


rule install_jupyter_kernel_default:
    params:
        name="default",
    conda:
        "conda/default.yaml"
    shell:
        """
        python -m ipykernel install --user --name={params.name} --env PATH $PATH
        """


# Any conda environment spec can be installed for Jupyter use
# just like the below:
use rule install_jupyter_kernel_default as install_jupyter_kernel_pymc with:
    params:
        name="pymc",
    conda:
        "conda/pymc.yaml"


# And then run `snakemake -j1 install_jupyter_kernel_pymc`.


rule start_jupyter:
    params:
        port=config["jupyter_port"],
    shell:
        "jupyter lab --port={params.port} --notebook-dir nb/"


rule start_ipython:
    shell:
        "ipython"


rule start_shell:
    container:
        config["container"]["toolz"]
    shell:
        "bash"


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
