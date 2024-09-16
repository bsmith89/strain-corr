rule sort_bib_from_raw:
    output:
        "doc/bibliography.bib",
    input:
        script="scripts/sort_bib.py",
        bib=["doc/bibliography_raw.bib"],
    conda: "conda/pandoc.yaml"
    shell:
        "python {input.script} {input.bib} > {output}"

rule render_pdf_to_png_imagemagick:
    output:
        "fig/{stem}.dpi{dpi}.png",
    input:
        "fig/{stem}.pdf",
    params:
        dpi=lambda w: int(w.dpi),
    conda: "conda/pandoc.yaml"
    shell:
        """
        convert -units PixelsPerInch -density {params.dpi} {input} {output}
        """


rule _temp_test_biblib:
    conda: "conda/pandoc.yaml"
    shell: 'python3 -c "import biblib"'

# rule render_figure_to_png:
#     output:
#         "fig/{stem}_figure.w{width}.png",
#     input:
#         "doc/static/{stem}_figure.svg",
#     params:
#         width=lambda w: int(w.width),
#     shell:
#         """
#         inkscape {input} --export-width={params.width} --export-filename {output}
#         """


rule render_figure_to_pdf:
    output:
        "fig/{stem}_figure.pdf",
    input:
        "doc/static/{stem}_figure.svg",
    conda: "conda/pandoc.yaml"
    shell:
        """
        inkscape {input} --export-filename {output}
        """



#
#
# rule render_pdf_to_tiff_imagemagick:
#     output:
#         "fig/{stem}.dpi{dpi}.tiff",
#     input:
#         "fig/{stem}.pdf",
#     params:
#         dpi=lambda w: int(w.dpi),
#     shell:
#         """
#         convert -units PixelsPerInch -density {params.dpi} {input} {output}
#         """
#
#
# rule link_static_pdf_figure:
#     output:
#         "fig/{stem}_figure.pdf",
#     input:
#         "doc/static/{stem}_figure.pdf",
#     shell:
#         alias_recipe
#
#
# rule pdf_to_eps:
#     output:
#         "fig/{stem}.eps",
#     input:
#         "fig/{stem}.pdf",
#     shell:
#         """
#         cd fig
#         pdf2ps {wildcards.stem}.pdf
#         ps2eps {wildcards.stem}.ps
#         rm {wildcards.stem}.ps
#         """


rule build_manuscript_docx:
    output:
        "build/manuscript.docx",
    input:
        source="doc/manuscript/manuscript.md",
        bib="doc/bibliography.bib",
        template="doc/manuscript/example_style.docx",
        csl="doc/manuscript/citestyle.csl",
        figures=[
            "fig/concept_diagram_figure.dpi200.png",
            "fig/benchmarking_figure.dpi200.png",
            "fig/hmp2_diversity_figure.dpi200.png",
            "fig/pangenomics_figure.dpi200.png",
            "fig/ucfmt_figure.dpi200.png",
            "fig/accuracy_by_depth_figure.dpi200.png",
            "fig/genome_fraction_refs_figure.dpi200.png",
        ],
    conda: "conda/pandoc.yaml"
    shell:
        """
        pandoc --from markdown --to docx \
               --embed-resources --standalone --reference-doc {input.template} \
               --filter pandoc-crossref --csl {input.csl} --citeproc \
               --bibliography={input.bib} -s {input.source} -o {output}
        """


localrules:
    build_manuscript_docx,


# rule build_manuscript_pdf:
#     output:
#         "build/{stem}.pdf",
#     input:
#         source="doc/{stem}.md",
#         bib="doc/bibliography.bib",
#         csl="doc/citestyle.csl",
#         figures=lambda w: config["figures"][w.stem],
#     shell:
#         """
#         pandoc --from markdown --to pdf \
#                --filter pandoc-crossref --csl {input.csl} --citeproc \
#                --pdf-engine=xelatex \
#                --bibliography={input.bib} -s {input.source} -o {output}
#         """
#
#
# localrules:
#     build_manuscript_pdf,


rule compile_manuscript_submission:
    output: directory("build/submission")
    input:
        docx="build/manuscript.docx",
        coverletter="doc/manuscript/coverletter-genome-research.docx",
        fig1="fig/concept_diagram_figure.pdf",
        fig2="fig/benchmarking_figure.pdf",
        fig3="fig/hmp2_diversity_figure.pdf",
        fig4="fig/pangenomics_figure.pdf",
        fig5="fig/ucfmt_figure.pdf",
        figS1="fig/accuracy_by_depth_figure.pdf",
        figS2="fig/genome_fraction_refs_figure.pdf",
        tableS1="fig/hmp2_inferred_strains_supplementary_table1.tsv",
        tableS2="fig/ucfmt_focal_strain_genes_supplementary_table2.tsv",
    shell:
        """
        mkdir -p {output}
        cp {input.docx} {output}/manuscript.docx
        cp {input.coverletter} {output}/coverletter.docx
        cp {input.fig1} {output}/Figure_1.pdf
        cp {input.fig2} {output}/Figure_2.pdf
        cp {input.fig3} {output}/Figure_3.pdf
        cp {input.fig4} {output}/Figure_4.pdf
        cp {input.fig5} {output}/Figure_5.pdf
        cp {input.figS1} {output}/Supplementary_Figure_S1.pdf
        cp {input.figS2} {output}/Supplementary_Figure_S2.pdf
        cp {input.tableS1} {output}/Supplementary_Table_S1.tsv
        cp {input.tableS2} {output}/Supplementary_Table_S2.tsv
        """
