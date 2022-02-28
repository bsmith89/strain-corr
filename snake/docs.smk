rule sort_bib_from_raw:
    output: 'doc/bibliography.bib'
    input:
        script='scripts/sort_bib.py',
        bib=['doc/bibliography_raw.bib'],
    shell:
        "{input.script} {input.bib} > {output}"


rule render_figure_to_png:
    output:
        "fig/{stem}_figure.w{width}.png",
    input:
        "doc/static/{stem}_figure.svg",
    params:
        width=lambda w: int(w.width),
    shell:
        """
        inkscape {input} --export-width={params.width} --export-filename {output}
        """

rule render_pdf_to_png_imagemagick:
    output:
        "fig/{stem}.dpi{dpi}.png",
    input:
        "fig/{stem}.pdf",
    params:
        dpi=lambda w: int(w.dpi),
    shell:
        """
        convert -units PixelsPerInch -density {params.dpi} {input} {output}
        """

rule render_pdf_to_tiff_imagemagick:
    output:
        "fig/{stem}.dpi{dpi}.tiff",
    input:
        "fig/{stem}.pdf",
    params:
        dpi=lambda w: int(w.dpi),
    shell:
        """
        convert -units PixelsPerInch -density {params.dpi} {input} {output}
        """

rule link_static_pdf_figure:
    output: 'fig/{stem}_figure.pdf'
    input: 'doc/static/{stem}_figure.pdf'
    shell:
        alias_recipe

rule pdf_to_eps:
    output: 'fig/{stem}.eps'
    input: 'fig/{stem}.pdf'
    shell:
        """
        cd fig
        pdf2ps {wildcards.stem}.pdf
        ps2eps {wildcards.stem}.ps
        rm {wildcards.stem}.ps
        """

rule render_figure_to_pdf:
    output:
        "fig/{stem}_figure.pdf",
    input:
        "doc/static/{stem}_figure.svg",
    shell:
        """
        inkscape {input} --export-filename {output}
        """

rule build_submission_docx:
    output:
        "build/{stem}.docx",
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        template="doc/static/example_style.docx",
        csl="doc/citestyle.csl",
        figures=lambda w: config['figures'][w.stem],
    shell:
        """
        pandoc --from markdown --to docx \
               --standalone --self-contained --reference-doc {input.template} \
               --filter pandoc-crossref --csl {input.csl} --citeproc \
               --bibliography={input.bib} -s {input.source} -o {output}
        """


localrules:
    build_submission_docx,


rule build_submission_pdf:
    output:
        "build/{stem}.pdf",
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        csl="doc/citestyle.csl",
        figures=lambda w: config['figures'][w.stem],
    shell:
        """
        pandoc --from markdown --to pdf \
               --filter pandoc-crossref --csl {input.csl} --citeproc \
               --pdf-engine=xelatex \
               --bibliography={input.bib} -s {input.source} -o {output}
        """


localrules:
    build_submission_pdf,
