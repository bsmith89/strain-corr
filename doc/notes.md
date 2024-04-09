# Analysis code and metadata for: "Accurate estimation of intraspecific microbial gene content variation in metagenomic data with MIDAS v3 and StrainPGC"

Our analysis workflow and downstream analyses are implemented as a Snakemake pipeline:

![The Workflow](doc/rulegraph.svg)

Code to reproduce our main text analyses and figures (as well as supplementary
results) are then integrated into analysis notebooks:

- Supplementary Results 1: Performance benchmarking in a synthetic community
    - [`nb/analyze_benchmarking_results.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_benchmarking_results.ipynb)
- Supplementary Results 2: Distribution of strains in HMP2
    - [`nb/analyze_distribution_of_hmp2_strains.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_distribution_of_hmp2_strains.ipynb)
- Supplementary Results 3: Diversity of inferred strains
    - [`nb/analyze_hmp2_strain_diversity.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_hmp2_strain_diversity.ipynb)
- Supplementary Results 4: Gene prevalence in inferred strains
    - [`nb/analyze_pangenome_fractions_hmp2_strains.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_pangenome_fractions_hmp2_strains.ipynb)
- Supplementary Results 5: Co-occurrence clustering
    - [`nb/analyze_gene_clusters_in_hmp2_strains.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_gene_clusters_in_hmp2_strains.ipynb)
- Supplementary Results 6: Tracking and comparison of E. coli strain gene content in UCFMT
    - [`nb/analyze_ucfmt_donor_strains_102506.ipynb`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/nb/analyze_ucfmt_donor_strains_102506.ipynb)

Supplementary tables:

- Supplementary Table 1: Details about all inferred strains in HMP2
    - [`hmp2_inferred_strains_supplementary_table1.tsv`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/fig/hmp2_inferred_strains_supplementary_table1.tsv)
- Supplementary Table 2: Details about gene content of E. coli strain-6 vs. strain-9 in UCFMT
    - [`ucfmt_focal_strain_genes_supplementary_table2.tsv`](https://github.com/bsmith89/StrainPGC-manuscript/blob/with-results/fig/ucfmt_focal_strain_genes_supplementary_table2.tsv)
