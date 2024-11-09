# Analysis code and metadata for: "Accurate estimation of intraspecific microbial gene content variation in metagenomic data with MIDAS v3 and StrainPGC"

For the tools---[StrainPGC](https://github.com/bsmith89/StrainPGC) and
[MIDAS](https://github.com/czbiohub-sf/MIDAS)---see their respective GitHub repositories.

The StrainPGC repository includes a user-friendly version of our integrated workflow.

Our complete analysis workflow and the downstream analyses for this paper are implemented as a Snakemake pipeline:

![The Workflow](doc/static/filegraph.svg)

Code to reproduce our main text analyses and figures (as well as supplementary
results) are then integrated into analysis notebooks:

- Performance benchmarking in a synthetic community
    - [`nb/analyze_benchmarking_results.ipynb`](nb/analyze_benchmarking_results.ipynb)
- Performance benchmarking with E. coli spike-in genomes
    - [`nb/analyze_spikein_benchmark.ipynb`](nb/analyze_spikein_benchmark.ipynb)
- Distribution of strains in HMP2
    - [`nb/analyze_distribution_of_hmp2_strains.ipynb`](nb/analyze_distribution_of_hmp2_strains.ipynb)
- Diversity of inferred strains
    - [`nb/analyze_hmp2_strain_diversity.ipynb`](nb/analyze_hmp2_strain_diversity.ipynb)
- Gene prevalence in inferred strains
    - [`nb/analyze_pangenome_fractions_hmp2_strains.ipynb`](nb/analyze_pangenome_fractions_hmp2_strains.ipynb)
- Co-occurrence clustering
    - [`nb/analyze_gene_clusters_in_hmp2_strains.ipynb`](nb/analyze_gene_clusters_in_hmp2_strains.ipynb)
- Tracking and comparison of E. coli strain gene content in UCFMT
    - [`nb/analyze_ucfmt_donor_strains_102506.ipynb`](nb/analyze_ucfmt_donor_strains_102506.ipynb)
