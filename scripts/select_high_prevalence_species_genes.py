#!/usr/bin/env python3

# USAGE: {input.script} {input.copy_number} {params.trim_quantile} {params.threshold} {output}

import xarray as xr
import sys
from lib.util import info
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    reference_copy_number_inpath = sys.argv[1]
    trim_quantile = float(sys.argv[2])
    threshold = float(sys.argv[3])
    outpath = sys.argv[4]

    info("Loading reference table.")
    reference_copy_number = xr.load_dataarray(reference_copy_number_inpath)

    # Representative genome selection (remove outliers based on gene number).
    gene_count = reference_copy_number.sum("gene_id")
    low_q, high_q = (
        gene_count.quantile([trim_quantile, 1 - trim_quantile]).astype(int).values
    )
    inter_quantile_range_genomes = idxwhere(
        ((gene_count > low_q) & (gene_count < high_q)).to_series()
    )
    num_inter_quantile_range_genomes = len(inter_quantile_range_genomes)
    info(
        f"Selected {num_inter_quantile_range_genomes} genomes "
        f"with between {low_q} and {high_q} genes "
        f"(respectively, the {trim_quantile} inter-quantile range)."
    )

    # Calculate frequency as a single-copy gene.
    frac_single_copy = (
        reference_copy_number.sel(genome_id=inter_quantile_range_genomes) == 1
    ).mean("genome_id")

    # Select single-copy genes.
    single_copy_genes = idxwhere((frac_single_copy > threshold).to_series())
    num_single_copy_genes = len(single_copy_genes)
    info(
        f"Found {num_single_copy_genes} single copy genes "
        f"in representative genomes."
    )

    # Write output.
    with open(outpath, "w") as f:
        for gene in single_copy_genes:
            print(gene, file=f)
