#!/usr/bin/env python3

import pandas as pd
import sys
import subprocess


def depth_estimator(
    mean_coverage,
    fraction_covered,
    gene_length,
    read_length,
    fragment_threshold,
    maximum_correction_factor=1e6,
):
    # The edge_length is the total length of the "shadowed" zone on either end
    # of a gene, due to dropping reads where less than fragment_threshold
    # of their length is aligned to the gene.
    edge_length = read_length * fragment_threshold * 2
    # The edge_ratio is the ratio of this shadowed region to the total
    # gene length. This will determine a correction factor to be applied to the
    # naive depth estimate.
    edge_ratio = edge_length / gene_length
    maximum_edge_ratio = 1 - 1 / maximum_correction_factor
    edge_ratio_clipped = edge_ratio.clip(lower=0, upper=maximum_edge_ratio)
    # Notice that we're clipping the edge ratio complement,
    # since this number can be close to or below zero
    # when the gene is shorter than the "shadowed" zone.
    # Unfortunately, this would mean that we will sometimes estimate an
    # infinite or negative depth due to the correction
    # (and especially when the read length is much shorter
    # than the nominal_read_length).
    # There is not an obvious way to deal with this, and it is
    # an unfortunately consequence of the "--aln_cov" parameter passed
    # to MIDAS2 run_genes.
    depth = mean_coverage * fraction_covered / (1 - edge_ratio_clipped)
    # Note that if fragment_threshold == 0 or read_length == 0, then the depth estimator
    # reduces to mean_coverage * fraction_covered, the naive estimate.
    return depth


if __name__ == "__main__":
    minimum_alignment_coverage = float(sys.argv[1])
    nominal_read_length = int(sys.argv[2])
    maximum_correction_factor = float(sys.argv[3])
    outpath = sys.argv[4]

    data = []
    for arg in sys.argv[5:]:
        sample, path = arg.split("=")
        with subprocess.Popen(["lz4", "-dc", path], stdout=subprocess.PIPE) as proc:
            data.append(
                pd.read_table(proc.stdout).assign(
                    sample=sample,
                    depth=lambda x: depth_estimator(
                        x.mean_coverage,
                        x.fraction_covered,
                        x.gene_length,
                        nominal_read_length,
                        minimum_alignment_coverage,
                        maximum_correction_factor,
                    ),
                )[["gene_id", "sample", "depth"]]
            )
    data = (
        pd.concat(data).set_index(["gene_id", "sample"]).squeeze().to_xarray().fillna(0)
    )
    data.to_netcdf(outpath)
