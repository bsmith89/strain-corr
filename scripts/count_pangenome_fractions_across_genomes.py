#!/usr/bin/env python3
# "{input.script} {input.prevalence} {input.gene} {params.core_vs_shell} {params.shell_vs_cloud} {output}"

import sys
import pandas as pd
import numpy as np
from warnings import warn


def classify_gene_prevalence(x, core_vs_shell, shell_vs_cloud):
    assert core_vs_shell > shell_vs_cloud
    if core_vs_shell == shell_vs_cloud:
        warn("Core/shell and shell/cloud thresholds are equal.")
    tag = np.where(
        x >= core_vs_shell, "core", np.where(x >= shell_vs_cloud, "shell", "cloud")
    )
    return pd.Series(tag, index=x.index)


if __name__ == "__main__":
    gene_presence_inpath = sys.argv[2]
    prevalence_inpath = sys.argv[1]
    core_vs_shell_thresh = float(sys.argv[3])
    shell_vs_cloud_thresh = float(sys.argv[4])
    outpath = sys.argv[5]

    prevalence = pd.read_table(
        prevalence_inpath,
        names=["gene_id", "prevalence"],
        index_col="gene_id",
    ).prevalence
    prevalence_class = classify_gene_prevalence(
        prevalence,
        core_vs_shell=core_vs_shell_thresh,
        shell_vs_cloud=shell_vs_cloud_thresh,
    )
    prevalence_class_dummy = prevalence_class.to_frame("label").assign(
        core=lambda x: x.label == "core",
        shell=lambda x: x.label == "shell",
        cloud=lambda x: x.label == "cloud",
    )[["core", "shell", "cloud"]]
    gene = (
        pd.read_table(gene_presence_inpath, index_col="gene_id")
        .rename_axis(columns="strain")
        .loc[prevalence_class_dummy.index]
    )
    result = gene.T @ prevalence_class_dummy
    result.to_csv(outpath, sep="\t")
