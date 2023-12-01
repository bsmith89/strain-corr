#!/usr/bin/env python3
# "{input.script} {input.gene} {input.filt} {input.mgen} {input.preparation} {input.stool} {input.subject} {output}"

import pandas as pd
import sys
import scipy as sp
import scipy.stats
from tqdm import tqdm
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    strain_comm_inpath = sys.argv[1]
    strain_gene_inpath = sys.argv[2]
    strain_filt_inpath = sys.argv[3]
    mgen_inpath = sys.argv[4]
    preparation_inpath = sys.argv[5]
    stool_inpath = sys.argv[6]
    subject_inpath = sys.argv[7]
    frac_thresh = float(sys.argv[8])
    num_thresh = int(sys.argv[9])
    outpath = sys.argv[10]

    # Load data
    strain_frac = (
        pd.read_table(strain_comm_inpath, index_col=["sample", "strain"])
        .community.unstack("strain")
        .rename_axis(
            index="mgen_id",
            columns="genome_id",
        )
        .rename(columns=str)
    )
    strain_list = idxwhere(
        pd.read_table(strain_filt_inpath, index_col="genome_id")
        .rename(str)
        .astype(bool)
        .passes_filter
    )
    strain_gene = pd.read_table(strain_gene_inpath, index_col="gene_id")

    # Load metadata
    mgen = pd.read_table(mgen_inpath, index_col="library_id")
    preparation = pd.read_table(preparation_inpath, index_col="preparation_id")
    stool = pd.read_table(stool_inpath, index_col="stool_id")
    subject = pd.read_table(subject_inpath, index_col="subject_id")
    mgen_meta = (
        mgen.join(
            preparation,
            on="preparation_id",
            lsuffix="_mgen",
            rsuffix="_preparation",
        )
        .join(stool, on="stool_id")
        .join(subject, on="subject_id")
    )

    # A strain is present in a subject if at least one of their samples
    # has it at >30% fraction.
    u = (strain_frac > 0.3).groupby(mgen_meta.subject_id).any()[strain_list]

    # A strain is present in a subject if at least two of their samples
    # have it at >30% fraction.
    u = (strain_frac[strain_list] > frac_thresh).groupby(
        mgen_meta.subject_id
    ).sum() >= num_thresh

    # Only consider subjects and strains with at least one of the SPGC strains.
    subject_list = idxwhere(u.any(1))
    u = u.loc[subject_list]

    v = strain_gene[u.columns]
    subject_strain_gene_content = (u @ v.T).T > 0
    subject_strain_gene_prevalence = subject_strain_gene_content.mean(1)
    num_strain_subjects = subject_strain_gene_content.shape[1]

    subject_ibd_diagnosis = subject.loc[
        subject_strain_gene_content.columns
    ].ibd_diagnosis

    result = {}
    for gene in tqdm(subject_strain_gene_content.index):
        contingency = (
            pd.DataFrame(
                dict(
                    diagnosis=subject_ibd_diagnosis,
                    gene=subject_strain_gene_content.loc[gene].map(
                        {True: "gene_present", False: "gene_absent"}
                    ),
                )
            )
            .value_counts()
            .unstack("diagnosis")
            .reindex(
                columns=["CD", "UC", "nonIBD"],
                index=["gene_present", "gene_absent"],
            )
            .fillna(0)
            .assign(IBD=lambda x: x.CD + x.UC)
            .astype(int)
        )
        result[gene] = contingency.stack()

    result = pd.DataFrame(result).T
    result.columns = result.columns.to_flat_index().map(
        {
            ("gene_present", "CD"): "present-CD",
            ("gene_present", "UC"): "present-UC",
            ("gene_present", "nonIBD"): "present-nonIBD",
            ("gene_present", "IBD"): "present-IBD",
            ("gene_absent", "CD"): "absent-CD",
            ("gene_absent", "UC"): "absent-UC",
            ("gene_absent", "nonIBD"): "absent-nonIBD",
            ("gene_absent", "IBD"): "absent-IBD",
        }
    )
    result = result.rename_axis(index="gene_id").assign(
        fisher_exact_pvalue_ibd=lambda x: x.apply(
            lambda y: sp.stats.fisher_exact(
                [
                    [y["present-IBD"], y["present-nonIBD"]],
                    [y["absent-IBD"], y["absent-nonIBD"]],
                ]
            )[1],
            axis=1,
        ),
        oddsratio_pc_ibd=lambda x: (
            ((x["present-IBD"] + 1) / (x["present-nonIBD"] + 1))
            / ((x["absent-IBD"] + 1) / (x["absent-nonIBD"] + 1))
        ),
    )

    result.to_csv(outpath, sep="\t")
