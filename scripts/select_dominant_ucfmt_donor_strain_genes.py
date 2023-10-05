#!/usr/bin/env python3

import pandas as pd
import sys
from lib.pandas_util import idxwhere


MIN_NUM_DONOR_SAMPLES = 1
MIN_NUM_SUBJECTS = 2


if __name__ == "__main__":
    frac_path = sys.argv[1]
    spgc_filt_path = sys.argv[2]
    detection_thresh = float(sys.argv[3])
    mgen_path = sys.argv[4]
    sample_path = sys.argv[5]
    subject_path = sys.argv[6]
    gene_hit_path = sys.argv[7]
    outpath = sys.argv[8]

    # Load/compile sample metadata
    mgen = pd.read_table(mgen_path, index_col="mgen_id")
    sample = pd.read_table(sample_path, index_col="sample_id")
    subject = pd.read_table(subject_path, index_col="subject_id")
    mgen_meta = mgen.join(sample, on="sample_id").join(subject, on="subject_id")

    comm = pd.read_table(
        frac_path,
        dtype={"sample": str, "strain": str, "community": float},
        index_col=["sample", "strain"],
    ).community.unstack()
    spgc_filt = (
        pd.read_table(spgc_filt_path, index_col="genome_id")
        .rename(str)
        .rename_axis(index="strain")
        .passes_filter
    )
    gene_hit = pd.read_table(gene_hit_path, index_col="gene_id")

    donor_subject_followup_tally = (
        comm.rename(index=dict(sample="mgen_id"))
        # Drop all but recipient follow-ups.
        .drop(
            idxwhere((~mgen_meta.recipient) | (mgen_meta.sample_type != "followup")),
            errors="ignore",
        )
        .gt(detection_thresh)
        .groupby(mgen_meta.subject_id)
        # .join(mgen_meta, on="mgen_id")
        .any()
        .groupby(subject.donor_subject_id)
        .sum()
        .stack()
        .sort_values(ascending=False)
    )

    donor_sample_tally = (
        comm.rename(index=dict(sample="mgen_id"))
        # Drop all but donor samples.
        .drop(
            idxwhere(mgen_meta.sample_type != "donor"),
            errors="ignore",
        )
        .gt(detection_thresh)
        .groupby(mgen_meta.subject_id)
        .sum()
        .rename_axis("donor_subject_id")
        .stack()
        .sort_values(ascending=False)
    )

    strain_selection = (
        pd.DataFrame(
            dict(
                num_subject=donor_subject_followup_tally,
                num_donor_sample=donor_sample_tally,
            )
        )
        .fillna(0)
        .sort_values("num_subject", ascending=False)
        .join(spgc_filt, on="strain")
        .fillna({"passes_filter": False})
    )

    donor_strain = dict(
        strain_selection
        # Require >1 recipients and >1 donor samples.
        # NOTE: This excludes all donors except D0097 and D0044.
        [
            lambda x: (x.num_donor_sample >= MIN_NUM_DONOR_SAMPLES)
            & (x.num_subject > MIN_NUM_SUBJECTS)
            & (x.passes_filter)
        ]
        .groupby("donor_subject_id")
        .head(1)
        .index.to_frame()[["strain", "donor_subject_id"]]
        .values
    )

    result = gene_hit.reindex(columns=donor_strain.keys()).rename(columns=donor_strain)[
        lambda x: x.sum(1) > 0  # Drop genes with no hits.
    ]
    result.to_csv(outpath, sep="\t")
