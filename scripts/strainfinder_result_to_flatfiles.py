#!/usr/bin/env python2

from __future__ import print_function
import pandas as pd
import numpy as np
import cPickle as cp
import sys
from sfinder import EM, Estimate, Data


if __name__ == "__main__":

    with open(sys.argv[1], "r") as f:
        em = cp.load(f)

    with open(sys.argv[2], "r") as f:
        lines = f.readlines()
        sample_list = lines[0].strip().split(",")
        position_list = lines[1].strip().split(",")
        input_allele_list = lines[2].strip().split(",")

    allele_list = ["A", "C", "G", "T"]

    pi = em.estimates[-1].z
    gamma_ = em.estimates[-1].p
    gamma = gamma_[:, :, : len(input_allele_list)]

    nstrain = pi.shape[1]
    strain_list = list(range(nstrain))

    # import pdb; pdb.set_trace()
    pi_df = pd.DataFrame(pi, index=sample_list, columns=strain_list).rename_axis(
        index="sample", columns="strain"
    ).stack().rename('community')
    gamma_df = pd.Series(
        gamma.flatten(),
        index=pd.MultiIndex.from_product(
            [strain_list, position_list, input_allele_list],
            names=["strain", "position", "allele"],
        ),
        name="genotype",
    )

    pi_df.to_csv(sys.argv[3], sep='\t', header=True)
    gamma_df.to_csv(sys.argv[4], sep='\t', header=True)
