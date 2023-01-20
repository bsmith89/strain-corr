#!/usr/bin/env python2
"""Convert TSV of metagenotypes to a format StrainFinder understands as an alignment.

Usage:

    python2 metagenotype_to_strainfinder_alignment.py <INPUT> <OUTPUT>


<INPUT> must be a TSV of four columns (and a matching header):

    sample | position | allele | metagenotype

"""

from __future__ import print_function
import pandas as pd
import cPickle as cp
import sys
from itertools import product

if __name__ == "__main__":
    print("Loading", file=sys.stderr)
    in_data = (
        pd.read_table(sys.argv[1], index_col=["sample", "position", "allele"])
        .squeeze()
        .astype(int)
        .sort_index()  # NOTE: This is a stop-gap. See below.
    )

    print("Parsing", file=sys.stderr)
    # TODO: These three *_list variables may be out of order
    # but apparently they work if the index was previously sorted
    # (see above). Ultimately, we want each of the lists sorted
    # like the indices of each axis of the multi-dimensional array.
    # This approach (pre-sorting the index) works (and we were
    # already doing that below), but I'm not 100% confident that Index.unique()
    # always results in a correctly sorted list.
    sample_list = in_data.index.unique("sample")
    position_list = in_data.index.unique("position")
    input_allele_list = in_data.index.unique("allele")
    n = len(sample_list)
    g = len(position_list)
    input_a = len(input_allele_list)
    if input_a < 4:
        allele_list = ["A", "C", "G", "T"]
    elif input_a == 4:
        allele_list = input_allele_list
    else:
        raise ValueError("Input TSV must not have more than 4 distinct alleles.")

    print("Reshaping", file=sys.stderr)
    full_size_data = (
        in_data.rename(
            {x: y for x, y in zip(input_allele_list, allele_list)}, level="allele"
        )
        .reindex(product(sample_list, position_list, allele_list))
        .sort_index()
        .fillna(0)
    )
    out_data = full_size_data.values.reshape((n, g, 4))

    print("Writing", file=sys.stderr)
    with open(sys.argv[2], "wb") as f:
        cp.dump(out_data, file=f)

    with open(sys.argv[3], "w") as f:
        print(','.join(str(x) for x in sample_list), file=f)
        print(','.join(str(x) for x in position_list), file=f)
        print(','.join(str(x) for x in input_allele_list), file=f)

    print("Checking", file=sys.stderr)
    out_index = full_size_data.index.to_numpy().reshape((g, n, 4))
    reconstruction = (
        pd.Series(
            out_data.flatten(),
            index=pd.MultiIndex.from_tuples(
                out_index.flatten(), names=["sample", "position", "allele"]
            ),
        )
        .unstack("allele")
        .rename(columns={"A": "alt", "C": "ref"})
        [['alt', 'ref']]
        .stack()
    )
    assert (in_data.sort_index() == reconstruction.sort_index()).all()

    print("Done", file=sys.stderr)
