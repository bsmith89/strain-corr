#!/usr/bin/env python3
# "{input.script} {input.pickle} {input.spgc_filt} {input.ref_filt} {params.thresh} {output}"
import pickle
import sys
import pandas as pd
from fastcluster import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from lib.pandas_util import idxwhere

if __name__ == "__main__":
    pickle_inpath = sys.argv[1]
    # spgc_filt_inpath = sys.argv[2]
    # ref_filt_inpath = sys.argv[3]
    thresh = float(sys.argv[2])
    outpath = sys.argv[3]

    # Import data
    with open(pickle_inpath, "rb") as f:
        d = pickle.load(f)
    # spgc_list = idxwhere(
    #     pd.read_table(spgc_filt_inpath, index_col="genome_id")
    #     .passes_filter.rename(str)
    #     .astype(bool)
    # )
    # ref_list = idxwhere(
    #     pd.read_table(ref_filt_inpath, index_col="genome_id").passes_filter.astype(bool)
    # )

    # Filter distance matrix to spgc and refs passing their respective filters.
    original_cdmat = d["cdmat"]
    original_labels = d["labels"]
    # labels = spgc_list + ref_list
    # assert (set(labels) - set(original_labels)) == set([])

    cdmat = squareform(
        pd.DataFrame(
            squareform(original_cdmat), index=original_labels, columns=original_labels
        )
        #.loc[labels, labels]
    )

    linkage = linkage(cdmat, method="average")
    cluster = fcluster(linkage, t=thresh, criterion="distance")

    # cluster = pd.Series(cluster, index=labels)
    cluster = pd.Series(cluster, index=original_labels)
    cluster.to_csv(outpath, sep="\t", header=False)
