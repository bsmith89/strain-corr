#!/usr/bin/env python3

import pandas as pd
import gzip
from io import StringIO
import sys
from lib.util import info
import numpy as np
import sfacts as sf

if __name__ == "__main__":
    snpdict_path = sys.argv[1]
    vcf_path = sys.argv[2]
    species_id = int(sys.argv[3])
    outpath = sys.argv[4]

    info("Scanning VCF for data start")
    # Count header rows prefixed by '##'.
    with gzip.open(vcf_path, mode="rt") as f:  # Assume gzipped input.
        for skiprows, line in enumerate(f):
            if line.startswith("##"):
                continue
            else:
                break
    info(f"Reading {vcf_path}")
    vcf = (
        pd.read_table(
            vcf_path,
            skiprows=skiprows,
            dtype={
                "#CHROM": str,
                "POS": int,
                "ID": int,
                "REF": str,
                "ALT": str,
                "QUAL": str,
                "FILTER": str,
                "INFO": str,
            },
            sep="\t",
        )
        .rename(columns={"#CHROM": "CHROM"})
        .drop(columns=["CHROM", "POS", "QUAL", "FILTER", "INFO", "FORMAT"])
        .rename(columns={"ID": "position"})
        .assign(species_id=species_id)
        .set_index(["species_id", "position"])
    )
    info(f"Reading {snpdict_path}")
    snpdict = pd.read_table(
        snpdict_path,
        names=["species_id", "position", "contig_id", "contig_position", "ref", "alt"],
        index_col=["species_id", "position"],
    )

    assert (
        vcf.xs(species_id).loc[snpdict.xs(species_id).index].ALT.apply(len) == 1
    ).all()

    info(f"Reshaping data")
    data = (
        (
            vcf.xs(species_id)
            .loc[snpdict.xs(species_id).index]
            .drop(columns=["REF", "ALT"])
            .rename_axis(columns="sample")
            .stack()
            == "1:0:0:0"
        )
        .astype(int)
        .to_frame(name="ref")
        .assign(alt=lambda x: 1 - x.ref)
        .rename_axis(columns="allele")
        .assign(species_id=species_id)
        .set_index("species_id", append=True)
        .stack()
        .to_xarray()
        .sel(species_id=species_id)
        .astype(np.ushort)
        .transpose("sample", "position", "allele")
    )

    data = sf.Metagenotype(data)

    info(f"Writing to {outpath}")
    data.dump(outpath)
    info("Finished")
