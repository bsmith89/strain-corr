#!/usr/bin/env python3

import pandas as pd
import sfacts as sf
import xarray as xr
import sys
import scipy as sp
from lib.pandas_util import idxwhere


if __name__ == "__main__":
    # NOTE (2023-09-03): I call it "spgc" everywhere, but this is NOT the StrainPGC pipeline.
    # "{input.script} {input.spgc_geno} {input.ref_geno} {params.min_geno_diss} {input.ref_gene} {input.ref_meta} {params.min_completeness} {params.max_contamination} {output}"
    spgc_geno_inpath = sys.argv[1]
    ref_geno_inpath = sys.argv[2]
    min_geno_diss = float(sys.argv[3])
    ref_gene_inpath = sys.argv[4]
    ref_meta_inpath = sys.argv[5]
    species_id = sys.argv[6]
    min_ref_completeness = float(sys.argv[7])
    max_ref_contamination = float(sys.argv[8])
    outpath = sys.argv[9]

    ref_meta = (
        pd.read_table(ref_meta_inpath, index_col="Genome")
        .rename_axis(index="genome_id")[
            lambda x: x.MGnify_accession == "MGYG-HGUT-" + species_id[1:]
        ]
        .rename(lambda s: "UHGG" + s[len("GUT_GENOME") :])
        .assign(
            Completeness=lambda x: x.Completeness / 100,
            Contamination=lambda x: x.Contamination / 100,
        )
    )
    ref_list = idxwhere(
        (ref_meta.Completeness > min_ref_completeness)
        & (ref_meta.Contamination < max_ref_contamination)
    )

    ref_gene_content = xr.load_dataarray(ref_gene_inpath) > 0

    spgc_agg_mgtp = sf.Metagenotype.load(spgc_geno_inpath)

    # NOTE (2023-09-03): This renaming is only necessary because I fudged the
    # MIDASDB labels.
    ref_geno = (
        sf.Metagenotype.load(ref_geno_inpath)
        .rename_coords(sample=lambda s: "UHGG" + s[len("GUT_GENOME") :])
        .sel(sample=ref_list)
    )
    ref_has_counts = ref_geno.total_counts().to_pandas().isna()
    ref_geno = ref_geno.to_estimated_genotype()
    assert (
        False
    ), "TODO: Fix this script. It's current broken, and assigns the same reference to every genotype."

    spgc_has_counts = (
        (spgc_agg_mgtp.total_counts() > 0).to_pandas().rename_axis(index="strain")
    )
    spgc_has_counts, ref_has_counts = spgc_has_counts.align(
        ref_has_counts, axis=1, fill_value=False
    )
    num_positions_compared = pd.DataFrame(
        sp.spatial.distance.cdist(
            spgc_has_counts,
            ref_has_counts,
            metric=lambda x, y: (x & y).sum(),
        ),
        columns=ref_has_counts.index,
        index=spgc_has_counts.index,
    ).astype(int)

    position_list = list(
        set(spgc_agg_mgtp.position.values) & set(ref_geno.position.values)
    )

    # Filter references by completeness and contamination.

    # FIXME (2023-09-03): Why does spgc_to_ref_cdist * num_positions_compared
    # not an integer when genotypes are discretized??
    spgc_geno = spgc_agg_mgtp.sel(position=position_list).to_estimated_genotype()
    ref_geno = ref_geno.sel(position=position_list)
    spgc_to_ref_cdist = spgc_geno.cdist(ref_geno)

    spgc_gene_content = {}
    for strain in spgc_geno.strain.values:
        strain_match = pd.DataFrame(
            dict(
                ref_diss=spgc_to_ref_cdist.loc[strain],
                npos=num_positions_compared.loc[strain],
            )
        ).assign(
            matched_position=lambda x: ((1 - x.ref_diss) * x.npos),
            ref_diss_pc=lambda x: x.ref_diss + (1 / x.npos),
        )
        select_genome = strain_match[
            lambda x: x.ref_diss >= min_geno_diss
        ].ref_diss_pc.idxmin()
        spgc_gene_content[strain] = ref_gene_content.sel(
            genome_id=select_genome
        ).to_series()
    spgc_gene_content = (
        pd.DataFrame(spgc_gene_content).rename_axis(columns="strain").astype(int)
    )[lambda x: x.sum(1) > 0]
    spgc_gene_content.to_csv(outpath, sep="\t")
