#!/usr/bin/env python3

import pandas as pd
import sfacts as sf
import sys
import scipy as sp


if __name__ == "__main__":
    strain_geno_inpath = sys.argv[1]
    genome_id = sys.argv[2]
    spgc_geno_inpath = sys.argv[3]
    species_depth_inpath = sys.argv[4]
    strain_sample_inpath = sys.argv[5]
    mgen_list_inpath = sys.argv[6]
    outpath = sys.argv[7]

    strain_geno = sf.Metagenotype.load(strain_geno_inpath).to_estimated_genotype()
    spgc_agg_mgtp = sf.Metagenotype.load(spgc_geno_inpath)
    spgc_geno = spgc_agg_mgtp.to_estimated_genotype()
    with open(mgen_list_inpath) as f:
        # Parse a list of strings excluding the first line.
        mgen_list = [x.strip() for x in f][1:]

    sample_to_strain = (
        pd.read_table(strain_sample_inpath, index_col=["sample"])
        .squeeze()
        .reindex(mgen_list, fill_value=-1)
    )
    species_depth = (
        pd.read_table(
            species_depth_inpath,
            names=["sample", "depth"],
            index_col="sample",
        )
        .depth.squeeze()
        .reindex(mgen_list, fill_value=0)
    )

    strain_sample_depths = (
        species_depth.groupby(sample_to_strain)
        .agg(["sum", "max"])
        .rename(columns={"sum": "strain_depth_sum", "max": "strain_depth_max"})
    ).rename(int)

    strain_to_spgc_geno_diss = strain_geno.sel(strain=[genome_id]).cdist(spgc_geno)

    _strain_geno, _spgc_geno = strain_geno.data.to_pandas().align(
        spgc_geno.data.to_pandas(), axis=1, fill_value=0.5
    )
    strain_to_spgc_geno_num_positions = pd.DataFrame(
        sp.spatial.distance.cdist(
            _strain_geno != 0.5,
            _spgc_geno != 0.5,
            metric=lambda x, y: (x & y).sum(),
        ),
        index=_strain_geno.index,
        columns=_spgc_geno.index,
    ).astype(int)

    result = (
        strain_to_spgc_geno_diss.stack()
        .to_frame("genotype_dissimilarity")
        .rename_axis(index=["genome_id", "strain"])
        .join(strain_sample_depths)
        )
    result = (result
        .join(sample_to_strain.value_counts().to_frame("num_strain_samples").rename_axis("strain"))
        .sort_values("genotype_dissimilarity")
        .assign(
            num_geno_positions_compared=strain_to_spgc_geno_num_positions.stack(),
            species_depth_sum=species_depth.sum(),
            species_depth_max=species_depth.max(),
        )
    )

    result.to_csv(outpath, sep="\t")
