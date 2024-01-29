#!/usr/bin/env python3

import sfacts as sf
import sys
import xarray as xr


if __name__ == "__main__":
    spgc_inpath = sys.argv[1]
    agg_mgtp_inpath = sys.argv[2]
    ambig_thresh = float(sys.argv[3])
    outpath = sys.argv[4]

    # Load data
    spgc = xr.load_dataset(spgc_inpath)

    # Strain genotype entropy
    strain_total_mgtp = sf.Metagenotype.load(agg_mgtp_inpath)
    strain_metagenotype_entropy = strain_total_mgtp.entropy()
    num_confident_geno_positions = (
        (strain_total_mgtp.dominant_allele_fraction() > (1 - ambig_thresh))
        .sum("position")
        .rename(sample="strain")
    )

    strain_sample_meta = (
        spgc[
            [
                "num_gene",
                "num_strain_sample",
                "sum_strain_depth",
                "max_strain_depth",
                "species_gene_frac",
                "log_selected_gene_depth_ratio_std",
            ]
        ]
        .to_dataframe()
        .assign(
            num_geno_positions=num_confident_geno_positions,
            strain_metagenotype_entropy=strain_metagenotype_entropy,
        )
    )

    # Write output
    strain_sample_meta.rename_axis(index="genome_id").to_csv(outpath, sep="\t")
