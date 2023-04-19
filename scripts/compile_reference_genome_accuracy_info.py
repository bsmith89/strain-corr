#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    (
        strain_samples_path,
        species_free_samples_path,
        species_depth_path,
        strain_thresh_path,
        *reference_genome_accuracy_path_map_list,
    ) = sys.argv[1:]

    sample_to_strain = pd.read_table(strain_samples_path, index_col="sample").strain
    num_strain_samples = sample_to_strain.value_counts()
    try:
        dominant_strain = num_strain_samples.idxmax()
    except ValueError:
        dominant_strain = -1

    species_free_samples = (
        pd.read_table(species_free_samples_path, names=["sample"]).squeeze().values
    )
    num_species_free_samples = len(species_free_samples)

    species_depth = pd.read_table(
        species_depth_path, names=["sample", "depth"], index_col="sample"
    ).depth
    strain_depth_std = species_depth.groupby(sample_to_strain).apply(
        lambda x: np.concatenate(
            [x.values, species_depth.loc[species_free_samples].values]
        ).std()
    )
    strain_depth_sum = species_depth.groupby(sample_to_strain).sum()
    strain_depth_max = species_depth.groupby(sample_to_strain).max()
    strain_depth_stats = pd.DataFrame(
        dict(
            num_strain_samples=num_strain_samples,
            strain_depth_sum=strain_depth_sum,
            strain_depth_max=strain_depth_max,
            strain_depth_std=strain_depth_std,
        )
    ).rename(int)

    num_reference_genomes = len(reference_genome_accuracy_path_map_list)
    num_xjin_strains = len(num_strain_samples.index)

    strain_thresh = pd.read_table(strain_thresh_path, index_col="strain")

    reference_genome_accuracy = []
    for key_value_pair in reference_genome_accuracy_path_map_list:
        genome, path = key_value_pair.split("=")
        data = pd.read_table(path, index_col="strain").assign(
            genome_id=genome,
            correlation_thresh=strain_thresh["correlation"],
            depth_thresh=strain_thresh["depth_low"],
        )
        reference_genome_accuracy.append(data)
    reference_genome_accuracy = pd.concat(reference_genome_accuracy).reset_index()

    result = reference_genome_accuracy.assign(
        dominant_strain=dominant_strain,
        total_num_reference_genomes=num_reference_genomes,
        total_num_xjin_strains=num_xjin_strains,
        num_species_free_samples=num_species_free_samples,
    ).join(strain_depth_stats, on="strain")

    result.to_csv(sys.stdout, sep="\t", index=False)
