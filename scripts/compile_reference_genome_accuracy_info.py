#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys


if __name__ == "__main__":
    (
        strain_meta_path,
        strain_samples_path,
        species_free_samples_path,
        species_depth_path,
        strain_thresh_path,
        *reference_genome_accuracy_path_map_list,
    ) = sys.argv[1:]

    strain_meta = pd.read_table(strain_meta_path, index_col="strain").rename(
        columns=dict(
            num_sample="num_strain_samples",
            sum_depth="strain_depth_sum",
            max_depth="strain_depth_max",
        )
    )
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
    ).join(strain_meta, on="strain")

    result.to_csv(sys.stdout, sep="\t", index=False)
