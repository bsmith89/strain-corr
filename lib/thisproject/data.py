import pandas as pd


def load_species_taxonomy(path):
    species_taxonomy = (
        pd.read_table(path, names=["genome_id", "species_id", "taxonomy_string"])
        .assign(species_id=lambda x: x.species_id.astype(str))
        .set_index("species_id")[["taxonomy_string"]]
        .assign(taxonomy_split=lambda x: x.taxonomy_string.str.split(";"))
    )
    for level_name, level_number in [
        ("p__", 1),
        ("c__", 2),
        ("o__", 3),
        ("f__", 4),
        ("g__", 5),
        ("s__", 6),
    ]:
        species_taxonomy = species_taxonomy.assign(
            **{
                level_name: species_taxonomy.taxonomy_split.apply(
                    lambda x: x[level_number]
                )
            }
        )
    species_taxonomy = species_taxonomy.drop(columns=["taxonomy_split"])
    return species_taxonomy


def load_species_depth(path):
    return (
        pd.read_table(
            path,
            index_col=["sample", "species_id"],
            dtype={"sample": str, "species_id": str, "depth": float},
        )
        .squeeze()
        .unstack("species_id", fill_value=0)
        .rename(str, axis="columns")
    )


def load_single_species_depth(path):
    return pd.read_table(
        path,
        names=["sample", "depth"],
        index_col="sample",
    ).squeeze()


def load_single_species_correlation(path):
    return pd.read_table(
        path, names=["sample", "correlation"], index_col="sample"
    ).squeeze()
