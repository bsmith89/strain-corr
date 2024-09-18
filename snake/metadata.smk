# {{{2 Data Configuration
from hashlib import md5

config["mgen"] = pd.read_table("meta/mgen_to_reads.tsv", index_col="mgen_id")
_mgen_group = pd.read_table("meta/mgen_group.tsv")

config["mgen_group"] = {}
for mgen_group_id, d in _mgen_group.groupby("mgen_group_id"):
    config["mgen_group"][mgen_group_id] = d.mgen_id.tolist()

config["species_group"] = (
    pd.read_table("meta/species_group.tsv")
    .astype(str)
    .groupby("species_group_id")
    .species_id.apply(lambda x: sorted(list(set(x))))
)


def hash_set_from_iterable(species_set):
    """A constant hash for any finite iterable with the same elements."""
    return md5(str(sorted(list(set(species_set)))).encode("utf-8")).hexdigest()


config["species_group_to_hash"] = {}
config["hash_to_species_set"] = {}
for species_group in config["species_group"].keys():
    species_set = sorted(list(set(config["species_group"][species_group])))
    _hash = str(hash_set_from_iterable(species_set))
    config["species_group_to_hash"][species_group] = _hash
    config["hash_to_species_set"][_hash] = species_set


config["species_group_to_sfacts_stem"] = (
    pd.read_table("meta/species_group.tsv")
    .astype(str)
    .set_index(["species_id", "species_group_id"])
    .sfacts_stem
)


def get_sfacts_stem(config, species, group):
    return config["species_group_to_sfacts_stem"][(species, group)]


config["species_group_to_spgc_stem"] = (
    pd.read_table("meta/species_group.tsv")
    .astype(str)
    .set_index(["species_id", "species_group_id"])
    .spgc_stem
)


def get_spgc_stem(config, species, group):
    return config["species_group_to_spgc_stem"][(species, group)]


rule download_uhgg_metadata:
    output:
        "ref/uhgg_genomes_all_v2.tsv",
    params:
        url="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv",
    shell:
        curl_recipe


# TODO (2024-06-10): Re-run GT-Pro for all v20 genomes?
if os.path.exists("ref/midasdb_uhgg_v20/genomes.tsv"):
    config["midasdb_uhgg_v20_species_genome"] = (
        pd.read_table(
            "ref/midasdb_uhgg_v20/genomes.tsv",
            dtype=str,
        )
        .groupby("species")
        .genome.apply(list)
    )
