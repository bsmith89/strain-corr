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
config["species_group_to_spgc_stem"] = (
    pd.read_table("meta/species_group.tsv")
    .astype(str)
    .set_index(["species_id", "species_group_id"])
    .spgc_stem
)


config["genome"] = pd.read_table("meta/genome.tsv", dtype=str, index_col=["genome_id"])

config["genome_group"] = (
    pd.read_table("meta/genome_group.tsv")
    .groupby("genome_group_id")
    .genome_id.apply(lambda x: sorted(list(set(x))))
)

# NOTE: This function is used, e.g. in snake/reference_genome.smk and
# snake/benchmark_strain_reconstruction.smk to gather a list of reference
# genomes for each species.
def species_group_genomes(species, group):
    genome_group_list = config["genome_group"][group]
    d = config["genome"].loc[genome_group_list]
    strain_list = idxwhere(d.species_id == species)
    # print(strain_list)
    return strain_list


rule debug_species_group_genomes:
    output:
        "data/group/{group}/species/sp-{species}/genomes.list",
    params:
        test=lambda w: species_group_genomes(w.species, w.group),
    shell:
        "echo {params.test}; false"


config["species_to_panphlan"] = pd.read_table(
    "meta/species_to_panphlan.tsv", dtype=str, index_col=["species_id"]
).panphlan_species_id.rename(str)

config["species_to_spanda"] = pd.read_table(
    "meta/species_to_panphlan.tsv", dtype=str, index_col=["species_id"]
).spanda_species_id.rename(str)


rule process_hmp2_metadata:
    output:
        subject="meta/subject.tsv",
        visit="meta/visit.tsv",
        stool="meta/stool.tsv",
        preparation="meta/preparation.tsv",
        mgen="meta/mgen.tsv",
        mtab="meta/mtab.tsv",
    input:
        script="scripts/parse_hmp2_metadata_tables.py",
        raw="raw/hmp2_metadata_2018-08-20.csv",
    shell:
        """
        cat {input.raw} | {input.script} {output.subject} {output.visit} {output.stool} {output.preparation} {output.mgen} {output.mtab}
        """


rule download_uhgg_metadata:
    output:
        "ref/uhgg_genomes_all_v2.tsv",
    params:
        url="http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/genomes-all_metadata.tsv",
    shell:
        curl_recipe


config["figures"]["submission"] = []


config["een_mgen_local_src"] = pd.read_table(
    "meta/een-mgen/gtpro_local.tsv", index_col="library_id"
)

config["midasdb_uhgg_new_species_genome"] = (
    pd.read_table(
        "ref/midasdb_uhgg_new/genomes.tsv",
        dtype=str,
    )
    .groupby("species")
    .genome.apply(list)
)
