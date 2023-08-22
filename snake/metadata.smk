# {{{2 Data Configuration

config["mgen"] = pd.read_table("meta/mgen_to_reads.tsv", index_col="mgen_id")
_mgen_group = pd.read_table("meta/mgen_group.tsv")

config["mgen_group"] = {}
for mgen_group, d in _mgen_group.groupby("mgen_group"):
    config["mgen_group"][mgen_group] = d.mgen_id.tolist()

config["species_group"] = (
    pd.read_table("meta/species_group.tsv")
    .astype(str)
    .groupby("species_group_id")
    .species_id.apply(list)
)

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


config["genome"] = pd.read_table(
    "meta/genome.tsv", dtype=str, index_col=["genome_id"]
).dropna(subset=["genome_path"])


# NOTE: This function is used, e.g. in snake/reference_genome.smk and
# snake/benchmark_strain_reconstruction.smk to gather a list of reference
# genomes for each species.
def species_genomes(species):
    strain_list = idxwhere(config["genome"].species_id == species)
    # assert len(strain_list) > 0
    return strain_list


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
