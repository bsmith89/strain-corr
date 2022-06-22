# {{{2 Data Configuration

config["species_group"] = (
    pd.read_table("meta/species_group.tsv", dtype=str)
    .groupby("species_group_id")
    .apply(lambda d: d.species_id.to_list())
)

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
