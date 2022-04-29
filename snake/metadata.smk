# {{{2 Data Configuration

config["mgen"] = (
    pd.read_table("meta/mgen.tsv")
    .rename(columns={"library_id": "mgen_id", "external_id": "stem"})
    .set_index("mgen_id")
    .assign(
        r1=lambda x: f"hmp2/" + x.stem + "_R1.fastq.gz",
        r2=lambda x: f"hmp2/" + x.stem + "_R2.fastq.gz",
    )[["r1", "r2"]]
)

config["mgen_group"] = (
    pd.read_table("meta/mgen_group.tsv")
    .groupby("mgen_group_id")
    .apply(lambda d: d.mgen_id.to_list())
)

config["species_list"] = pd.read_table(
    "meta/species_group.tsv", names=["species_id", "species_group_id"]
).species_id.unique()

# _mgen_meta = "meta/mgen.tsv"
# if path.exists(_mgen_meta):
#     # FIXME: library_id should be mgen_id
#     _mgen = pd.read_table(_mgen_meta, index_col="library_id")
#     config["mgen"] = {}
#     for mgen_id, row in _mgen.iterrows():
#         config["mgen"][mgen_id] = {}
#         # FIXME: filename_r* should be in the table
#         # config["mgen"][mgen_id]["r1"] = row["filename_r1"]
#         # config["mgen"][mgen_id]["r2"] = row["filename_r2"]
# else:
#     warn(
#         dd(
#             f"""
#             Could not load config from `{_mgen_meta}`.
#             Check that path is defined and file exists.
#             """
#         )
#     )
#     config["mgen"] = {}
#
# _mgen_x_mgen_group_meta = "meta/mgen_x_mgen_group.tsv"
# if path.exists(_mgen_x_mgen_group_meta):
#     _mgen_x_mgen_group = pd.read_table(_mgen_x_mgen_group_meta)
#     config["mgen_group"] = {}
#     for mgen_group, d in _mgen_x_mgen_group.groupby("mgen_group"):
#         config["mgen_group"][mgen_group] = d.mgen_id.tolist()
# else:
#     warn(
#         dd(
#             f"""
#             Could not load config from `{_mgen_x_mgen_group_meta}`.
#             Check that path is defined and file exists.
#             """
#         )
#     )
#     config["mgen_group"] = {}

config["figures"]["submission"] = []


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
