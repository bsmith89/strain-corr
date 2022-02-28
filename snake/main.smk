# {{{1 Preamble

# {{{2 Imports

from lib.snake import (
    alias_recipe,
    alias_recipe_norelative,
    noperiod_wc,
    integer_wc,
    single_param_wc,
    limit_numpy_procs_to_1,
    curl_recipe,
    limit_numpy_procs,
    resource_calculator,
    nested_defaultdict,
    nested_dictlookup,
)
from lib.pandas_util import idxwhere
import pandas as pd
import math
from itertools import product
from textwrap import dedent as dd
from scripts.build_db import DatabaseInput
import os.path as path
from warnings import warn
import snakemake.utils

# {{{1 Configuration

# {{{2 General Configuration

snakemake.utils.min_version("6.7")


config = nested_defaultdict()

configfile: "config.yaml"


local_config_path = "config_local.yaml"
if path.exists(local_config_path):

    configfile: local_config_path


if "container" in config:

    container: config["container"]


if "MAX_THREADS" in config:
    MAX_THREADS = config["MAX_THREADS"]
else:
    MAX_THREADS = 99
if "USE_CUDA" in config:
    USE_CUDA = config["USE_CUDA"]
else:
    USE_CUDA = 0

# {{{2 Data Configuration

_mgen_meta = "meta/mgen.tsv"
if path.exists(_mgen_meta):
    _mgen = pd.read_table(_mgen_meta, index_col="mgen_id")
    config["mgen"] = {}
    for mgen_id, row in _mgen.iterrows():
        config["mgen"][mgen_id] = {}
        config["mgen"][mgen_id]["r1"] = row["filename_r1"]
        config["mgen"][mgen_id]["r2"] = row["filename_r2"]
else:
    warn(
        dd(
            f"""
            Could not load config from `{_mgen_meta}`.
            Check that path is defined and file exists.
            """
        )
    )
    config["mgen"] = {}

_mgen_x_mgen_group_meta = "meta/mgen_x_mgen_group.tsv"
if path.exists(_mgen_x_mgen_group_meta):
    _mgen_x_mgen_group = pd.read_table(_mgen_x_mgen_group_meta)
    config["mgen_group"] = {}
    for mgen_group, d in _mgen_x_mgen_group.groupby("mgen_group"):
        config["mgen_group"][mgen_group] = d.mgen_id.tolist()
else:
    warn(
        dd(
            f"""
            Could not load config from `{_mgen_x_mgen_group_meta}`.
            Check that path is defined and file exists.
            """
        )
    )
    config["mgen_group"] = {}

config['figures']['submission'] = [
]


# {{{2 Sub-pipelines


include: "snake/template.smk"
include: "snake/util.smk"
include: "snake/general.smk"
include: "snake/docs.smk"
include: "snake/mgen_preprocess.smk"


if path.exists("snake/local.smk"):

    include: "snake/local.smk"


wildcard_constraints:
    r="r|r1|r2|r3",
    group=noperiod_wc,
    mgen=noperiod_wc,
    species=noperiod_wc,
    strain=noperiod_wc,
    compound_id=noperiod_wc,
    hmm_cutoff="XX",
    model=noperiod_wc,
    param=single_param_wc,
    params=noperiod_wc,


# {{{1 Default actions


rule all:
    input:
        ["sdata/database.db"],


# {{{1 Database


database_inputs = [
    # Metadata
    DatabaseInput("subject", "smeta/subject.tsv", True),
    DatabaseInput("sample", "meta/sample.tsv", True),
    # Metagenomes
    DatabaseInput("mgen", "meta/mgen.tsv", True),
    DatabaseInput("mgen_x_mgen_group", "meta/mgen_x_mgen_group.tsv", True),
]


rule build_db:
    output:
        "sdata/database.db",
    input:
        script="scripts/build_db.py",
        schema="schema.sql",
        inputs=[entry.path for entry in database_inputs],
    params:
        args=[entry.to_arg() for entry in database_inputs],
    shell:
        dd(
            r"""
        {input.script} {output} {input.schema} {params.args}
        """
        )
