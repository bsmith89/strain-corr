# {{{1 Preamble

# {{{2 Imports

from lib.snake import (
    alias_recipe,
    alias_recipe_norelative,
    noperiod_wc,
    integer_wc,
    single_param_wc,
    endswith_period_or_slash_wc,
    curl_recipe,
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

    container: config["container"]["default"]


# {{{2 Sub-pipelines


include: "snake/metadata.smk"
include: "snake/template.smk"
include: "snake/util.smk"
include: "snake/general.smk"
include: "snake/docs.smk"
include: "snake/include.smk"


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
        config["ALL"],


# {{{1 Database


database_inputs = [
    # # Metagenomes
    # DatabaseInput("mgen", "meta/mgen.tsv", True),
]


rule build_db:
    output:
        "data/database.db",
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
