from collections import defaultdict
from warnings import warn
from numpy import ceil
import snakemake.io
import re

alias_recipe = "ln -rs {input} {output}"
alias_recipe_norelative = "ln -sT {input} {output}"
alias_fmt = lambda input, output: alias_recipe.format(input=input, output=output)
curl_recipe = "curl '{params.url}' > {output}"
curl_unzip_recipe = "curl '{params.url}' | zcat > {output}"

independent_theano_compiledir = """
        # Circumvent theano compiledir locking.
        compiledir=$(mktemp -d)
        # tar --strip-components 1 -xzf raw/theano.tgz -C $compiledir
        export THEANO_FLAGS="base_compiledir=$compiledir"
        echo $THEANO_FLAGS

        """

# Utility wildcard constrains
noperiod_wc = "[^.]+"
no_period_or_slash_wc = "[^./]+"
integer_wc = "[0-9]+"
float_noperiod_wc = "[0-9]+(e[0-9]+)?"
single_param_wc = "[^.-]+"
params_wc = noperiod_wc
endswith_period_wc = ".*\."
endswith_period_or_slash_wc = ".*[./]"


def nested_defaultdict():
    return defaultdict(nested_defaultdict)


def nested_dictlookup(mapping, *args):
    value = mapping
    for key in args:
        value = value[key]
    return value


def get_checkpoint_by_path(checkpoint, path, output_idx=0):
    regex = snakemake.io.regex(checkpoint.rule.output[output_idx])
    match = re.match(
        regex,
        path,
    )
    return checkpoint.get(**match.groupdict()), match.groupdict()
