from collections import defaultdict
from warnings import warn
from numpy import ceil

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
integer_wc = "[0-9]+"
float_noperiod_wc = "[0-9]+(e[0-9]+)?"
single_param_wc = "[^.-]+"
params_wc = noperiod_wc


def nested_defaultdict():
    return defaultdict(nested_defaultdict)


def nested_dictlookup(mapping, *args):
    value = mapping
    for key in args:
        value = value[key]
    return value


def resource_calculator(
    baseline=1,
    threads_exponent=0,
    attempt_base=1,
    input_size_exponent=None,
    agg=None,
    bad_inputs='replace',
    **input_size_multipliers
):
    if agg is None:
        agg = sum
    if input_size_exponent is None:
        input_size_exponent = {}
    _input_size_exponent = defaultdict(lambda: 1)
    _input_size_exponent.update(input_size_exponent)
    input_size_exponent = _input_size_exponent

    def func(wildcards, input, threads, attempt):
        input_sizes = {}
        for k in input_size_multipliers:
            if hasattr(wildcards, k) and hasattr(input, k):
                warn(f"{k} is both a wildcard and an input file. Input file size takes precedence.")
            if hasattr(input, k):
                input_sizes[k] = getattr(input, k).size / 1024 / 1024
            elif hasattr(wildcards, k):
                input_sizes[k] = float(getattr(wildcards, k))

            # FIXME: Can't find a way to get around the `TypeError: '>' not supported between instances of 'TBDString' and 'int'` problem...
            # if not isinstance(input_sizes[k], (int, float)):
            #     if bad_inputs == 'fail':
            #         raise ValueError(f"Size of input '{k}' is invalid: {input_sizes[k]}")
            #     elif bad_inputs == 'warn':
            #         warn(f"Size of input '{k}' is invalid: {input_sizes[k]}")
            #     elif bad_inputs == 'ignore':
            #         pass
            #     elif bad_inputs == 'replace':
            #         warn(f"Size of input '{k}' is invalid: {input_sizes[k]}. Replacing with 1")
            #         input_sizes[k] = 1
            #     else:
            #         raise ValueError(f"{bad_inputs} is not a valid option for *bad_inputs*.")
            # else:
            #     print(input_sizes[k], file=sys.stderr)

            # # print(type(getattr(input, k).size))
            # try:
            # except AttributeError as err:
            #     warn(str(err))
            #     input_sizes[k] = 0
        base_estimate = agg(
            [baseline] + [
                input_size_multipliers[k] * input_sizes[k] ** input_size_exponent[k]
                for k in input_size_multipliers
            ]
        )
        outvalue = base_estimate * threads ** threads_exponent * attempt_base ** attempt
        # print(weighted_input_size, threads, threads_exponent, attempt_base, attempt, outvalue)
        return ceil(outvalue)

    return func
