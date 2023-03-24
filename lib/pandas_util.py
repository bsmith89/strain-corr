import pandas as pd
import subprocess


def idxwhere(condition, x=None):
    if x is None:
        x = condition
    return list(x[condition].index)


def normalize_rows(df):
    return df.divide(df.sum(1), axis=0)


def aligned_index(*args, axis="index", how="inner"):
    a0, *aa = args
    idx = set(getattr(a0, axis))
    for a in aa:
        if how == "inner":
            idx &= set(getattr(a, axis))
        elif how == "outer":
            idx |= set(getattr(a, axis))
        else:
            raise ValueError("*how* parameter must be one of {'inner', 'outer'}.")

    return idx


def align_indexes(*args, axis="index", how="inner"):
    idx = aligned_index(*args, axis=axis, how=how)
    assert idx
    return [a.reindex(idx, axis=axis) for a in args]


def repeated(df, subset=None):
    if subset is None:
        subset = df.columns
    out = pd.Series(
        (df[subset].values[1:] == df[subset].values[:-1]).all(1), index=df.index[1:]
    )
    return out.reindex(df.index, fill_value=False)


def invert_mapping(x):
    assert not x.duplicated().any()
    assert x.index.is_unique

    x_value_name = x.name
    x_index_name = x.index.name

    return x.reset_index().set_index(x_value_name)[x_index_name]


def read_table_lz4(path, *args, **kwargs):
    with subprocess.Popen(["lz4", "-dc", path], stdout=subprocess.PIPE) as proc:
        out = pd.read_table(proc.stdout, *args, **kwargs)
    return out
