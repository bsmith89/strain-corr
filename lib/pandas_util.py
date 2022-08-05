import pandas as pd


def idxwhere(condition, x=None):
    if x is None:
        x = condition
    return list(x[condition].index)


def normalize_rows(df):
    return df.divide(df.sum(1), axis=0)


def align_indexes(*args, axis="index"):
    a0, *aa = args
    idx = set(getattr(a0, axis))
    for a in aa:
        idx &= set(getattr(a, axis))

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
