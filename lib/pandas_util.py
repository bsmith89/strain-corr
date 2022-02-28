def idxwhere(condition, x=None):
    if x is None:
        x = condition
    return list(x[condition].index)


def normalize_rows(df):
    return df.divide(df.sum(1), axis=0)


def align_indexes(*args):
    a0, *aa = args
    idx = set(a0.index)
    for a in aa:
        idx &= set(a.index)

    assert idx
    return [a.reindex(idx) for a in args]
