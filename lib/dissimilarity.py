import pandas as pd
from scipy.spatial.distance import squareform, pdist
from lib.pandas_util import idxwhere
import numpy as np
import rpy2
from rpy2 import robjects
from rpy2.robjects import packages, pandas2ri
import warnings


rpy2.robjects.packages.importr("vegan")
rpy2.robjects.packages.importr("dplyr")
r_mantel_partial = rpy2.robjects.r["mantel.partial"]
r_mantel = rpy2.robjects.r["mantel"]


def dmatrix(x, **kwargs):
    return pd.DataFrame(squareform(pdist(x, **kwargs)), index=x.index, columns=x.index)


def align_dmatrices(*args, drop_labels=[]):
    master_idx = args[0].index
    arg_idx = []
    for x in args:
        x_index = set(x.index)
        x_columns = set(x.columns)
        assert x_index == x_columns
        arg_idx.append(x_index)
    shared_idx = set.intersection(*arg_idx) - set(drop_labels)
    shared_idx = idxwhere(master_idx.to_series().isin(shared_idx))
    return [x.loc[shared_idx, shared_idx] for x in args]


def mask_dmatrix(x, mask):
    "Set elements to NaN where mask is True."
    x = x.copy()
    assert (x.index == mask.index).all()
    assert (x.index == x.columns).all()
    assert (mask.index == mask.columns).all()
    x[mask] = np.nan
    return x


def triu_stack(df, k=0, columns_name_suffix="", index_name_suffix=""):
    if df.columns.name == df.index.name:
        if columns_name_suffix == index_name_suffix:
            warnings.warn(
                "Both levels of the multi-index have the same name. See *_name_suffix kwarg to fix this."
            )
        else:
            df = df.rename_axis(
                index=df.index.name + index_name_suffix,
                columns=df.columns.name + columns_name_suffix,
            )
    return df.where(np.triu(np.ones(df.shape), k=k).astype(np.bool)).stack()


def diss_table(columns_name_suffix="", index_name_suffix="", **kwargs):
    series = {}
    for k, v in kwargs.items():
        series[k] = triu_stack(
            v,
            k=1,
            columns_name_suffix=columns_name_suffix,
            index_name_suffix=index_name_suffix,
        )
    return pd.DataFrame(series)


def partial_mantel_test(
    dx, dy, dz, strata=None, method="pearson", na_rm=True, permutations=999
):
    with rpy2.robjects.conversion.localconverter(
        rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter
    ):
        x = r_mantel_partial(
            dx,
            dy,
            dz,
            strata=strata,
            method=method,
            na_rm=na_rm,
            permutations=permutations,
        )
    return x[2][0], x[3][0]


def mantel_test(dx, dy, strata=None, method="pearson", na_rm=True, permutations=999):
    with rpy2.robjects.conversion.localconverter(
        rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter
    ):
        x = r_mantel(
            dx, dy, strata=strata, method=method, na_rm=na_rm, permutations=permutations
        )
    return x[2][0], x[3][0]
