import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.manifold import MDS, TSNE, Isomap
from scipy.spatial.distance import pdist, squareform
import numpy as np
from itertools import cycle
from sklearn.decomposition import PCA
from collections import defaultdict
import seaborn as sns
import scipy as sp
from lib.pandas_util import align_indexes


DEFAULT_MARKER_LIST = ["o", "v", "s", ">", "D", "X", "h", "^"]
DEFAULT_COLOR_LIST = [
    "black",
    "blue",
    "green",
    "orange",
    "purple",
    "red",
    "pink",
]
DEFAULT_LINESTYLE_LIST = ["-", "--", "-.", ":"]


def construct_ordered_palette(x, cm="Spectral", other="grey", extend=None, vmin=0, vmax=1, desaturate_levels=None):
    if desaturate_levels is None:
        desaturate_levels = [1.0]
    labels = pd.Series(x).unique()
    cm = mpl.cm.get_cmap(cm)
    colormap = defaultdict(lambda: other)
    if extend:
        colormap.update(extend)
    for i, (s, desat) in enumerate(zip(labels, cycle(desaturate_levels))):
        if s in colormap:
            continue
        if len(labels) == 1:
            coord = 0
        else:
            coord = i / (len(labels) - 1)
        colormap[s] = sns.desaturate(sns.saturate(cm(coord * (vmax - vmin) + vmin)), desat)
    return colormap


def construct_ordered_palette_from_list(x, colors=[], other="grey", extend=None):
    labels = pd.Series(x).unique()
    colormap = defaultdict(lambda: other)
    if extend:
        colormap = colormap.update(extend)
    for s, c in zip(labels, colors):
        if s in colormap:
            continue
        colormap[s] = c
    return colormap


def subplots_grid(ncols, naxes, ax_width=3, ax_height=2.5, **kwargs):
    nrows = int(np.ceil(naxes / ncols))
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(ax_width * ncols, ax_height * nrows), **kwargs
    )
    axs = np.reshape(axs, (nrows, ncols))
    return fig, axs


def subplots_grid_fixed_rows(nrows, naxes, ax_width=3, ax_height=2.5, **kwargs):
    ncols = int(np.ceil(naxes / nrows))
    fig, axs = plt.subplots(
        nrows, ncols, figsize=(ax_width * ncols, ax_height * nrows), **kwargs
    )
    axs = np.reshape(axs, (nrows, ncols))
    return fig, axs


def demo_palette(palette):
    df = pd.DataFrame(palete)
    print(df)


def pca_ordination(data):
    pca = PCA().fit(data)
    d1 = pca.transform(data)

    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(pca.explained_variance_ratio_, index=d1.columns)
    return d1, frac_explained, {}


def mds_ordination(data, is_dmat=False, mds_kwargs=None, pdist_kwargs=None):
    # PCoA  # TODO: Modularize out?
    if mds_kwargs is None:
        mds_kwargs = {}
    if pdist_kwargs is None:
        pdist_kwargs = {}

    if is_dmat:
        dmat = data.loc[:, data.index]  # Ensure symmetric
    else:
        dmat = pd.DataFrame(
            squareform(pdist(data, **pdist_kwargs)),
            index=data.index,
            columns=data.index,
        )
    d1 = MDS(
        n_components=2,
        max_iter=3000,
        eps=1e-12,
        random_state=1,
        dissimilarity="precomputed",
        n_jobs=1,
        **mds_kwargs,
    ).fit_transform(dmat)

    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def nmds_ordination(
    data, is_dmat=False, mds_kwargs=None, nmds_kwargs=None, pdist_kwargs=None
):
    # PCoA  # TODO: Modularize out?
    if mds_kwargs is None:
        mds_kwargs = {}
    if nmds_kwargs is None:
        nmds_kwargs = {}
    if pdist_kwargs is None:
        pdist_kwargs = {}

    if is_dmat:
        dmat = data.loc[:, data.index]  # Ensure symmetric
    else:
        dmat = pd.DataFrame(
            squareform(pdist(data, **pdist_kwargs)),
            index=data.index,
            columns=data.index,
        )

    _mds_kwargs = dict(
        max_iter=3000,
        eps=1e-9,
        random_state=1,
        dissimilarity="precomputed",
        n_jobs=1,
    )
    _mds_kwargs.update(mds_kwargs)
    init = MDS(
        n_components=2,
        **_mds_kwargs,
    ).fit_transform(dmat)

    _nmds_kwargs = dict(
        max_iter=3000,
        eps=1e-12,
        dissimilarity="precomputed",
        random_state=1,
        n_jobs=1,
        n_init=1,
    )
    _nmds_kwargs.update(nmds_kwargs)
    nmds = MDS(
        metric=False,
        n_components=2,
        **_nmds_kwargs,
    )
    d1 = nmds.fit_transform(dmat, init=init)

    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def nmds2_ordination(
    data, is_dmat=False, pca_kwargs=None, nmds_kwargs=None, pdist_kwargs=None
):
    # PCoA  # TODO: Modularize out?
    if pca_kwargs is None:
        pca_kwargs = {}
    if nmds_kwargs is None:
        nmds_kwargs = {}
    if pdist_kwargs is None:
        pdist_kwargs = {}

    if is_dmat:
        dmat = data.loc[:, data.index]  # Ensure symmetric
    else:
        dmat = pd.DataFrame(
            squareform(pdist(data, **pdist_kwargs)),
            index=data.index,
            columns=data.index,
        )

    _pca_kwargs = dict(
        tol=1e-9,
        random_state=1,
    )
    _pca_kwargs.update(pca_kwargs)
    pca = PCA(**_pca_kwargs).fit(data)
    init = pca.transform(data)

    _nmds_kwargs = dict(
        max_iter=3000,
        eps=1e-12,
        dissimilarity="precomputed",
        random_state=1,
        n_jobs=1,
        n_init=1,
    )
    _nmds_kwargs.update(nmds_kwargs)
    nmds = MDS(
        metric=False,
        n_components=2,
        **_nmds_kwargs,
    )
    d1 = nmds.fit_transform(dmat, init=init)

    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def tsne_ordination(data, is_dmat=False, tsne_kwargs=None, pdist_kwargs=None):
    # PCoA  # TODO: Modularize out?
    if tsne_kwargs is None:
        tsne_kwargs = {}
    if pdist_kwargs is None:
        pdist_kwargs = {}

    if is_dmat:
        dmat = data.loc[:, data.index]  # Ensure symmetric
    else:
        dmat = pd.DataFrame(
            squareform(pdist(data, **kwargs)),
            index=data.index,
            columns=data.index,
        )
    d1 = TSNE(n_components=2, metric="precomputed", **tsne_kwargs).fit_transform(dmat)

    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def isomap_ordination(data, is_dmat=False, ordin_kwargs=None, pdist_kwargs=None):
    if ordin_kwargs is None:
        ordin_kwargs = {}
    if pdist_kwargs is None:
        pdist_kwargs = {}
    if is_dmat:
        dmat = data.loc[:, data.index]  # Ensure symmetric
    else:
        dmat = pd.DataFrame(
            squareform(pdist(data, **pdist_kwargs)),
            index=data.index,
            columns=data.index,
        )

    d1 = Isomap(metric="precomputed", **ordin_kwargs).fit_transform(dmat)
    d1 = pd.DataFrame(
        d1,
        index=data.index,
        columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1],
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def scatterplot(
    x,
    y,
    data,
    subset=None,
    showby="__none__",
    colorby="__none__",
    color_palette=None,
    colorby_order=None,
    markerby="__none__",
    marker_palette=None,
    markerby_order=None,
    markersizeby="__none__",
    # markersize_palette=None,
    # markersizeby_order=None,
    edgecolorby="__none__",
    edgecolor_palette=None,
    edgecolorby_order=None,
    edgestyleby="__none__",
    edgestyle_palette=None,
    edgestyleby_order=None,
    zorderby="__none__",
    ax=None,
    scatter_kws={},
    fill_legend=True,
):

    data = data.copy()

    # Set up plotting defaults
    if not ax:
        fig, ax = plt.subplots()

    # Add constant column "__none__":
    data["__none__"] = "__none__"

    if colorby_order is None:
        colorby_order = data[colorby].sort_values().unique()
    if color_palette is None:
        color_palette = {k: v for k, v in zip(colorby_order, cycle(DEFAULT_COLOR_LIST))}

    if markerby_order is None:
        markerby_order = data[markerby].sort_values().unique()
    if marker_palette is None:
        marker_palette = {
            k: v for k, v in zip(markerby_order, cycle(DEFAULT_MARKER_LIST))
        }

    # if markersizeby_order is None:
    #     markersizeby_order = data[markersizeby].sort_values().unique()
    # if markersize_palette is None:
    #     markersize_palette = {k: v for k, v in zip(markersizeby_order, cycle([None]))}

    if edgecolorby_order is None:
        edgecolorby_order = data[edgecolorby].sort_values().unique()
    if edgecolor_palette is None:
        edgecolor_palette = {
            k: v for k, v in zip(edgecolorby_order, cycle(DEFAULT_COLOR_LIST))
        }

    if edgestyleby_order is None:
        edgestyleby_order = data[edgestyleby].sort_values().unique()
    if edgestyle_palette is None:
        edgestyle_palette = {
            k: v for k, v in zip(edgestyleby_order, cycle(DEFAULT_LINESTYLE_LIST))
        }

    # Default scatter options, replaced if also in scatter_kws.
    scatter_kws_ = dict(lw=2)
    scatter_kws_.update(scatter_kws)

    for (
        (
            feat_show,
            feat_color,
            feat_marker,
            feat_edgecolor,
            feat_edgestyle,
            feat_zorder,
        ),
        d,
    ) in data.groupby([showby, colorby, markerby, edgecolorby, edgestyleby, zorderby]):
        if (
            (not feat_show)
            or (feat_color not in colorby_order)
            or (feat_marker not in markerby_order)
            or (feat_edgecolor not in edgecolorby_order)
            or (feat_edgestyle not in edgestyleby_order)
        ):
            continue

        if feat_zorder == "__none__":
            feat_zorder = 0
        if markersizeby == "__none__":
            feat_markersize = None
        else:
            feat_markersize = d[markersizeby]

        ax.scatter(
            x,
            y,
            data=d,
            marker=marker_palette[feat_marker],
            c=[color_palette[feat_color]],
            label="__nolegend__",
            edgecolor=[edgecolor_palette[feat_edgecolor]],
            linestyle=edgestyle_palette[feat_edgestyle],
            s=feat_markersize,
            zorder=feat_zorder,
            **scatter_kws_,
        )

    if fill_legend:
        for feat_color in colorby_order:
            ax.scatter(
                [],
                [],
                marker="o",
                c=[color_palette[feat_color]],
                label=feat_color,
                **scatter_kws_,
            )
        for feat_marker in markerby_order:
            ax.scatter(
                [],
                [],
                marker=marker_palette[feat_marker],
                c="grey",
                label=feat_marker,
                **scatter_kws_,
            )
        # for feat_markersize in markersizeby_order:
        #     ax.scatter(
        #         [],
        #         [],
        #         marker="o",
        #         s=markersize_palette[feat_markersize],
        #         c="grey",
        #         label=feat_markersize,
        #         **scatter_kws_,
        #     )
        for feat_edgecolor in edgecolorby_order:
            ax.scatter(
                [],
                [],
                marker="o",
                c="lightgrey",
                edgecolor=edgecolor_palette[feat_edgecolor],
                label=feat_edgecolor,
                **scatter_kws_,
            )
        for feat_edgestyle in edgestyleby_order:
            ax.scatter(
                [],
                [],
                marker="o",
                c="none",
                edgecolor="black",
                linestyle=edgestyle_palette[feat_edgestyle],
                label=feat_edgestyle,
                **scatter_kws_,
            )

    return ax


def ordination_plot(
    data,
    meta,
    subset=None,
    ordin=nmds_ordination,
    ordin_kws=dict(metric="cosine"),
    xy=("PC1", "PC2"),
    frac_explained_label=True,
    **kwargs,
):
    """Plot nMDS ordination with markers colored/shaped by metadata features.

    TODO: Sane/minimal defaults for kwarg=None.
    TODO: Align this function to be similar to Seaborn plotting functions.
    """
    x, y = xy
    if subset is None:
        subset = pd.Series(True, index=data.index)
    data, meta, subset = align_indexes(data, meta, subset)
    data = data.loc[subset].copy()
    meta = meta.loc[subset].copy()

    d1, frac_explained, ordin_out = ordin(data, **ordin_kws)
    assert not d1.columns.isin(
        meta.columns
    ).any(), "Overlap between ordination dimension names and some metadata column."
    d1 = d1.join(meta)

    ax = scatterplot(x, y, d1, **kwargs)

    if frac_explained_label:
        ax.set_xlabel(f"{x} ({frac_explained[x]:0.1%})")
        ax.set_ylabel(f"{y} ({frac_explained[y]:0.1%})")
    else:
        ax.set_xlabel(f"{x}")
        ax.set_ylabel(f"{y}")
    ax.legend(bbox_to_anchor=(1, 1))
    return ax, d1, ordin_out


def rotate_xticklabels(ax=None, rotation=45, ha="right", **kwargs):
    if ax is None:
        ax = plt.gca()
    ax.set_xticklabels(
        [x.get_text() for x in ax.get_xticklabels()],
        rotation=rotation,
        ha=ha,
        **kwargs,
    )


def rotate_yticklabels(ax=None, rotation=-45, va="bottom", **kwargs):
    if ax is None:
        ax = plt.gca()
    ax.set_yticklabels(
        [y.get_text() for y in ax.get_yticklabels()],
        rotation=rotation,
        va=va,
        **kwargs,
    )


def boxplot_with_points(
    x,
    y,
    hue=None,
    data=None,
    legend=True,
    ax=None,
    dist_plotter=sns.boxplot,
    dist_kwargs={},
    points_plotter=sns.stripplot,
    points_kwargs={},
    **kwargs,
):
    if not ax:
        fig, ax = plt.subplots()

    if data is None:
        data = {}
        data["x"] = x
        data["y"] = y
        x = "x"
        y = "y"
        if hue:
            data["hue"] = hue
            hue = "hue"
        data = pd.DataFrame(data)

    data_kwargs = dict(data=data, x=x, y=y, hue=hue, **kwargs)

    distkw = dict(showfliers=False, ax=ax, saturation=0.35, whis=0)
    distkw.update(data_kwargs)
    distkw.update(dist_kwargs)
    dist_plotter(**distkw)
    handles, labels = ax.get_legend_handles_labels()

    pointskw = dict(dodge=True, linewidth=1, ax=ax)
    pointskw.update(data_kwargs)
    pointskw.update(points_kwargs)
    points_plotter(**pointskw)
    if legend:
        ax.legend(handles=handles)

    return ax


def residuals_plot(fit, ax=None, data=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()

    x = fit.predict()
    y = fit.resid

    # Add trend line, but x must be strictly increasing for UnivariateSpline
    # so here's a serious hack. :-/
    delta = lambda x: x[-1:] - x[:-1]
    min_delta = delta(np.unique(x)).min()
    x = x + np.random.uniform(-min_delta / 1000, min_delta / 1000, len(x))
    x, y = zip(*sorted(zip(x, y)))
    spline = sp.interpolate.UnivariateSpline(x, y)

    _data = fit.model.data.row_labels.to_series().to_frame()[[]].assign(x=x, y=y)
    if data is not None:
        _data = _data.join(data, rsuffix="_")

    ax.scatter("x", "y", data=_data, **kwargs)
    ax.plot(x, spline(x))
    ax.set_ylabel("Standardized Residuals")
    ax.set_xlabel("Fitted Value")


def orderplot(
    y,
    data,
    lineby=None,
    orderby=None,
    orderby_order=None,
    colorby=None,
    colorby_palette=None,
    colorby_order=None,
    xjitter=0,
    yjitter=0,
    ax=None,
    plot_kwargs={},
    legend_kwargs={},
):
    d0 = data.copy()
    if not ax:
        ax = plt.gca()

    if orderby_order is None:
        orderby_order = d0[orderby].unique()

    orderby_idx = (
        pd.Series(orderby_order, name="orderby")
        .reset_index()
        .set_index("orderby")
        .squeeze()
    )
    x = "_orderby__idx"
    d0[x] = d0[orderby].map(orderby_idx)

    if colorby_order is None:
        colorby_order = d0[colorby].unique()

    for colorby_feat, d1 in d0.groupby(colorby):
        for lineby_feat, d2 in d1.groupby(lineby):
            x_jit = np.random.uniform(0, xjitter)
            y_jit = np.random.uniform(0, yjitter)
            d3 = (
                d2.reset_index()
                .set_index(orderby)
                .reindex(orderby_order)
                .dropna(subset=[y])
                .assign(x=lambda v: v[x] + x_jit, y=lambda v: v[y] + y_jit)
            )
            ax.plot(
                "x",
                "y",
                "",
                data=d3,
                c=colorby_palette[colorby_feat],
                label="__nolegend__",
                **plot_kwargs,
            )

    for colorby_feat in colorby_order:
        ax.plot([], [], "", c=colorby_palette[colorby_feat], label=colorby_feat)
    ax.legend(**legend_kwargs, title=colorby)
    ax.set_ylabel(y)
    ax.set_xlabel(orderby)
    ax.set_xticks(orderby_idx)
    ax.set_xticklabels(orderby_idx.index)
    rotate_xticklabels(ax=ax)
    return ax

def hide_axes_and_spines(ax=None):
    if ax is None:
        ax = plt.gca()

    ax.set_axis_off()

    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    ax.yaxis.set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

def hide_axes_and_spines_except_yaxis(ax=None):
    if ax is None:
        ax = plt.gca()

    ax.xaxis.set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


def plot_stacked_barplot(data, x_var, order, palette=None, ax=None, **kwargs):
    if ax is None:
        ax = plt.subplot()
    if palette is None:
        palette = construct_ordered_palette(order)

    # Bar styles
    bar_kwargs = dict(
        width=1.0,
        alpha=1.0,
        edgecolor="k",
        lw=1,
    )
    bar_kwargs.update(kwargs)

    # Plot each bar segment
    _last_top = 0
    for y_var in order:
        ax.bar(
            x=data[x_var],
            height=data[y_var],
            bottom=_last_top,
            label=y_var,
            color=palette[y_var],
            **bar_kwargs,
        )
        _last_top += data[y_var]
    ax.set_xticks(data[x_var].values)
    return ax
