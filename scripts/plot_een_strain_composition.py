#!/usr/bin/env python3


import pandas as pd
import sfacts as sf
import matplotlib.pyplot as plt
import scipy as sp
import sys
import lib.plot
import numpy as np

SUBJECT_ORDER = ["A", "B", "H"]
SAMPLE_TYPE_ORDER = ["human", "Fermenter", "mouse"]
DROP_STRAINS_THRESH = 0.5
YLINTHRESH = 1e-4


def linkage_order(linkage, labels):
    return labels[sp.cluster.hierarchy.to_tree(linkage).pre_order(lambda x: x.id)]


def label_een_experiment_sample(x):
    if x.sample_type == "human":
        label = f"[{x.name}] {x.collection_date_relative_een_end} {x.diet_or_media}"
    elif x.sample_type in ["Fermenter_inoculum"]:
        label = f"[{x.name}] {x.source_samples} inoc {x.diet_or_media}"
    elif x.sample_type in ["Fermenter"]:
        # label = f"[{x.name}] {x.source_samples} frmnt {x.diet_or_media}"
        label = f"[{x.name}] {x.source_samples} {x.diet_or_media}"
    elif x.sample_type in ["mouse"]:
        if x.status_mouse_inflamed == "Inflamed":
            # label = f"[{x.name}] {x.source_samples} 🐭 {x.mouse_genotype} {x.diet_or_media} inflam"
            label = f"[{x.name}] {x.source_samples} {x.diet_or_media} inflam"
        elif x.status_mouse_inflamed == "not_Inflamed":
            # label = f"[{x.name}] {x.source_samples} 🐭 {x.mouse_genotype} {x.diet_or_media} not_inf"
            label = f"[{x.name}] {x.source_samples} {x.diet_or_media} not_inf"
        else:
            raise ValueError(f"sample type {x.status_mouse_inflamed} not understood")
    else:
        raise ValueError(f"sample type {x.sample_type} not understood")
    return label


if __name__ == "__main__":
    sample_inpath = sys.argv[1]
    species_depth_inpath = sys.argv[2]
    sfacts_inpath = sys.argv[3]
    species = sys.argv[4]
    outpath = sys.argv[5]

    sample = (
        pd.read_table(sample_inpath)
        .assign(
            label=lambda x: x[
                ["collection_date_relative_een_end", "diet_or_media", "sample_id"]
            ].apply(tuple, axis=1)
        )
        .set_index("sample_id")
        .assign(full_label=lambda d: d.apply(label_een_experiment_sample, axis=1))
    )

    motu_depth = (
        pd.read_table(
            species_depth_inpath,
            names=["sample", "species_id", "depth"],
            index_col=["sample", "species_id"],
        )
        .depth.unstack(fill_value=0)
        .rename(columns=str, index=lambda x: "CF_" + str(int(x.split("_")[1])))
        .rename({"CF_15": "CF_11", "CF_11": "CF_15"})  # Sample swap
    )
    motu_rabund = motu_depth.divide(motu_depth.sum(1), axis=0)

    sf_fit = (
        sf.data.World.load(
            sfacts_inpath,
        )
        .rename_coords(sample=lambda s: "CF_{}".format(int(s.split("_")[1])))
        .rename_coords(sample={"CF_11": "CF_15", "CF_15": "CF_11"})
        .drop_low_abundance_strains(0.01)
        .rename_coords(strain=str)
    )

    # Genotype similarity ordered palette:
    strain_linkage = sf_fit.genotype.linkage(optimal_ordering=True)
    strain_order = list(
        linkage_order(
            strain_linkage,
            sf_fit.strain.values,
        )
    )
    strain_order.remove("-1")  # Drop "other" strain.
    strain_palette = lib.plot.construct_ordered_palette(
        strain_order,
        cm="rainbow",
    )

    d0 = (
        sample[
            [
                "subject_id",
                "collection_date_relative_een_end",
                "sample_type",
                "diet_or_media",
                "mouse_genotype",
                "source_samples",
                "status_mouse_inflamed",
                "full_label",
            ]
        ]
        .sort_values(
            [
                "subject_id",
                "collection_date_relative_een_end",
                "sample_type",
                "diet_or_media",
                "mouse_genotype",
                "source_samples",
                "status_mouse_inflamed",
            ]
        )
        .assign(
            rabund=motu_rabund[species],
        )
        .dropna(subset=["rabund"])
    )

    _grid_sample_counts = (
        d0[["subject_id", "sample_type"]]
        .value_counts()
        .unstack()
        .reindex(columns=SAMPLE_TYPE_ORDER)
    )
    fig, axs = plt.subplots(
        nrows=len(SUBJECT_ORDER),
        ncols=len(SAMPLE_TYPE_ORDER),
        figsize=(20, 5 * len(SUBJECT_ORDER)),
        width_ratios=_grid_sample_counts.max().values,
        sharey=True,
        gridspec_kw=dict(wspace=0.1, hspace=1.5),
    )

    for subject_id, ax_row in zip(SUBJECT_ORDER, axs):
        d1 = d0[lambda x: (x.subject_id == subject_id)]
        strain_frac_sample_list = list(set(d1.index) & set(sf_fit.sample.values))
        if len(strain_frac_sample_list) == 0:
            print(f"No strain analysis for {subject_id}.")
            comm = []
            _strain_order = []
        else:
            w = (
                sf_fit.sel(sample=strain_frac_sample_list)
                .drop_low_abundance_strains(DROP_STRAINS_THRESH)
                .rename_coords(strain=str)
            )
            _strain_order = [s for s in strain_order if s in w.strain] + ["-1"]
            comm = w.community.to_pandas()
        d2 = d1.join(comm)
        for sample_type, ax in zip(SAMPLE_TYPE_ORDER, ax_row):
            d3 = d2[lambda x: (x.sample_type == sample_type)].assign(
                xpos=lambda x: np.arange(len(x.index))
            )
            ax.set_title((species, subject_id, sample_type))
            ax.scatter("xpos", "rabund", data=d3, color="k", s=10, label="__nolegend__")
            # ax.set_aspect(700, adjustable="datalim", anchor="NW")
            ax.set_ylim(-1e-3, 1)
            ax.set_yscale("symlog", linthresh=1e-4, linscale=0.1)
            ax.set_xlim(-0.5, _grid_sample_counts[sample_type].max())
            ax.set_xticks(d3.xpos)
            ax.set_xticklabels(d3.full_label)

            # Plot stacked barplot
            ax1 = ax.twinx()
            top_last = 0
            for strain in _strain_order:
                ax1.bar(
                    x="xpos",
                    height=strain,
                    data=d3,
                    bottom=top_last,
                    width=1.0,
                    alpha=1.0,
                    color=strain_palette[strain],
                    edgecolor="k",
                    lw=1,
                    label="__nolegend__",
                )
                top_last += d3[strain]
                ax.scatter(
                    [], [], color=strain_palette[strain], label=strain, marker="s", s=80
                )
            ax1.set_yticks([])
            # Put strains behind points:
            ax.set_zorder(ax1.get_zorder() + 1)  # put ax in front of ax1
            ax.patch.set_visible(False)  # hide the 'canvas'
            ax1.patch.set_visible(True)  # show the 'canvas'

            ax1.set_ylim(0, 1)
            lib.plot.rotate_xticklabels(ax=ax)
            if sample_type == "mouse":
                ax.legend(bbox_to_anchor=(1, 1))

    fig.savefig(outpath, bbox_inches="tight")
