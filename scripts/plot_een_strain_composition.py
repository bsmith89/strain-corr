#!/usr/bin/env python3


import pandas as pd
import sfacts as sf
import matplotlib.pyplot as plt
import scipy as sp
import sys
import lib.plot
from lib.pandas_util import idxwhere
import numpy as np

# SUBJECT_ORDER = ["A", "B", "H"]
# SAMPLE_TYPE_ORDER = ["human", "Fermenter", "mouse"]
DROP_STRAINS_THRESH = 0.2
OFFSET_WIDTH = 0.4

PLOT_DESIGN_NUM_PANELS = (3, 5)
FIGSIZE = (PLOT_DESIGN_NUM_PANELS[1] * 5, PLOT_DESIGN_NUM_PANELS[0] * 3.5)
SUBJECT_TRANSFER_SAMPLE_DESIGN = {
    # B – PreEEN (active disease)
    # CF-7 with direct transfer CF 379,380,381,384,385,386
    # CF7 with ex vivo samples CF 97,98,99,100 (first FR, then EEN)
    # Ex vivo with post ex vivo transfer
    # CF 97,98 with 397,406,408,409,395,402
    # CF 99,100 with 430,426,427,428,429,431
    "B": {
        "Subject B:\nTime Series": (
            0,
            [
                "CF_7",
                "CF_8",
                "CF_9",
                "CF_10",
                "CF_11",
                "CF_12",
                "CF_13",
                "CF_14",
            ],
        ),
        "Direct Transfer": (
            1,
            [
                "CF_7",
                "CF_379",
                "CF_380",
                "CF_381",
                "CF_384",
                "CF_385",
                "CF_386",
            ],
        ),
        "Ex Vivo": (
            1,
            [
                "CF_7",
                "CF_97",
                "CF_98",
                "CF_99",
                "CF_100",
            ],
        ),
        "Post Ex Vivo Transfer #1": (
            2,
            [
                "CF_97",
                "CF_98",
                "CF_397",
                "CF_406",
                "CF_408",
                "CF_409",
                "CF_395",
                "CF_402",
            ],
        ),
        "Post Ex Vivo Transfer #2": (
            2,
            [
                "CF_99",
                "CF_100",
                "CF_430",
                "CF_426",
                "CF_427",
                "CF_428",
                "CF_429",
                "CF_431",
            ],
        ),
    },
    # A – EEN (remission)
    # CF-3 with direct transfer CF 140,141,142, 149,150,151
    # CF-3 with ex vivo samples CF_103,104, 101,102 (started with EEN, then FR)
    # Ex vivo with post ex vivo transfer
    # cf 103,104 with 107-175
    # cf 101,102 with 152-157
    "A": {
        "Subject A:\nTime Series": (
            0,
            [
                "CF_1",
                "CF_2",
                "CF_3",
                "CF_6",
            ],
        ),
        "Direct Transfer": (
            1,
            [
                "CF_3",
                "CF_140",
                "CF_141",
                "CF_142",
                "CF_149",
                "CF_150",
                "CF_151",
            ],
        ),
        "Ex Vivo": (
            1,
            [
                "CF_3",
                "CF_103",
                "CF_104",
                "CF_101",
                "CF_102",
            ],
        ),
        "Post Ex Vivo Transfer #1": (
            2,
            [
                "CF_103",
                "CF_104",
                "CF_170",
                "CF_171",
                "CF_172",
                "CF_173",
                "CF_174",
                "CF_175",
            ],
        ),
        "Post Ex Vivo Transfer #2": (
            2,
            [
                "CF_101",
                "CF_102",
                "CF_152",
                "CF_153",
                "CF_154",
                "CF_155",
                "CF_156",
                "CF_157",
            ],
        ),
    },
    # H – EEN (persistent inflammation)
    # CF48 with direct transfer CF 115-120
    # CF48 with ex vivo samples 107,108 (109,110) (first EEN, then FR (no FR transfer in mice)
    # Ex vivo with post ex vivo transfer
    # 107,108 with 127-133 (chow) and 667-672 (EEN-like mouse diet; PD)
    "H": {
        "Subject H:\nTime Series": (
            0,
            [
                "CF_46",
                "CF_47",
                "CF_48",
                "CF_49",
                "CF_50",
                "CF_51",
            ],
        ),
        "Direct Transfer": (
            1,
            [
                "CF_48",
                "CF_115",
                "CF_116",
                "CF_117",
                "CF_118",
                "CF_119",
                "CF_120",
            ],
        ),
        "Ex Vivo": (
            1,
            [
                "CF_48",
                "CF_107",
                "CF_108",
                "CF_109",
                "CF_110",
            ],
        ),
        "Post Ex Vivo - Chow Diet": (
            2,
            [
                "CF_107",
                "CF_108",
                "CF_127",
                "CF_128",
                # "CF_129",  # ???  # In Debbie's email, but seems to be wrong.
                "CF_130",
                "CF_131",
                "CF_132",
                "CF_133",
            ],
        ),
        "Post Ex Vivo - EEN-Like Diet": (
            2,
            [
                "CF_107",
                "CF_108",
                "CF_667",
                "CF_668",
                "CF_669",
                "CF_670",
                "CF_671",
                "CF_672",
            ],
        ),
    },
}


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

    strain_fit = (
        sf.data.World.load(
            f"data/group/een/species/sp-{species}/r.proc.gtpro.filt-poly05-cvrg05.ss-g10000-block0-seed0.fit-sfacts48-s85-seed0.world.nc"
        )
        .rename_coords(sample=lambda s: "CF_{}".format(int(s.split("_")[1])))
        .rename_coords(sample={"CF_11": "CF_15", "CF_15": "CF_11"})
        .drop_low_abundance_strains(DROP_STRAINS_THRESH)
        .rename_coords(strain=str)
    )

    # Genotype similarity ordered palette:
    strain_linkage = strain_fit.genotype.linkage(optimal_ordering=True)
    strain_order = list(
        linkage_order(
            strain_linkage,
            strain_fit.strain.values,
        )
    )
    strain_order.remove("-1")  # Drop "other" strain.
    strain_order.append("-1")  # Add to end of list
    strain_palette = lib.plot.construct_ordered_palette(
        strain_order,
        cm="rainbow",
        extend={"-1": "lightgrey"},
    )

    d0 = (
        sample[
            [
                "subject_id",
                "collection_date_relative_een_end",
                "sample_type",
                "diet_or_media",
                "mouse_genotype",
                "timepoint",
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
            simple_label=lambda x: np.where(
                x.sample_type == "human",
                x.timepoint,
                np.where(
                    x.sample_type == "Fermenter",
                    x.diet_or_media,
                    np.where(
                        x.status_mouse_inflamed == "Inflamed",
                        x.diet_or_media + " / inflam",
                        x.diet_or_media + " / not",
                    ),
                ),
            ),
            rabund=motu_rabund[species],
        )
        .dropna(subset=["rabund"])
    )

    fig, axs = plt.subplots(
        *PLOT_DESIGN_NUM_PANELS,
        figsize=FIGSIZE,
        squeeze=False,
        sharey=True,
        gridspec_kw=dict(hspace=1.5, wspace=0),
    )

    for subject, axs_row in zip(SUBJECT_TRANSFER_SAMPLE_DESIGN, axs):
        subject_comm_sample_list = list(
            set(idxwhere(sample.subject_id == subject)) & set(strain_fit.sample.values)
        )

        try:
            subject_comm = (
                strain_fit.sel(sample=subject_comm_sample_list)
                .drop_low_abundance_strains(0.2, agg_strain_coord="-1")
                .community.to_pandas()
            )
        except ValueError:
            subject_comm = pd.DataFrame([], columns=["-1"])

        for (sample_list_label, (num_offset_samples, sample_list)), ax in zip(
            SUBJECT_TRANSFER_SAMPLE_DESIGN[subject].items(), axs_row
        ):
            d1 = (
                d0.reindex(sample_list)
                .assign(
                    t=lambda x: range(len(x)),
                )
            ).join(subject_comm)
            d1.loc[d1.index[:num_offset_samples], "t"] -= OFFSET_WIDTH  # Offset width

            lib.plot.plot_stacked_barplot(
                data=d1,
                x_var="t",
                order=[s for s in strain_order if s in subject_comm.columns],
                palette=strain_palette,
                ax=ax,
                width=0.8,
                lw=0.5,
            )

            ax.set_title(sample_list_label)
            ax.set_xticklabels(d1.simple_label)
            ax.set_aspect(4, anchor="NW")
            lib.plot.rotate_xticklabels(rotation=45, ax=ax)
            ax.set_yticks([0, 0.5, 1.0])
            ax.set_xlim(d1.t.min() - 0.5, d1.t.max() + 0.5)
        ax.legend(bbox_to_anchor=(1, 1), ncols=2)
    fig.savefig(outpath, bbox_inches="tight")

    # _grid_sample_counts = (
    #     d0[["subject_id", "sample_type"]]
    #     .value_counts()
    #     .unstack()
    #     .reindex(columns=SAMPLE_TYPE_ORDER)
    # )
    # fig, axs = plt.subplots(
    #     nrows=len(SUBJECT_ORDER),
    #     ncols=len(SAMPLE_TYPE_ORDER),
    #     figsize=(20, 5 * len(SUBJECT_ORDER)),
    #     width_ratios=_grid_sample_counts.max().values,
    #     sharey=True,
    #     gridspec_kw=dict(wspace=0.1, hspace=1.5),
    # )

    # for subject_id, ax_row in zip(SUBJECT_ORDER, axs):
    #     d1 = d0[lambda x: (x.subject_id == subject_id)]
    #     strain_frac_sample_list = list(set(d1.index) & set(sf_fit.sample.values))
    #     if len(strain_frac_sample_list) == 0:
    #         print(f"No strain analysis for {subject_id}.")
    #         comm = []
    #         _strain_order = []
    #     else:
    #         w = (
    #             sf_fit.sel(sample=strain_frac_sample_list)
    #             .drop_low_abundance_strains(DROP_STRAINS_THRESH)
    #             .rename_coords(strain=str)
    #         )
    #         _strain_order = [s for s in strain_order if s in w.strain] + ["-1"]
    #         comm = w.community.to_pandas()
    #     d2 = d1.join(comm)
    #     for sample_type, ax in zip(SAMPLE_TYPE_ORDER, ax_row):
    #         d3 = d2[lambda x: (x.sample_type == sample_type)].assign(
    #             xpos=lambda x: np.arange(len(x.index))
    #         )
    #         ax.set_title((species, subject_id, sample_type))
    #         ax.scatter("xpos", "rabund", data=d3, color="k", s=10, label="__nolegend__")
    #         # ax.set_aspect(700, adjustable="datalim", anchor="NW")
    #         ax.set_ylim(-1e-3, 1)
    #         ax.set_yscale("symlog", linthresh=1e-4, linscale=0.1)
    #         ax.set_xlim(-0.5, _grid_sample_counts[sample_type].max())
    #         ax.set_xticks(d3.xpos)
    #         ax.set_xticklabels(d3.full_label)

    #         # Plot stacked barplot
    #         ax1 = ax.twinx()
    #         top_last = 0
    #         for strain in _strain_order:
    #             ax1.bar(
    #                 x="xpos",
    #                 height=strain,
    #                 data=d3,
    #                 bottom=top_last,
    #                 width=1.0,
    #                 alpha=1.0,
    #                 color=strain_palette[strain],
    #                 edgecolor="k",
    #                 lw=1,
    #                 label="__nolegend__",
    #             )
    #             top_last += d3[strain]
    #             ax.scatter(
    #                 [], [], color=strain_palette[strain], label=strain, marker="s", s=80
    #             )
    #         ax1.set_yticks([])
    #         # Put strains behind points:
    #         ax.set_zorder(ax1.get_zorder() + 1)  # put ax in front of ax1
    #         ax.patch.set_visible(False)  # hide the 'canvas'
    #         ax1.patch.set_visible(True)  # show the 'canvas'

    #         ax1.set_ylim(0, 1)
    #         lib.plot.rotate_xticklabels(ax=ax)
    #         if sample_type == "mouse":
    #             ax.legend(bbox_to_anchor=(1, 1))

    # fig.savefig(outpath, bbox_inches="tight")
