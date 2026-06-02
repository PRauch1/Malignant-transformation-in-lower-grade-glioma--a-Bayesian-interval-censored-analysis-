#!/usr/bin/env python3
"""Full-cohort swimmer plot for Figure 1.

Shows all 155 patients: 77 transformers (MT, top panel) + 78 non-transformers
(right-censored, bottom panel). Preserves the editorial look of the original
'Heterogeneity of Disease Evolution in LGG' figure while extending the argument
to the whole cohort.

Data: data/ttm_baseline.csv + data/event_timing.csv
  - event == 1: MT observed in interval [L, R]
  - event == 0: no MT; L = last MRI (right-censored); R = inf
  - t_surgery from event_timing.csv marks surgical intervention
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

# --------------------------------------------------------------------------- #
# Style
# --------------------------------------------------------------------------- #
mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.titlesize": 12,
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

INK = "#1F2733"
MUTED = "#6B7280"
INDOLENT = "#2F5A8C"      # cool blue
TRANSFORM = "#C14B36"     # warm crimson
CENSOR = "#8A9BA8"        # cool grey-blue for censored bars
SURGERY = "#1B1F23"

# Diverging cmap: short therapeutic window = warm, long = cool
tw_cmap = LinearSegmentedColormap.from_list(
    "therapeutic_window",
    ["#9A1A1A", "#D85B3A", "#F2B872", "#F5F2E0", "#9FC6D8", "#3E7FA5", "#1F3E66"],
    N=256,
)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)


def load_data():
    base = pd.read_csv(os.path.join(DATA, "ttm_baseline.csv"))
    tim = pd.read_csv(os.path.join(DATA, "event_timing.csv"))
    df = base.merge(tim, on="patient_number", how="left")
    df["R_num"] = df["R"].replace("inf", np.inf).astype(float)
    # For non-transformers, R = inf; replace with L for plotting (drawn to L + arrow)
    return df


def main():
    df = load_data()
    mt = df[df["event"] == 1.0].copy()
    nc = df[df["event"] == 0.0].copy()

    # ---- MT panel: sort by R descending so indolent (long-latency) are at top,
    #       early transformers at bottom. Color by therapeutic window.
    mt["tw"] = mt["R_num"] - mt["t_surgery"]
    mt = mt.sort_values("R_num", ascending=False).reset_index(drop=True)

    # ---- Non-MT panel: sort by follow-up length L descending so longest at top
    nc["follow_up"] = nc["L"]
    nc = nc.sort_values("follow_up", ascending=False).reset_index(drop=True)

    # --------------------------------------------------------------------- #
    # Figure layout: two vertically stacked swimmer panels with shared x.
    # Heights proportional to N so rows are visually comparable.
    # --------------------------------------------------------------------- #
    FIG_W = 12.0
    ROW_H = 0.072       # compressed rows — figure prints on one page
    H_MT = ROW_H * len(mt) + 0.45
    H_NC = ROW_H * len(nc) + 0.45
    HEAD_PAD = 1.15
    GAP_PAD = 0.85
    FOOT_PAD = 1.85
    FIG_H = H_MT + H_NC + HEAD_PAD + GAP_PAD + FOOT_PAD

    fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor="white")

    # Axes rectangles in figure-fraction
    top_h = H_MT / FIG_H
    bot_h = H_NC / FIG_H
    head_h = HEAD_PAD / FIG_H
    gap_h = GAP_PAD / FIG_H
    foot_h = FOOT_PAD / FIG_H

    LEFT = 0.11
    RIGHT = 0.81
    WIDTH = RIGHT - LEFT

    top_y = 1.0 - head_h - top_h
    bot_y = foot_h + 0.005
    ax_top = fig.add_axes([LEFT, top_y, WIDTH, top_h])
    ax_bot = fig.add_axes([LEFT, bot_y, WIDTH, bot_h])

    XMIN = 0.02   # ~1 week
    XMAX = 22.0   # years
    for ax in (ax_top, ax_bot):
        ax.set_xscale("log")
        ax.set_xlim(XMIN, XMAX)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.spines["left"].set_color("#CFCFCF")
        ax.spines["bottom"].set_color(INK)
        ax.tick_params(axis="y", length=0)
        ax.grid(which="major", axis="x", color="#ECECEC", lw=0.6, zorder=0)

    # -------- TOP PANEL: transformers ------------------------------------- #
    norm = Normalize(vmin=-1.0, vmax=8.0)
    for i, row in mt.iterrows():
        y = len(mt) - 1 - i  # i=0 is longest R (indolent), goes to top
        L = max(row["L"], XMIN * 0.9)
        R = row["R_num"]
        if not np.isfinite(R):
            continue
        col = tw_cmap(norm(row["tw"]))
        # Bar: diagnosis (X=XMIN) to R, with highlighted [L, R] interval end
        ax_top.plot([XMIN, R], [y, y], color=col, lw=2.2, alpha=0.95,
                    solid_capstyle="round", zorder=3)
        # Interval bracket [L, R] subtly distinguished (slightly thicker)
        ax_top.plot([L, R], [y, y], color=col, lw=3.6, alpha=0.8,
                    solid_capstyle="round", zorder=4)
        # Surgery dot
        if np.isfinite(row["t_surgery"]) and row["t_surgery"] > 0:
            ax_top.plot(row["t_surgery"], y, marker="o", markersize=3.2,
                        markerfacecolor=SURGERY, markeredgecolor="white",
                        markeredgewidth=0.5, zorder=6)
        # R marker (MT observed here)
        ax_top.plot(R, y, marker="|", markersize=7, color=col, zorder=5,
                    markeredgewidth=1.6)

    ax_top.set_ylim(-1.2, len(mt))
    ax_top.set_yticks([])
    ax_top.set_ylabel(f"Transformers  (n = {len(mt)})",
                      fontsize=10.5, color=INK, labelpad=8)
    ax_top.set_xticklabels([])
    ax_top.tick_params(axis="x", length=0)

    # Phenotype call-outs on the right (indolent at top, early-transforming at bottom)
    call_x = XMAX * 1.10
    ax_top.text(call_x, len(mt) - 3, "Indolent phenotype",
                fontsize=10.5, fontweight="bold", color=INDOLENT,
                va="center", ha="left", clip_on=False)
    ax_top.text(call_x, len(mt) - 6.0,
                "extended stability\nbefore late transformation",
                fontsize=8.8, color=MUTED, va="top", ha="left",
                clip_on=False)
    ax_top.text(call_x, 5.0, "Early transforming",
                fontsize=10.5, fontweight="bold", color=TRANSFORM,
                va="center", ha="left", clip_on=False)
    ax_top.text(call_x, 2.8,
                "short latency\nto malignancy",
                fontsize=8.8, color=MUTED, va="top", ha="left",
                clip_on=False)

    # Panel letter + title
    fig.text(LEFT - 0.075, top_y + top_h - 0.008, "a",
             fontsize=17, fontweight="bold", color=INK)
    fig.text(LEFT, top_y + top_h + 0.008,
             "Patients with documented malignant transformation",
             fontsize=12, fontweight="bold", color=INK)

    # -------- BOTTOM PANEL: non-transformers ------------------------------ #
    # Last-MRI bar with right-censor arrow after the last observation
    arrow_end = XMAX * 0.95
    for i, row in nc.iterrows():
        y = len(nc) - 1 - i  # plot longest-followed at top
        L = max(row["L"], XMIN * 0.9)
        ax_bot.plot([XMIN, L], [y, y], color=CENSOR, lw=2.0, alpha=0.90,
                    solid_capstyle="round", zorder=3)
        # last MRI tick
        ax_bot.plot(L, y, marker="|", markersize=6.5, color=CENSOR,
                    markeredgewidth=1.4, zorder=5)
        # right-censor arrow — short arrow pointing right to suggest continuation
        arr = FancyArrowPatch(
            (L, y), (min(L * 1.55, arrow_end), y),
            arrowstyle="-|>", mutation_scale=7,
            color=CENSOR, lw=1.0, alpha=0.75, zorder=4,
            shrinkA=0, shrinkB=0,
        )
        ax_bot.add_patch(arr)
        # Surgery dot
        if np.isfinite(row["t_surgery"]) and row["t_surgery"] > 0:
            ax_bot.plot(row["t_surgery"], y, marker="o", markersize=3.0,
                        markerfacecolor=SURGERY, markeredgecolor="white",
                        markeredgewidth=0.5, zorder=6)

    ax_bot.set_ylim(-1.2, len(nc))
    ax_bot.set_yticks([])
    ax_bot.set_ylabel(f"Non-transformers  (n = {len(nc)})",
                      fontsize=10.5, color=INK, labelpad=8)

    # Custom x ticks in "clinical time" on the bottom axis only
    xticks = [1/12, 3/12, 1, 3, 10]
    xlabels = ["1 mo", "3 mo", "1 yr", "3 yrs", "10 yrs"]
    ax_bot.set_xticks(xticks)
    ax_bot.set_xticklabels(xlabels)
    ax_bot.set_xlabel("Time from diagnosis  (log scale)",
                      fontsize=10.5, color=INK, labelpad=10)
    # Drop a "Diagnosis" annotation out in the left margin, not over the axis
    fig.text(LEFT - 0.06, bot_y - 0.012, "Diagnosis",
             fontsize=9, color=MUTED, ha="left", va="top")

    # Callout for censored panel
    ax_bot.text(call_x, len(nc) * 0.62,
                "Stable disease",
                fontsize=10.5, fontweight="bold", color=INDOLENT,
                va="center", ha="left", clip_on=False)
    ax_bot.text(call_x, len(nc) * 0.62 - 4.5,
                "right-censored\nat last surveillance MRI",
                fontsize=8.8, color=MUTED, va="top", ha="left",
                clip_on=False)

    fig.text(LEFT - 0.075, bot_y + bot_h - 0.008, "b",
             fontsize=17, fontweight="bold", color=INK)
    fig.text(LEFT, bot_y + bot_h + 0.008,
             "Patients without documented transformation (right-censored)",
             fontsize=12, fontweight="bold", color=INK)

    # -------- Legend strip at bottom (spans x) ---------------------------- #
    legend_y = 0.010
    legend_h = 0.045
    legend_ax = fig.add_axes([LEFT, legend_y, WIDTH, legend_h])
    legend_ax.set_xlim(0, 1); legend_ax.set_ylim(0, 1)
    legend_ax.axis("off")

    # Colour scale (therapeutic window)
    scale_l = 0.00; scale_w = 0.46
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    legend_ax.imshow(grad, aspect="auto", cmap=tw_cmap,
                     extent=(scale_l, scale_l + scale_w, 0.30, 0.70),
                     zorder=2)
    legend_ax.add_patch(Rectangle((scale_l, 0.30), scale_w, 0.40,
                                  fill=False, edgecolor="#BBB", lw=0.8))
    for frac, lab in [(0.0, "reactive (<0)"), (0.2, "0"),
                      (0.40, "1 yr"),
                      (1.0, "> 8 yrs")]:
        x = scale_l + frac * scale_w
        legend_ax.plot([x, x], [0.30, 0.70], color=INK, lw=0.6)
        legend_ax.text(x, 0.15, lab, fontsize=8.5, ha="center", va="top",
                       color=INK)
    legend_ax.text(scale_l + scale_w / 2, 0.95,
                   "Therapeutic window (surgery → malignancy)",
                   fontsize=9, ha="center", va="bottom",
                   color=INK, fontstyle="italic")

    # Surgery / censor legend glyphs
    sl_x = 0.58
    legend_ax.plot(sl_x + 0.01, 0.55, marker="o", markersize=4.5,
                   markerfacecolor=SURGERY, markeredgecolor="white",
                   markeredgewidth=0.6)
    legend_ax.text(sl_x + 0.03, 0.55, "Surgical intervention",
                   fontsize=9, va="center", color=INK)

    legend_ax.plot([sl_x + 0.005, sl_x + 0.045], [0.15, 0.15],
                   color=CENSOR, lw=2.2)
    arr2 = FancyArrowPatch((sl_x + 0.045, 0.15), (sl_x + 0.075, 0.15),
                           arrowstyle="-|>", mutation_scale=6,
                           color=CENSOR, lw=1.0)
    legend_ax.add_patch(arr2)
    legend_ax.text(sl_x + 0.085, 0.15, "Right-censored (ongoing follow-up)",
                   fontsize=9, va="center", color=INK)

    # -------- Figure title ----------------------------------------------- #
    title_y = 1.0 - 0.018
    fig.text(0.5, title_y,
             "Heterogeneity of disease evolution across the full LGG cohort",
             fontsize=14.5, fontweight="bold", color=INK,
             ha="center", va="top")
    fig.text(0.5, title_y - 0.023,
             f"Transformers (n = {len(mt)}) and non-transformers (n = {len(nc)}) on a shared log time axis",
             fontsize=10.5, color=MUTED, ha="center", va="top",
             fontstyle="italic")

    out_png = os.path.join(OUT, "fig0_kinetic_trajectories.png")
    out_jpg = os.path.join(OUT, "fig0_kinetic_trajectories.jpg")
    out_pdf = os.path.join(OUT, "fig0_kinetic_trajectories.pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(out_jpg, dpi=300, bbox_inches="tight", facecolor="white")
    try:
        fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
        print("Wrote:", out_pdf)
    except PermissionError:
        print("Skipped PDF (existing file locked):", out_pdf)
    plt.close(fig)
    print("Wrote:", out_png)
    print("      ", out_jpg)


if __name__ == "__main__":
    main()
