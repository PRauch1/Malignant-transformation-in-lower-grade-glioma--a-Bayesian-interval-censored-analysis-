"""
Regenerate Figure 2 (CE hazard ratio across four models) and Figure 4 (EOR contrasts forest)
with formatting fixes:
  - Header band ("Hazard ratio" + direction hints) given more breathing room
  - "← lower hazard" / "higher hazard →" hints positioned cleanly,
    rendered only on the side of the null where data actually live
  - STR-vs-Biopsy lower CI in Figure 4 corrected from 0.33 → 0.32
    (matches table3_tvc_treatment_effects.csv: 0.79 (0.32, 1.92))
  - Slightly increased per-row height to reduce vertical crowding
  - Tightened column widths so labels don't truncate

Outputs: figures_corrected/Figure2_CE_forest.{png,pdf}
         figures_corrected/Figure4_EOR_forest.{png,pdf}
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures_corrected")
os.makedirs(OUT, exist_ok=True)

ACCENT = "#0F3A5F"
ACCENT_SEC = "#8C1D40"
MUTED = "#7F7F7F"
HEADER_BG = "#E8EEF4"
RULE = "#888888"

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
})

def draw_tabular_forest(
    rows,
    *,
    xlim,
    xticks,
    xlabel="Hazard ratio  (log scale, 95% CI)",
    figsize=(9.0, 0.46),         # slightly taller per-row
    label_width=0.36,
    numeric_width=0.22,
    pval_width=0.10,
    ref=1.0,
    filename=None,
):
    n = len(rows)
    # +1.7 to give the header band more vertical room
    fig = plt.figure(figsize=(figsize[0], 1.7 + figsize[1] * n))
    gs = gridspec.GridSpec(
        1, 4,
        width_ratios=[label_width, numeric_width,
                      1 - label_width - numeric_width - pval_width, pval_width],
        wspace=0.02,
    )
    ax_label = fig.add_subplot(gs[0, 0])
    ax_num = fig.add_subplot(gs[0, 1])
    ax_plot = fig.add_subplot(gs[0, 2])
    ax_pval = fig.add_subplot(gs[0, 3])

    # Use a slightly larger top margin so column headers + direction hints fit
    HEADER_TOP = -1.85
    DIVIDER_Y = -0.65
    DIRECTION_Y = -0.30

    for ax in (ax_label, ax_num, ax_pval):
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlim(0, 1); ax.set_ylim(n + 0.4, HEADER_TOP)

    ax_plot.set_xscale("log")
    ax_plot.set_xlim(*xlim)
    ax_plot.set_xticks(list(xticks))
    ax_plot.set_xticklabels([str(x) if x >= 1 else str(x) for x in xticks])
    ax_plot.set_ylim(n + 0.4, HEADER_TOP)
    ax_plot.set_yticks([])
    for s in ("top", "left", "right"):
        ax_plot.spines[s].set_visible(False)
    ax_plot.spines["bottom"].set_color("#222222")

    # Reference line
    ax_plot.axvline(ref, color=RULE, lw=0.8, linestyle="--", zorder=1)

    # --- Column headers (top, well above direction hints) ---
    HEADER_Y = -1.45
    ax_label.text(0.01, HEADER_Y, "Model", fontsize=9.8, fontweight="bold", color="#222")
    ax_num.text(0.5, HEADER_Y, "HR (95% CI)", fontsize=9.8, fontweight="bold", color="#222", ha="center")
    ax_pval.text(0.5, HEADER_Y, "p", fontsize=9.8, fontweight="bold", color="#222", ha="center")
    ax_plot.text(np.sqrt(xlim[0] * xlim[1]), HEADER_Y, "Hazard ratio",
                 fontsize=9.8, fontweight="bold", color="#222", ha="center",
                 transform=ax_plot.transData)

    # Header bottom rule (a bit below the direction hints)
    for ax in (ax_label, ax_num, ax_pval):
        ax.axhline(DIVIDER_Y, color="#222", lw=0.8, xmin=0.02, xmax=0.98)
    ax_plot.axhline(DIVIDER_Y, color="#222", lw=0.8, xmin=0.00, xmax=1.00,
                    clip_on=False)

    # --- Direction hints (between header label and divider rule) ---
    lo, hi = xlim
    # Render direction hints only when both sides of the null are inside the plot
    if lo < ref:
        ax_plot.text(np.sqrt(lo * ref), DIRECTION_Y, "← lower hazard",
                     fontsize=8.6, color="#2C6E49", ha="center", va="center",
                     transform=ax_plot.transData, clip_on=False, style="italic")
    if hi > ref:
        # If lo >= 1, the higher-hazard arrow goes at the geometric mean of the entire range
        anchor = np.sqrt(max(lo, ref) * hi)
        ax_plot.text(anchor, DIRECTION_Y, "higher hazard →",
                     fontsize=8.6, color="#8C1D40", ha="center", va="center",
                     transform=ax_plot.transData, clip_on=False, style="italic")

    # --- Rows ---
    for i, row in enumerate(rows):
        y = i
        if row.get("type") == "header":
            for ax in (ax_label, ax_num, ax_plot, ax_pval):
                ax.axhspan(y - 0.45, y + 0.45,
                           xmin=0.0, xmax=1.0,
                           facecolor=HEADER_BG, edgecolor="none", zorder=0)
            ax_label.text(0.01, y, row["label"],
                          fontsize=10.2, fontweight="bold", color="#0F3A5F",
                          va="center")
            if "n" in row:
                ax_num.text(0.5, y, row["n"], fontsize=9.2, color="#333",
                            ha="center", va="center", style="italic")
            continue
        ax_label.text(0.02, y, row["label"], fontsize=9.7, color="#222",
                      va="center")
        ax_num.text(0.5, y, row["hr_ci"], fontsize=9.7, color="#222",
                    ha="center", va="center", family="DejaVu Sans")
        col = row.get("color", ACCENT)
        ms = row.get("ms", 7.5)
        ax_plot.hlines(y, row["lo"], row["hi"], colors=col, linewidth=1.6, zorder=3)
        ax_plot.plot([row["lo"], row["lo"]], [y - 0.14, y + 0.14], color=col, lw=1.2)
        ax_plot.plot([row["hi"], row["hi"]], [y - 0.14, y + 0.14], color=col, lw=1.2)
        ax_plot.plot(row["hr"], y, marker="s", markersize=ms,
                     markerfacecolor=col, markeredgecolor="white",
                     markeredgewidth=0.7, zorder=4)
        ax_pval.text(0.5, y, row.get("pval", ""), fontsize=9.7,
                     ha="center", va="center", color="#222")

    ax_plot.set_xlabel(xlabel, fontsize=9.7, color="#222", labelpad=8)

    plt.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.10)
    if filename:
        fig.savefig(os.path.join(OUT, filename + ".png"), dpi=300,
                    bbox_inches="tight", facecolor="white")
        fig.savefig(os.path.join(OUT, filename + ".pdf"),
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


# === Figure 2 — CE hazard ratio across the four complementary models ===
def fig2_ce_concordance():
    rows = [
        {"type": "header", "label": "Full cohort",   "n": "n = 155, 77 events"},
        {"type": "row", "label": "Bayesian Weibull (horseshoe) \u2014 primary",
         "hr": 22.1, "lo": 10.5, "hi": 53.2, "hr_ci": "22.1  (10.5 \u2013 53.2)", "pval": "< 0.001", "color": "#0F3A5F"},
        {"type": "row", "label": "Fine-Gray competing risks \u2014 midpoint",
         "hr": 10.5, "lo": 5.62, "hi": 19.5, "hr_ci": "10.5  (5.6 \u2013 19.5)",  "pval": "< 0.001", "color": "#2E7D32"},
        {"type": "row", "label": "Cox-TVC \u2014 Rubin-pooled multiple imputation",
         "hr": 13.3, "lo": 6.83, "hi": 26.0, "hr_ci": "13.3  (6.8 \u2013 26.0)",  "pval": "< 0.001", "color": "#C8742B"},
        {"type": "row", "label": "Fine-Gray \u2014 native interval-censored (intccr)",
         "hr": 2.94, "lo": 1.27, "hi": 6.79, "hr_ci": "2.94  (1.27 \u2013 6.79)", "pval": "0.012", "color": "#8C1D40"},
        {"type": "header", "label": "IDH-mutant", "n": "n = 109, 48 events"},
        {"type": "row", "label": "Bayesian Weibull (horseshoe) \u2014 primary",
         "hr": 92.2, "lo": 26.9, "hi": 416.0, "hr_ci": "92.2  (26.9 \u2013 416)", "pval": "< 0.001", "color": "#0F3A5F"},
        {"type": "row", "label": "Cox-TVC \u2014 Rubin-pooled multiple imputation",
         "hr": 24.1, "lo": 9.66, "hi": 59.9, "hr_ci": "24.1  (9.7 \u2013 59.9)",  "pval": "< 0.001", "color": "#C8742B"},
    ]
    draw_tabular_forest(
        rows,
        xlim=(1.0, 600),
        xticks=(1, 3, 10, 30, 100, 300),
        xlabel="Hazard ratio  (log scale)",
        figsize=(9.4, 0.46),
        filename="Figure2_CE_forest",
    )

# === Figure 4 — Extent of resection ===
def fig4_eor():
    rows = [
        {"type": "header", "label": "Full cohort",   "n": "n = 155, 77 events"},
        {"type": "row", "label": "GTR vs Biopsy — primary (3-level)",
         "hr": 0.317, "lo": 0.130, "hi": 0.773, "hr_ci": "0.32  (0.13 – 0.77)", "pval": "0.012"},
        {"type": "row", "label": "GTR vs Biopsy — Model A",
         "hr": 0.251, "lo": 0.0973, "hi": 0.647, "hr_ci": "0.25  (0.10 – 0.65)", "pval": "0.004"},
        {"type": "row", "label": "GTR vs Biopsy — Model B",
         "hr": 0.369, "lo": 0.144, "hi": 0.945, "hr_ci": "0.37  (0.14 – 0.95)", "pval": "0.038"},
        {"type": "row", "label": "STR vs Biopsy — primary",
         "hr": 0.789, "lo": 0.3248, "hi": 1.92,  "hr_ci": "0.79  (0.32 – 1.92)", "pval": "0.60", "color": MUTED},
        {"type": "row", "label": "GTR vs non-GTR — binary sensitivity",
         "hr": 0.370, "lo": 0.180, "hi": 0.762, "hr_ci": "0.37  (0.18 – 0.76)", "pval": "0.007"},
        {"type": "header", "label": "IDH-mutant", "n": "n = 109, 48 events"},
        {"type": "row", "label": "GTR vs Biopsy — primary (3-level)",
         "hr": 0.159, "lo": 0.0474, "hi": 0.534, "hr_ci": "0.16  (0.05 – 0.53)", "pval": "0.003"},
        {"type": "row", "label": "STR vs Biopsy — primary",
         "hr": 0.399, "lo": 0.129, "hi": 1.23, "hr_ci": "0.40  (0.13 – 1.23)",  "pval": "0.11", "color": MUTED},
        {"type": "row", "label": "GTR vs non-GTR — binary sensitivity",
         "hr": 0.306, "lo": 0.122, "hi": 0.767, "hr_ci": "0.31  (0.12 – 0.77)", "pval": "0.012"},
    ]
    draw_tabular_forest(
        rows,
        xlim=(0.04, 3.5),
        xticks=(0.05, 0.1, 0.2, 0.5, 1, 2),
        xlabel="Hazard ratio  (log scale)",
        figsize=(9.4, 0.46),
        filename="Figure4_EOR_forest",
    )

if __name__ == "__main__":
    fig2_ce_concordance()
    fig4_eor()
    print("Wrote: Figure2_CE_forest and Figure4_EOR_forest")
