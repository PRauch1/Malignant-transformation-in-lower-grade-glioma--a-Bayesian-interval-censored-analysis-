#!/usr/bin/env python3
"""
Ridgeline visualization of the Bayesian horseshoe posterior over 99 radiomic
features (par_ratio = 0.10, full cohort, interval-censored Weibull).

Reconstructs each feature's posterior as a split-normal on the log-HR scale
using (beta_median, log(HR_q025), log(HR_q975)) — this preserves the
asymmetry characteristic of horseshoe posteriors (narrow spike near zero,
longer tail on the "alive" side).

Output: figures/fig3_horseshoe_ridgeline.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

# ----- Paths -----
HERE = Path(__file__).resolve().parent
CSV  = HERE.parent / "scripts" / "05_radiomic_posterior_full.csv"
OUT  = HERE / "figures" / "fig3_horseshoe_ridgeline.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

# ----- Load -----
df = pd.read_csv(CSV)
df = df[df["cohort"] == "Full"].copy()
assert len(df) == 99, f"Expected 99 features, got {len(df)}"

df["beta_q025"] = np.log(df["HR_q025"])
df["beta_q975"] = np.log(df["HR_q975"])
df["ci_width"]  = df["beta_q975"] - df["beta_q025"]
df["abs_med"]   = df["beta_median"].abs()

# ----- Pre-specified two-feature signature -----
SIG = {
    "original_gldm_dependencevariance_FLAIR",
    "original_glszm_largearealowgraylevelemphasis_T1c",
}

def pretty_name(f):
    """Pretty-print a PyRadiomics feature name."""
    # original_gldm_dependencevariance_FLAIR  ->  GLDM DependenceVariance · FLAIR
    parts = f.split("_")
    if parts[0] == "original":
        parts = parts[1:]
    family = parts[0].upper()
    seq    = parts[-1]
    middle = "_".join(parts[1:-1])
    # CamelCase the middle
    camel = "".join(w.capitalize() for w in middle.split("_"))
    # Insert space before internal capitals for readability
    import re
    camel_spaced = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", camel)
    return f"{family}  ·  {camel_spaced}  ·  {seq}"

# Order: widest CI (most "alive") at the top, most-shrunk at the bottom
df = df.sort_values("ci_width", ascending=True).reset_index(drop=True)
# i=0 is the most-shrunk, i=98 is the widest — we place i=0 at the BOTTOM visually

from scipy.stats import gaussian_kde
DRAWS = pd.read_csv(HERE / "horseshoe_posterior_draws.csv")
# Reorder draws to match df's feature order (already sorted by ci_width)
DRAWS = DRAWS[df["feature"].tolist()]

X_MIN, X_MAX = -0.9, 0.9
N_POINTS = 600
xs = np.linspace(X_MIN, X_MAX, N_POINTS)

densities = np.zeros((len(df), N_POINTS))
for i, feat in enumerate(df["feature"]):
    samples = DRAWS[feat].to_numpy()
    # Tiny bandwidth for spikes, scott for wider posteriors
    kde = gaussian_kde(samples, bw_method="scott")
    densities[i] = kde(xs)
densities_norm = densities / densities.max(axis=1, keepdims=True)

# ----- Colours -----
# Perceptually-uniform gradient keyed to CI width (proxy for "liveliness")
cmap_main = mpl.colormaps["viridis"]
# Re-sample so bottom (shrunk) = deep indigo, top (wide) = warm gold
# Viridis goes dark purple -> teal -> yellow, which is perfect: shrunk=dark, alive=bright
norm_rank = np.linspace(0.05, 0.92, len(df))
colors = [cmap_main(v) for v in norm_rank]

SIG_COLOR = "#E6351A"  # signal red for the pre-specified pair

# ----- Figure -----
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.edgecolor": "#2A2A2A",
    "axes.linewidth": 0.6,
    "xtick.color": "#2A2A2A",
    "ytick.color": "#2A2A2A",
})

FIG_W = 9.0
FIG_H = 12.6
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=220)
fig.patch.set_facecolor("#FBFAF7")

# Main ridgeline axes — narrow enough to leave ~30% for labels + gradient bar
ax = fig.add_axes([0.050, 0.085, 0.48, 0.895])
ax.set_facecolor("#FBFAF7")

# Marginal gradient-bar axes far to the right, with feature labels between
ax_bar = fig.add_axes([0.925, 0.085, 0.016, 0.895])
ax_bar.set_facecolor("#FBFAF7")

ROW_H = 0.88   # vertical spacing between ridge baselines
RIDGE_SCALE = 1.30  # visual height multiplier for each ridge

for i, (_, row) in enumerate(df.iterrows()):
    y0 = i * ROW_H
    dens = densities_norm[i] * RIDGE_SCALE
    is_sig = row["feature"] in SIG
    face_col = SIG_COLOR if is_sig else colors[i]
    edge_col = "#8A1A0A" if is_sig else mpl.colors.to_rgba(face_col, 0.95)
    alpha_fill = 0.88 if is_sig else 0.78

    # Fill
    ax.fill_between(
        xs, y0, y0 + dens,
        color=face_col, alpha=alpha_fill,
        linewidth=0.0, zorder=2 + (0.5 if is_sig else 0),
    )
    # Edge
    ax.plot(
        xs, y0 + dens,
        color=edge_col, linewidth=0.55 if not is_sig else 0.9,
        alpha=0.95, zorder=2.4 + (0.5 if is_sig else 0),
    )
    # Baseline
    ax.plot(
        [X_MIN, X_MAX], [y0, y0],
        color="#FBFAF7", linewidth=0.45, zorder=1.8,
    )

# Zero reference — the horseshoe's "gravity well"
ax.axvline(
    0, color="#2A2A2A", linewidth=0.8, linestyle=(0, (4, 3)),
    alpha=0.55, zorder=3,
)

# ----- Axes -----
y_top = len(df) * ROW_H + 2.2  # extra headroom so the top ridge breathes
ax.set_ylim(-0.6, y_top)
ax.set_xlim(X_MIN, X_MAX)

# X axis: show as hazard ratio, log spacing
hr_ticks = [0.5, 0.7, 1.0, 1.5, 2.0]
x_ticks = [np.log(h) for h in hr_ticks]
ax.set_xticks(x_ticks)
ax.set_xticklabels([f"{h:g}" for h in hr_ticks], fontsize=9)
ax.set_xlabel(
    "Hazard ratio for malignant transformation  (log scale)",
    fontsize=10.5, labelpad=9, color="#2A2A2A",
)

# Y axis: hide ticks, render feature names only for key rows
ax.set_yticks([])
for spine in ("top", "right", "left"):
    ax.spines[spine].set_visible(False)
ax.spines["bottom"].set_color("#A8A29E")
ax.spines["bottom"].set_linewidth(0.6)

# Label only (a) the two signature features and (b) the top 5 widest ridges.
# Use callouts with smart vertical de-overlap so the top labels never collide.
label_idx = []
for idx in df.index[-5:].tolist():
    label_idx.append(idx)
for idx in df.index[df["feature"].isin(SIG)].tolist():
    if idx not in label_idx:
        label_idx.append(idx)
label_idx = sorted(label_idx)

# Compute target y for each label (ridge's centre-of-mass)
MIN_LABEL_GAP = 1.5  # label-coords in data y
targets = [idx * ROW_H + 0.35 for idx in label_idx]
placed = []
for t in targets:
    if not placed:
        placed.append(t)
    else:
        y = max(t, placed[-1] + MIN_LABEL_GAP)
        placed.append(y)

# Connect each label to its ridge with a thin leader line
LABEL_X = X_MAX + 0.04
for idx, y_label in zip(label_idx, placed):
    row = df.iloc[idx]
    y0 = idx * ROW_H + 0.35
    is_sig = row["feature"] in SIG
    txt_col = "#8A1A0A" if is_sig else "#2A2A2A"
    lead_col = "#B0B0B0" if not is_sig else "#8A1A0A"
    prefix = "\u2605  " if is_sig else ""
    # Leader: from ridge edge to label
    ax.plot(
        [X_MAX + 0.005, LABEL_X - 0.005],
        [y0, y_label],
        color=lead_col, linewidth=0.5, alpha=0.55, zorder=4,
        clip_on=False,
    )
    ax.text(
        LABEL_X, y_label,
        prefix + pretty_name(row["feature"]),
        fontsize=7.8,
        color=txt_col,
        fontweight="bold" if is_sig else "regular",
        va="center", ha="left", zorder=5,
        clip_on=False,
    )

# ----- Marginal gradient bar showing the shrinkage axis -----
bar_grad = np.linspace(0.05, 0.92, 256).reshape(-1, 1)
ax_bar.imshow(
    bar_grad, aspect="auto", cmap=cmap_main, origin="lower",
    extent=[0, 1, -0.4, y_top], interpolation="bilinear",
)
ax_bar.set_xlim(0, 1)
ax_bar.set_ylim(-0.4, y_top)
ax_bar.set_xticks([])
ax_bar.set_yticks([])
for spine in ("top", "right", "left", "bottom"):
    ax_bar.spines[spine].set_visible(False)

# Annotations to the right of the bar
bar_txt_x = 1.35
ax_bar.text(
    bar_txt_x, y_top - 0.8,
    "Horseshoe-\nreleased",
    fontsize=8.6, color="#B8860B", ha="left", va="top",
    fontweight="bold", clip_on=False,
)
ax_bar.text(
    bar_txt_x, y_top - 3.3,
    "wider posterior,\nweaker shrinkage",
    fontsize=7.2, color="#5A5A5A", ha="left", va="top",
    style="italic", clip_on=False,
)
ax_bar.text(
    bar_txt_x, 1.6,
    "Strongly\nshrunk",
    fontsize=8.6, color="#3B1A5C", ha="left", va="bottom",
    fontweight="bold", clip_on=False,
)
ax_bar.text(
    bar_txt_x, -0.05,
    "posterior pinned\nto log(HR) = 0",
    fontsize=7.2, color="#5A5A5A", ha="left", va="bottom",
    style="italic", clip_on=False,
)

# Legend for the signature (clean, top-left inside the ridgeline axes)
leg_handles = [
    Line2D([0], [0], color=SIG_COLOR, lw=6, alpha=0.88,
           label="Pre-specified two-feature signature"),
    Line2D([0], [0], color=cmap_main(0.6), lw=6, alpha=0.85,
           label="Other radiomic features (n = 97)"),
    Line2D([0], [0], color="#2A2A2A", lw=0.9, linestyle=(0, (4, 3)),
           label="Null reference  (HR = 1)"),
]
leg = ax.legend(
    handles=leg_handles,
    loc="upper left",
    bbox_to_anchor=(0.002, 0.998),
    fontsize=8.4, frameon=False,
    handlelength=1.6, handletextpad=0.6,
    labelcolor="#2A2A2A",
)

# ----- Summary caption beneath the x-axis -----
fig.text(
    0.055, 0.018,
    "99 features  ·  4 000 posterior draws  ·  regularised-horseshoe prior  (Piironen–Vehtari, par_ratio = 0.10)",
    fontsize=7.8, color="#6A6A6A", ha="left", va="bottom",
    fontstyle="italic",
)

fig.savefig(OUT, dpi=220, facecolor=fig.get_facecolor())
print(f"Saved: {OUT}")
