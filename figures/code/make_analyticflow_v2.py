"""
Analytic-flow diagram (Supplementary Figure S2): four complementary survival models —
Primary Bayesian, Competing risks (Fine-Gray), Time-varying Cox, and Native interval-censored
(intccr) — with a cross-model CONCORDANCE banner. Matches the submitted Figure S2.
Regenerates figures_corrected/Figure_AnalyticFlow.{png,pdf}.
"""
import os
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures_corrected")
os.makedirs(OUT, exist_ok=True)
ACCENT = "#0F3A5F"; HEADER_BG = "#E8EEF4"

def box(ax, cx, cy, w, h, text, *, fc="#FFFFFF", ec=ACCENT, fs=9.0, weight="normal", color="#222"):
    r = FancyBboxPatch((cx-w/2, cy-h/2), w, h, boxstyle="round,pad=0.0,rounding_size=0.8",
                       linewidth=1.2, edgecolor=ec, facecolor=fc)
    ax.add_patch(r)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fs, fontweight=weight,
            color=color, multialignment="center")

def varrow(ax, x, y1, y2):
    ax.annotate("", xy=(x, y2), xytext=(x, y1),
                arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.1, shrinkA=4, shrinkB=4))

def oarrow(ax, x1, y1, x2, y2, via_y):
    ax.plot([x1, x1, x2, x2], [y1, via_y, via_y, y2], color="#555", lw=1.1, zorder=1)
    ax.annotate("", xy=(x2, y2), xytext=(x2, via_y+0.25*(y2-via_y)),
                arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.1, shrinkA=0, shrinkB=0))

def fig_flow():
    fig, ax = plt.subplots(figsize=(12.4, 7.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    box(ax, 50, 92.5, 62, 11,
        "Adult lower-grade glioma cohort   (n = 155)\n"
        "77 malignant-transformation events  ·  6 non-MT deaths  ·  72 censored",
        fc=ACCENT, ec=ACCENT, fs=10.5, weight="bold", color="white")

    arm_w, arm_h = 21.5, 17
    arm_y = 70.5
    arm_xs = [13.5, 38.0, 62.0, 86.5]
    arms = [
        ("Primary Bayesian", "Interval-censored Weibull\nRegularized horseshoe\n(102 coefficients)"),
        ("Competing risks", "Fine-Gray subdistribution\nmidpoint  ·  20-fold MI"),
        ("Time-varying Cox", "Counting-process\n20-fold multiple imputation\nRubin pooling"),
        ("Native interval-censored", "Fine-Gray sieve ML (intccr)\nlikelihood integrated\nover the bracket [L, R]"),
    ]
    head_h = 4.6
    for x, (head, body) in zip(arm_xs, arms):
        box(ax, x, arm_y+arm_h/2-head_h/2, arm_w, head_h, head, fc=HEADER_BG, ec=ACCENT, fs=8.3, weight="bold", color="#0F3A5F")
        box(ax, x, arm_y-head_h/2, arm_w, arm_h-head_h, body, fc="white", ec=ACCENT, fs=7.7)
    for x in arm_xs:
        oarrow(ax, 50, 86.5, x, arm_y+arm_h/2+0.5, via_y=82)

    res_y, res_h = 46, 16
    res_bodies = [
        "CE HR  22.1  (10.5 – 53.2)\nAge & volume null\n0 / 99 radiomic CrI excludes 1\nα = 0.55  (decreasing hazard)",
        "CE sHR  10.5  (5.6 – 19.5)\nGLDM-FLAIR sHR 1.54 (p = 0.006)\nGLSZM-T1c sHR 1.17  (ns)\nGray's test MT  p < 0.0001",
        "CE HR  13.3  (6.8 – 26.0)\nGTR vs Biopsy  0.32 (0.13–0.77)\nGTR vs non-GTR  0.37 (0.18–0.76)\np = 0.007",
        "CE sHR  2.94  (1.27 – 6.79)\nGLDM-FLAIR sHR 1.04  (ns)\nGLSZM-T1c sHR 1.05  (ns)\nsignature collapses to null",
    ]
    for x, body in zip(arm_xs, res_bodies):
        box(ax, x, res_y, arm_w, res_h, body, fc="#F9FAFB", ec="#BBBBBB", fs=6.9)
    for x in arm_xs:
        varrow(ax, x, arm_y-arm_h/2-0.5, res_y+res_h/2+0.5)

    tri_y = 26.5
    box(ax, 50, tri_y, 88, 9.5,
        "CONCORDANCE\n"
        "CE hazard ratio is order-of-magnitude consistent across four independent frameworks;\n"
        "the IDH-mutant subgroup amplifies the CE effect (Bayesian HR 92.2)",
        fc="#FFF7E6", ec="#B88018", fs=9.0, weight="bold", color="#6E4C00")
    for x in arm_xs:
        oarrow(ax, x, res_y-res_h/2-0.3, 50, tri_y+4.75, via_y=35)

    concl_y = 10.5
    box(ax, 50, concl_y, 94, 9.5,
        "CONCLUSIONS\n"
        "Contrast enhancement is the dominant baseline predictor  ·  GTR confers an approximately two-thirds (~68%) hazard reduction\n"
        "no single radiomic feature is credible under horseshoe  ·  two-feature signature detectable only in midpoint Fine-Gray",
        fc=ACCENT, ec=ACCENT, fs=9.0, weight="bold", color="white")
    varrow(ax, 50, tri_y-4.75-0.3, concl_y+4.75+0.3)

    box(ax, 50, 2.4, 94, 3.6,
        "Sensitivity analyses:  par_ratio grid  ·  adult-only cohort  ·  era split at 2010   →   all conclusions unchanged",
        fc="#F0EEF7", ec="#4B3F88", fs=8.5, color="#3A2F6B")

    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    fig.savefig(os.path.join(OUT, "Figure_AnalyticFlow.png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(OUT, "Figure_AnalyticFlow.pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("regenerated 4-arm analytic flow")

fig_flow()
