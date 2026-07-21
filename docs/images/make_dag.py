#!/usr/bin/env python3
"""Render the AutoAFIDs pipeline overview (docs/images/dag.png).

A hand-laid, all-features-enabled view of the Snakemake rulegraph that mirrors
the rules in ``autoafids/workflow`` (Snakefile + rules/*.smk). Colours and the
rounded-box style follow Snakemake's own ``--rulegraph`` output. Regenerate with:

    python docs/images/make_dag.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).with_name("dag.png")

# node -> (x, y, edge colour). y increases upward; solid rule boxes get a
# saturated outline, input files stay black-on-white.
NODES = {
    # inputs
    "T2w":               (1.3, 9.0, "#111111"),
    "FLAIR":             (3.2, 9.0, "#111111"),
    "ct":                (5.0, 9.0, "#111111"),
    "T1w":               (6.9, 9.0, "#111111"),
    # rules
    "synthsr":           (3.2, 7.5, "#8E44AD"),
    "n4_bias_correction":(6.4, 7.5, "#C4B400"),
    "norm_im":           (3.2, 6.0, "#2C7FB8"),
    "regmni2sub":        (6.4, 6.0, "#5FB316"),
    "resample_im":       (3.2, 4.5, "#B01E2E"),
    "mni2subfids":       (6.2, 4.5, "#CE7B1E"),
    "download_cnn_model":(9.2, 4.5, "#22B455"),
    "applyfidmodel":     (5.0, 3.0, "#17BEBB"),
    "regqc":             (2.9, 1.5, "#9CBE1F"),
    "stereotaxy":        (5.0, 1.5, "#C0157F"),
    "fidqc":             (7.1, 1.5, "#B4C11E"),
    "all":               (5.0, 0.0, "#22C7E0"),
}

# per-node half-width, sized to the label so text never overflows the box
HALFW = {n: max(0.72, 0.24 + 0.058 * len(n)) for n in NODES}

INPUTS = {"T2w", "FLAIR", "ct", "T1w"}

# (src, dst, style): "solid" | "dashed"
EDGES = [
    ("T2w", "synthsr", "dashed"),
    ("FLAIR", "synthsr", "dashed"),
    ("ct", "synthsr", "dashed"),
    ("T1w", "n4_bias_correction", "dashed"),
    ("synthsr", "norm_im", "solid"),
    ("synthsr", "regmni2sub", "solid"),
    ("n4_bias_correction", "norm_im", "solid"),
    ("n4_bias_correction", "regmni2sub", "solid"),
    ("norm_im", "resample_im", "solid"),
    ("regmni2sub", "mni2subfids", "solid"),
    ("resample_im", "applyfidmodel", "solid"),
    ("mni2subfids", "applyfidmodel", "solid"),
    ("download_cnn_model", "applyfidmodel", "solid"),
    ("applyfidmodel", "regqc", "solid"),
    ("applyfidmodel", "stereotaxy", "solid"),
    ("applyfidmodel", "fidqc", "solid"),
    ("regqc", "all", "solid"),
    ("stereotaxy", "all", "solid"),
    ("fidqc", "all", "solid"),
]

BOX_H = 0.42                       # half-height in data units
ARROW = "#8A8A8A"

fig, ax = plt.subplots(figsize=(8.2, 9.4))
ax.set_xlim(-0.3, 10.8)
ax.set_ylim(-0.9, 9.8)
ax.axis("off")


def anchor(node, side):
    x, y, _ = NODES[node]
    w = HALFW[node]
    return {"top": (x, y + BOX_H), "bottom": (x, y - BOX_H),
            "left": (x - w, y), "right": (x + w, y)}[side]


def draw_edge(src, dst, style, rad=0.0):
    p0 = anchor(src, "bottom")
    p1 = anchor(dst, "top")
    ax.add_patch(FancyArrowPatch(
        p0, p1, connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-|>", mutation_scale=13, lw=1.6,
        color=ARROW, linestyle=style, shrinkA=1, shrinkB=3,
        zorder=1, capstyle="round"))


# straight edges
for src, dst, style in EDGES:
    draw_edge(src, dst, style)

# applyfidmodel -> all, curving left past the QC row
draw_edge("applyfidmodel", "all", "solid", rad=0.42)

# the "if SynthSR" branch: norm_im feeds the detector directly for non-T1w
p0 = anchor("norm_im", "left")
p1 = anchor("applyfidmodel", "left")
ax.add_patch(FancyArrowPatch(
    p0, p1, connectionstyle="arc3,rad=-0.55", arrowstyle="-|>",
    mutation_scale=13, lw=1.6, color=ARROW, linestyle="dashed",
    shrinkA=3, shrinkB=3, zorder=1, capstyle="round"))
ax.text(0.55, 5.15, "if SynthSR", fontsize=12, style="italic",
        color="#333333", ha="center", va="center")

# nodes on top
for name, (x, y, ec) in NODES.items():
    w = HALFW[name]
    lw = 2.4 if name in INPUTS else 2.6
    ax.add_patch(FancyBboxPatch(
        (x - w, y - BOX_H), 2 * w, 2 * BOX_H,
        boxstyle="round,pad=0.02,rounding_size=0.18",
        mutation_aspect=1.0, fc="white", ec=ec, lw=lw, zorder=2))
    ax.text(x, y, name, ha="center", va="center", zorder=3,
            fontsize=11.5, color="#1a1a1a", family="sans-serif")

fig.savefig(OUT, dpi=220, bbox_inches="tight", facecolor="white")
print(f"wrote {OUT}")
