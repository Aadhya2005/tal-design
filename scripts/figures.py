"""
figures.py
----------
Plasmid map and sequence-metric figures for the design report.
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge

from codon_tools import gc_windows

INK = "#1b1b1f"
MUTED = "#6b6b76"
ACCENT = "#1f6f8b"
WARM = "#c75b39"
GREEN = "#4a7c59"
LIGHT = "#e8e6e1"


def circular_map(
    total_bp: int,
    features: list[tuple[str, int, int, str]],
    title: str,
    out_path: str,
) -> str:
    """
    features: list of (label, start_bp, end_bp, colour) with 1-based inclusive
    coordinates on the finished plasmid.
    """
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.set_aspect("equal")
    ax.axis("off")

    r_out, r_in = 1.0, 0.86

    ax.add_patch(Wedge((0, 0), r_out, 0, 360, width=r_out - r_in,
                       facecolor=LIGHT, edgecolor="none", zorder=1))

    def bp_to_deg(bp: float) -> float:
        # 12 o'clock origin, clockwise
        return 90.0 - 360.0 * (bp / total_bp)

    # Draw the arcs first, collecting label anchors.
    anchors: list[tuple[float, str]] = []
    for label, start, end, colour in features:
        a1, a2 = bp_to_deg(start - 1), bp_to_deg(end)
        ax.add_patch(Wedge((0, 0), r_out, a2, a1, width=r_out - r_in,
                           facecolor=colour, edgecolor="white",
                           linewidth=0.8, zorder=2))
        anchors.append((bp_to_deg((start - 1 + end) / 2.0), label))

    # Small features sit almost on top of each other. Pushing their labels
    # further out does not help when they share an angle, because the text
    # just slides along the same line. They are separated in ANGLE instead,
    # with a leader line back to the true position on the ring.
    MIN_SEP = 9.0
    LABEL_R = 1.16

    right = sorted([a for a in anchors if np.cos(np.deg2rad(a[0])) >= 0],
                   key=lambda t: -t[0])
    left = sorted([a for a in anchors if np.cos(np.deg2rad(a[0])) < 0],
                  key=lambda t: t[0])

    placed: list[tuple[float, float, str]] = []  # (true_deg, label_deg, text)
    for group in (right, left):
        last: float | None = None
        for deg, label in group:
            d = deg
            if last is not None:
                if group is right and d > last - MIN_SEP:
                    d = last - MIN_SEP
                elif group is left and d < last + MIN_SEP:
                    d = last + MIN_SEP
            placed.append((deg, d, label))
            last = d

    for true_deg, label_deg, label in placed:
        tr = np.deg2rad(true_deg)
        lr = np.deg2rad(label_deg)
        ha = "left" if np.cos(lr) >= 0 else "right"
        ax.plot(
            [r_out * np.cos(tr), 1.045 * np.cos(tr), (LABEL_R - 0.015) * np.cos(lr)],
            [r_out * np.sin(tr), 1.045 * np.sin(tr), (LABEL_R - 0.015) * np.sin(lr)],
            color=MUTED, linewidth=0.7, zorder=2,
        )
        ax.text(LABEL_R * np.cos(lr), LABEL_R * np.sin(lr), label,
                ha=ha, va="center", fontsize=9.5, color=INK, zorder=3)

    # tick marks every 1000 bp
    for bp in range(0, total_bp, 1000):
        rad = np.deg2rad(bp_to_deg(bp))
        ax.plot([0.82 * np.cos(rad), 0.86 * np.cos(rad)],
                [0.82 * np.sin(rad), 0.86 * np.sin(rad)],
                color=MUTED, linewidth=0.8, zorder=2)
        if bp:
            ax.text(0.76 * np.cos(rad), 0.76 * np.sin(rad), f"{bp//1000}k",
                    ha="center", va="center", fontsize=7.5, color=MUTED)

    ax.text(0, 0.06, title, ha="center", va="center",
            fontsize=13, color=INK, weight="bold")
    ax.text(0, -0.06, f"{total_bp:,} bp", ha="center", va="center",
            fontsize=10.5, color=MUTED)

    ax.set_xlim(-1.85, 1.85)
    ax.set_ylim(-1.55, 1.45)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out_path


def gc_profile(naive: str, optimised: str, out_path: str, window: int = 50) -> str:
    """Sliding-window GC before and after optimisation."""
    fig, ax = plt.subplots(figsize=(9.5, 3.4))

    for seq, label, colour in ((naive, "Naive reverse translation", MUTED),
                               (optimised, "Codon-optimised", ACCENT)):
        vals = gc_windows(seq, window=window, step=5)
        xs = [p for p, _ in vals]
        ys = [g for _, g in vals]
        ax.plot(xs, ys, color=colour, linewidth=1.4, label=label)

    ax.axhspan(30, 70, color=GREEN, alpha=0.07, zorder=0)
    ax.axhline(30, color=GREEN, linewidth=0.8, linestyle="--", alpha=0.6)
    ax.axhline(70, color=GREEN, linewidth=0.8, linestyle="--", alpha=0.6)
    ax.text(len(optimised) * 0.995, 71.5, "constraint band 30-70%",
            ha="right", va="bottom", fontsize=8, color=GREEN)

    ax.set_xlabel(f"Position in CDS (nt), {window} nt sliding window")
    ax.set_ylabel("GC (%)")
    ax.set_title("GC content before and after optimisation", fontsize=11.5, color=INK)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(15, 85)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, facecolor="white")
    plt.close(fig)
    return out_path


def codon_usage_bars(naive_stats: dict, opt_stats: dict, out_path: str) -> str:
    """Side-by-side CAI and GC comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.2))

    labels = ["Naive", "Optimised"]
    colours = [MUTED, ACCENT]

    cai_vals = [naive_stats["CAI"], opt_stats["CAI"]]
    axes[0].bar(labels, cai_vals, color=colours, width=0.55)
    axes[0].set_ylim(0, 1.0)
    axes[0].set_ylabel("CAI")
    axes[0].set_title("Codon Adaptation Index", fontsize=10.5)
    axes[0].axhline(0.8, color=WARM, linestyle="--", linewidth=0.9)
    axes[0].text(1.45, 0.81, "0.8 target", fontsize=8, color=WARM,
                 ha="right", va="bottom")
    for i, v in enumerate(cai_vals):
        axes[0].text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=9.5)

    gc_vals = [naive_stats["GC_percent"], opt_stats["GC_percent"]]
    axes[1].bar(labels, gc_vals, color=colours, width=0.55)
    axes[1].set_ylim(0, 80)
    axes[1].set_ylabel("GC (%)")
    axes[1].set_title("Overall GC content", fontsize=10.5)
    for i, v in enumerate(gc_vals):
        axes[1].text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=9.5)

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, facecolor="white")
    plt.close(fig)
    return out_path
