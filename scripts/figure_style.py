"""Shared publication-figure styling for project-001.

The module changes presentation only. Statistical calculations and source data
remain in the figure-specific scripts.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


COLORS = {
    "ink": "#243746",
    "muted": "#687782",
    "grid": "#D9E0E4",
    "smc1": "#2878A6",
    "smc1_light": "#DCECF4",
    "smc2": "#D97932",
    "smc2_light": "#F8E7D8",
    "teal": "#2A8C82",
    "signal": "#A64B57",
    "neutral": "#AAB4BA",
    "pale": "#F4F6F7",
}


def apply_publication_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.titleweight": "semibold",
            "axes.labelsize": 7.5,
            "axes.labelcolor": COLORS["ink"],
            "axes.edgecolor": COLORS["ink"],
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "xtick.color": COLORS["ink"],
            "ytick.color": COLORS["ink"],
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "legend.fontsize": 6.8,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.transparent": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax, label: str, x: float = -0.10, y: float = 1.10) -> None:
    fig = ax.figure
    y_fig = fig.transFigure.inverted().transform(ax.transAxes.transform((0, y)))[1]
    x_fig = 0.015 if ax.get_position().x0 < 0.50 else 0.55
    fig.text(x_fig, y_fig, label.upper(), fontsize=9, fontweight=700,
             color=COLORS["ink"], va="bottom", ha="left")


def figure_title(fig, text: str) -> None:
    """Use a compact claim title; figure numbering belongs in the legend."""
    fig.suptitle(text, x=0.01, y=0.995, ha="left", va="top",
                 fontsize=10.5, fontweight="bold", color=COLORS["ink"])


def quiet_axis(ax) -> None:
    ax.tick_params(direction="out", length=3, width=0.7)
    ax.grid(False)


def save_figure(fig, outdir: Path, stem: str, *, png_dpi: int = 300) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=png_dpi, bbox_inches="tight")
    # LZW prevents the enormous uncompressed TIFF files produced previously.
    fig.savefig(outdir / f"{stem}.tiff", dpi=600, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


apply_publication_style()
