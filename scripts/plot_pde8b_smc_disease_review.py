from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from project_paths import PROJECT_DATA_ROOT


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

DONOR = "#4169A1"
PAH = "#F21A24"
INK = "#263746"
GRID = "#E5E7E9"


def fixed_jitter(n: int) -> np.ndarray:
    return np.linspace(-0.065, 0.065, n) if n > 1 else np.zeros(1)


def main() -> None:
    source = PROJECT_DATA_ROOT / "outputs" / "scrna-patient-pseudobulk-audit"
    out = PROJECT_DATA_ROOT / "outputs" / "pde8b-smc-disease-review-v1"
    out.mkdir(parents=True, exist_ok=True)
    values = pd.read_csv(source / "gse210248_pde8b_smc_patient_values.csv")
    stats = pd.read_csv(source / "pde8b_patient_level_effect_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(6.1, 3.45), sharey=True)
    colors = {"Donor": DONOR, "PAH": PAH}
    for label, ax, cell_type in zip("ab", axes, ["SMC 1", "SMC 2"]):
        z = values[values["cell_type"].eq(cell_type)].copy()
        for xpos, condition in enumerate(["Donor", "PAH"]):
            q = z[z["disease"].eq(condition)].sort_values("sample")
            x = xpos + fixed_jitter(len(q))
            ax.scatter(x, q["log_cp10k"], s=58, color=colors[condition],
                       edgecolor="white", linewidth=0.8, alpha=0.9, zorder=3)
            mean = q["log_cp10k"].mean()
            ax.hlines(mean, xpos - 0.20, xpos + 0.20, color=INK, lw=1.7, zorder=4)
        stat = stats[stats["analysis"].eq(f"GSE210248 {cell_type}")].iloc[0]
        ax.set_xticks([0, 1], ["Donor\n(n=3)", "PAH\n(n=3)"])
        ax.set_xlim(-0.45, 1.45)
        ax.set_title(cell_type, pad=8, fontweight="semibold")
        ax.grid(axis="y", color=GRID, lw=0.65, zorder=0)
        ax.text(0.04, 0.96,
                f"Hedges' g = {stat['hedges_g']:.2f}\nExact P = {stat['exact_two_sided_permutation_p']:.2f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=7.5)
        ax.text(-0.14, 1.05, label, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")
    axes[0].set_ylabel("PDE8B expression (log CP10K)")
    fig.subplots_adjust(left=0.13, right=0.98, top=0.88, bottom=0.18, wspace=0.22)

    stem = out / "Figure_PDE8B_SMC1_SMC2_disease_review_v1"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white",
                pil_kwargs={"compression": "tiff_lzw"})
    values.to_csv(out / "Figure_PDE8B_SMC1_SMC2_disease_source_data.csv", index=False)
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
