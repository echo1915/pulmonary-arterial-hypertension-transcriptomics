from __future__ import annotations

import json
import sys
from pathlib import Path
import sys

_LOCAL_DEPS = Path(__file__).resolve().parents[1] / ".tmp" / "figure5-python-deps"
if _LOCAL_DEPS.exists():
    sys.path.insert(0, str(_LOCAL_DEPS))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from project_paths import OUTPUT_ROOT, PROJECT_DATA_ROOT  # noqa: E402


CANDIDATES = ["PDE8B", "PIEZO2", "SLC16A12"]
AUDIT = OUTPUT_ROOT / "pah_audit"
OUT = PROJECT_DATA_ROOT / "outputs" / "candidate-comparison"
OUT.mkdir(parents=True, exist_ok=True)

COLORS = {"PDE8B": "#C95F4B", "PIEZO2": "#3F78A8", "SLC16A12": "#38988A"}
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def load_inputs():
    priority = pd.read_csv(AUDIT / "lung_mechanism_candidate_prioritization.csv")
    priority = priority[priority.symbol.isin(CANDIDATES)].copy()
    raw = pd.read_csv(AUDIT / "gse53408_raw_cel_candidate_validation.csv")
    raw = raw[raw.symbol.isin(CANDIDATES)].copy()
    sc1 = pd.read_csv(AUDIT / "scrna_gse210248_candidate_celltype_effects.csv")
    sc1 = sc1[sc1.gene.isin(CANDIDATES)].copy()
    sc2 = pd.read_csv(AUDIT / "scrna_gse293580_candidate_effects.csv")
    sc2 = sc2[sc2.gene.isin(CANDIDATES)].copy()
    return priority, raw, sc1, sc2


def build_evidence_matrix(priority, raw, sc1, sc2):
    rows = []
    for gene in CANDIDATES:
        p = priority.loc[priority.symbol.eq(gene)].iloc[0]
        r = raw.loc[raw.symbol.eq(gene)].iloc[0]
        s = sc1.loc[sc1.gene.eq(gene)].copy()
        s["abs_g"] = s.hedges_g.abs()
        top = s.sort_values(["abs_g", "delta_log_cp10k"], ascending=[False, False]).iloc[0]
        smc = s[s.cell_type.isin(["SMC 1", "SMC 2"])].sort_values("abs_g", ascending=False)
        smc_best = smc.iloc[0]
        rep = sc2.loc[sc2.gene.eq(gene)].set_index("contrast")
        rows.append({
            "gene": gene,
            "main_meta_g": p.meta_g,
            "main_meta_fdr": p.fdr,
            "main_I2": p.I2,
            "four_bulk_direction_consistent": bool(p.same_direction),
            "PAH_vs_IPF_PH_g": p.g_PAH_vs_IPF_PH_GSE15197,
            "GSE53408_raw_g": r.hedges_g,
            "GSE53408_raw_fdr": r["adj.P.Val"],
            "GSE210248_top_cell_type": top.cell_type,
            "GSE210248_top_celltype_g": top.hedges_g,
            "GSE210248_best_SMC_subtype": smc_best.cell_type,
            "GSE210248_best_SMC_g": smc_best.hedges_g,
            "GSE210248_SMC_concordant": bool(np.sign(smc_best.hedges_g) == np.sign(p.meta_g)),
            "GSE293580_IPAH_g": rep.loc["IPAH_vs_Donor", "hedges_g"],
            "GSE293580_IPAH_concordant": bool(rep.loc["IPAH_vs_Donor", "direction_concordant_with_bulk"]),
            "GSE293580_SSc_PAH_g": rep.loc["SSc-PAH_vs_Donor", "hedges_g"],
            "GSE293580_SSc_PAH_concordant": bool(rep.loc["SSc-PAH_vs_Donor", "direction_concordant_with_bulk"]),
        })
    return pd.DataFrame(rows)


def draw_figure(evidence, priority):
    contexts = ["GSE15197", "GSE113439", "GSE254617", "GSE208592", "GSE53408 raw"]
    fig = plt.figure(figsize=(11.5, 6.8))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1.05, 1], hspace=0.52, wspace=0.34,
                            left=.08, right=.97, top=.90, bottom=.12)

    def label(panel_ax, letter):
        y_fig = fig.transFigure.inverted().transform(panel_ax.transAxes.transform((0, 1.10)))[1]
        x_fig = .015 if panel_ax.get_position().x0 < .50 else .55
        fig.text(x_fig, y_fig, letter.upper(), fontsize=10, fontweight=700,
                 va="bottom", ha="left")

    ax = fig.add_subplot(grid[:, 0])
    y0 = np.arange(len(contexts))[::-1]
    offsets = {"PDE8B": 0.22, "PIEZO2": 0.0, "SLC16A12": -0.22}
    for gene in CANDIDATES:
        p = priority.loc[priority.symbol.eq(gene)].iloc[0]
        e = evidence.loc[evidence.gene.eq(gene)].iloc[0]
        values = [p.g_GSE15197, p.g_GSE113439, p.g_GSE254617, p.g_GSE208592, e.GSE53408_raw_g]
        ax.scatter(values, y0 + offsets[gene], s=36, color=COLORS[gene], label=gene, zorder=3)
        for value, yy in zip(values, y0 + offsets[gene]):
            ax.plot([0, value], [yy, yy], lw=0.7, color=COLORS[gene], alpha=0.45)
    ax.axvline(0, color="#6B7280", lw=0.8)
    ax.set_yticks(y0, contexts)
    ax.set_xlabel("Hedges g (PAH vs control)")
    ax.set_title("Parallel lung-tissue evidence", loc="left", fontweight="bold")
    label(ax, "a")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.14))

    ax = fig.add_subplot(grid[0, 1])
    heat_cols = ["GSE210248_best_SMC_g", "GSE293580_IPAH_g", "GSE293580_SSc_PAH_g"]
    heat_names = ["GSE210248\nbest SMC", "GSE293580\nIPAH", "GSE293580\nSSc-PAH"]
    mat = evidence.set_index("gene").loc[CANDIDATES, heat_cols].to_numpy(float)
    lim = max(1.5, np.nanmax(np.abs(mat)))
    image = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(3), heat_names)
    ax.set_yticks(range(3), CANDIDATES)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(mat[i, j]) > lim * 0.55 else "#111827")
    fig.colorbar(image, ax=ax, shrink=0.72, label="Hedges g")
    ax.set_title("Patient-level single-cell evidence", loc="left", fontweight="bold")
    label(ax, "b")

    ax = fig.add_subplot(grid[1, 1])
    criteria = [
        ("4 bulk cohorts\nconcordant", evidence.four_bulk_direction_consistent),
        ("PAH vs IPF-PH\nconcordant", np.sign(evidence.PAH_vs_IPF_PH_g) == np.sign(evidence.main_meta_g)),
        ("Raw GSE53408\nFDR < 0.05", evidence.GSE53408_raw_fdr < 0.05),
        ("SMC effect\nconcordant", evidence.GSE210248_SMC_concordant),
        ("IPAH replication\nconcordant", evidence.GSE293580_IPAH_concordant),
        ("SSc-PAH replication\nconcordant", evidence.GSE293580_SSc_PAH_concordant),
    ]
    binary = np.column_stack([np.asarray(v, dtype=bool) for _, v in criteria])
    ax.imshow(binary, cmap=mpl.colors.ListedColormap(["#E5E7EB", "#2A9D8F"]), vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(criteria)), [x[0] for x in criteria], rotation=28, ha="right")
    ax.set_yticks(range(3), evidence.gene)
    for i in range(binary.shape[0]):
        for j in range(binary.shape[1]):
            ax.text(j, i, "Yes" if binary[i, j] else "No", ha="center", va="center",
                    color="white" if binary[i, j] else "#6B7280", fontweight="bold")
    ax.set_title("Predefined evidence dimensions (no composite score)", loc="left", fontweight="bold")
    label(ax, "c")

    fig.savefig(OUT / 'Supplementary_Figure_candidate_prioritization_bridge.png', dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / 'Supplementary_Figure_candidate_prioritization_bridge.tiff', dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / 'Supplementary_Figure_candidate_prioritization_bridge.pdf', bbox_inches='tight', facecolor='white')
    fig.savefig(OUT / 'Supplementary_Figure_candidate_prioritization_bridge.svg', bbox_inches='tight', facecolor='white')
    plt.close(fig)


def main():
    priority, raw, sc1, sc2 = load_inputs()
    evidence = build_evidence_matrix(priority, raw, sc1, sc2)
    evidence.to_csv(OUT / "three_candidate_evidence_matrix.csv", index=False)
    draw_figure(evidence, priority)
    summary = {
        "candidate_count": 3,
        "candidates": CANDIDATES,
        "quantitative_score_used": False,
        "reason": "Avoid post-hoc weighting; evidence dimensions are reported separately.",
        "key_observation": "All three validate in raw GSE53408; PDE8B has SMC-concordant effects and two etiologic replications, while SLC16A12 fails IPAH directional replication.",
    }
    (OUT / "three_candidate_comparison_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(evidence.to_string(index=False))


if __name__ == "__main__":
    main()

