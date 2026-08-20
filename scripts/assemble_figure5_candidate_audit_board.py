"""Assemble a source-level audit board for Figure 5 candidates.

This is a review artifact, not a submission figure.  It deliberately preserves each
source image in full so that panel provenance and duplication can be assessed before
the final Figure 5 is redesigned.
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

from project_paths import PROJECT_DATA_ROOT


DATA_ROOT = PROJECT_DATA_ROOT
OUTPUT_ROOT = DATA_ROOT / "outputs"
AUDIT_DIR = OUTPUT_ROOT / "figure5-candidate-audit-v1"


CANDIDATES = [
    {
        "id": "S1",
        "title": "Prior F5: coexpression, pathways, etiologic replication",
        "category": "prior_F5",
        "path": OUTPUT_ROOT / "submission-draft-figures" / "Figure_5_coexpression_pathways_replication.png",
    },
    {
        "id": "S2",
        "title": "Prior dense F5: vascular program and exploratory interactions",
        "category": "prior_F5",
        "path": OUTPUT_ROOT / "manuscript-figures-dense-v2" / "Figure_5_PDE8B_vascular_program_and_interactions.png",
    },
    {
        "id": "S3",
        "title": "New candidate: robust-signature to PDE8B bridge",
        "category": "new_candidate",
        "path": OUTPUT_ROOT / "signature-pde8b-bridge" / "Figure_signature_PDE8B_bridge_draft.png",
    },
    {
        "id": "S4",
        "title": "New candidate: pathway-level bridge audit",
        "category": "new_candidate",
        "path": OUTPUT_ROOT / "pathway-bridge" / "Figure_pathway_bridge_draft.png",
    },
    {
        "id": "S5",
        "title": "New candidate: GSE293580 whole-atlas annotation QC",
        "category": "new_candidate_QC",
        "path": OUTPUT_ROOT / "gse293580-reanalysis" / "gse293580_phase1_umap_draft.png",
    },
    {
        "id": "S6",
        "title": "New candidate: GSE293580 mural-cell annotation QC",
        "category": "new_candidate_QC",
        "path": OUTPUT_ROOT / "gse293580-reanalysis" / "gse293580_phase2_mural_umap_draft.png",
    },
    {
        "id": "S7",
        "title": "New candidate: integrated coexpression and SMC review",
        "category": "new_candidate_mixed",
        "path": OUTPUT_ROOT / "pde8b-smc-integrated-review-v1" / "Figure_PDE8B_SMC_integrated_review_v1.png",
    },
    {
        "id": "S8",
        "title": "Prior F5 state-remodeling draft (now overlaps tentative F4)",
        "category": "overlaps_F4",
        "path": OUTPUT_ROOT / "six-main-figure-draft-v1" / "Figure_5_SMC_state_program_remodeling.png",
    },
]


MODULE_AUDIT = [
    {
        "module_id": "M1",
        "evidence_module": "Cross-cohort PDE8B coexpression landscape",
        "source_ids": "S1;S2",
        "statistical_unit": "Subjects within four lung cohorts; cohort-level random-effects meta-analysis",
        "key_evidence": "15,015 common genes; 543 stable mechanism-module genes; 201 high-confidence core genes",
        "key_limit": "Association and program coupling only; no causal regulation",
        "provisional_destination": "MAIN",
        "redesign_note": "Retain meta-correlation/FDR landscape; show uncertainty and keep labels sparse",
    },
    {
        "module_id": "M2",
        "evidence_module": "PDE8B core-module pathway enrichment",
        "source_ids": "S1;S2",
        "statistical_unit": "Genes/pathway sets (not biological replicates)",
        "key_evidence": "Contractile, cytoskeletal, cyclic-nucleotide and ciliary themes",
        "key_limit": "Database-dependent and partly redundant; ciliary prominence requires narrative justification",
        "provisional_destination": "MAIN_IF_SPACE",
        "redesign_note": "Separate positive and negative modules; use restrained dot plot and fewer terms",
    },
    {
        "module_id": "M3",
        "evidence_module": "Connections to established PAH/cAMP vascular genes",
        "source_ids": "S2",
        "statistical_unit": "Subjects within cohorts; cohort-level meta-correlations",
        "key_evidence": "Includes PRKG1, PDE5A and other established vascular/PAH genes",
        "key_limit": "Selected-gene panel must have an a priori selection rule",
        "provisional_destination": "MAIN",
        "redesign_note": "Replace decorative lollipop with estimates and confidence intervals or cohort-consistency display",
    },
    {
        "module_id": "M4",
        "evidence_module": "Independent GSE293580 patient-level replication",
        "source_ids": "S1;S5;S6",
        "statistical_unit": "Patients; strict SMC has donor n=2, IPAH n=3, SSc-PAH n=4",
        "key_evidence": "PDE8B direction is elevated in strict SMC; program shifts are direction-level only",
        "key_limit": "Very low n; exact tests and adjusted module tests are non-significant; donor PDE8B values are zero",
        "provisional_destination": "MAIN_DESCRIPTIVE",
        "redesign_note": "Show every patient point and detection status; do not use bars alone or emphasize unstable Hedges g",
    },
    {
        "module_id": "M5",
        "evidence_module": "Robust 188-gene disease signature to PDE8B bridge",
        "source_ids": "S3",
        "statistical_unit": "Genes",
        "key_evidence": "Disease effect and PDE8B coexpression direction correlate (rho=0.365); 67.7% direction concordance",
        "key_limit": "No membership enrichment; only 6 of 186 mapped robust genes overlap stable PDE8B module",
        "provisional_destination": "SUPPLEMENT",
        "redesign_note": "Keep the direction-coupling scatter; remove flow-card styling and avoid inherited-gene-set language",
    },
    {
        "module_id": "M6",
        "evidence_module": "Pathway-level bridge between robust signature and PDE8B program",
        "source_ids": "S4",
        "statistical_unit": "Enriched pathway terms",
        "key_evidence": "Descriptive convergence at broad biological-theme level",
        "key_limit": "Zero exact shared pathway terms; theme matching is interpretive",
        "provisional_destination": "SUPPLEMENT_COMBINE_WITH_M5",
        "redesign_note": "Do not retain as a standalone main panel",
    },
    {
        "module_id": "M7",
        "evidence_module": "GSE293580 whole-atlas and mural annotation QC",
        "source_ids": "S5;S6",
        "statistical_unit": "Cells for visualization; patients for downstream inference",
        "key_evidence": "Shows atlas integration, mural refinement, disease-group and patient mixing",
        "key_limit": "UMAP is annotation/QC evidence, not disease-effect inference",
        "provisional_destination": "SUPPLEMENT",
        "redesign_note": "Retain only after marker-score and donor-mixing QC are finalized",
    },
    {
        "module_id": "M8",
        "evidence_module": "Exploratory SMC-pericyte ligand-receptor expression proxies",
        "source_ids": "S2",
        "statistical_unit": "Patients; primary comparison has donor n=2 and PAH n=6",
        "key_evidence": "Selected axes show directional differences",
        "key_limit": "No primary family reaches FDR significance (best q=0.417); proxy is not communication inference",
        "provisional_destination": "SUPPLEMENT_OPTIONAL_OR_DELETE",
        "redesign_note": "Exclude from main F5 unless the manuscript retains an explicitly exploratory interaction section",
    },
    {
        "module_id": "M9",
        "evidence_module": "Integrated 12-gene network, patient/cohort heatmaps and SMC trends",
        "source_ids": "S7",
        "statistical_unit": "Mixed: cohorts, patients and cells",
        "key_evidence": "Cohort consistency and patient-level pseudobulk patterns are visible",
        "key_limit": "Network selection is decorative; cell trends overlap tentative F4 and PDE8B trend is weak",
        "provisional_destination": "SUPPLEMENT_PARTIAL",
        "redesign_note": "Keep cohort/patient consistency only; delete network and duplicated cell-level state trends",
    },
    {
        "module_id": "M10",
        "evidence_module": "Prior F5 SMC state-remodeling composite",
        "source_ids": "S8",
        "statistical_unit": "Patients/subjects in source analyses",
        "key_evidence": "Contractile loss and remodeling-state gain",
        "key_limit": "Scientific content is already absorbed by tentative F4; schematic is conclusion-card styling",
        "provisional_destination": "REMOVE_FROM_F5",
        "redesign_note": "Do not duplicate in F5 or create another supplementary copy",
    },
]


def _check_inputs() -> None:
    missing = [str(item["path"]) for item in CANDIDATES if not item["path"].exists()]
    if missing:
        raise FileNotFoundError("Missing source images:\n" + "\n".join(missing))


def _write_manifest() -> Path:
    manifest_path = AUDIT_DIR / "Figure_5_candidate_source_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_id", "category", "title", "source_path"])
        writer.writeheader()
        for item in CANDIDATES:
            writer.writerow(
                {
                    "source_id": item["id"],
                    "category": item["category"],
                    "title": item["title"],
                    "source_path": str(item["path"]),
                }
            )
    return manifest_path


def _write_module_audit() -> Path:
    audit_path = AUDIT_DIR / "Figure_5_candidate_module_audit_v1.csv"
    fieldnames = list(MODULE_AUDIT[0])
    with audit_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(MODULE_AUDIT)
    return audit_path


def _font(size: int, bold: bool = False) -> ImageFont.ImageFont:
    font_name = "arialbd.ttf" if bold else "arial.ttf"
    font_path = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts" / font_name
    if font_path.exists():
        return ImageFont.truetype(str(font_path), size=size)
    return ImageFont.load_default()


def _assemble_board() -> tuple[Path, Path]:
    canvas_width, canvas_height = 2600, 3600
    margin_x, top_y, bottom_margin = 60, 230, 70
    gap_x, gap_y = 55, 45
    cell_width = (canvas_width - 2 * margin_x - gap_x) // 2
    cell_height = (canvas_height - top_y - bottom_margin - 3 * gap_y) // 4
    title_height = 58
    border_width = 7
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    draw.text(
        (margin_x, 50),
        "Figure 5 candidate-source audit board",
        font=_font(40, bold=True),
        fill="#111820",
    )
    draw.text(
        (margin_x, 105),
        "Review artifact; not submission artwork. Whole source figures are shown without cropping.",
        font=_font(23),
        fill="#3f4a54",
    )
    draw.text(
        (margin_x, 145),
        "Use source IDs to decide main figure, supplement, or removal.",
        font=_font(23),
        fill="#3f4a54",
    )

    border_colors = {
        "prior_F5": "#5B8DB8",
        "new_candidate": "#D08C3F",
        "new_candidate_QC": "#6A9C78",
        "new_candidate_mixed": "#8D6AA8",
        "overlaps_F4": "#9B9B9B",
    }

    for index, item in enumerate(CANDIDATES):
        row, column = divmod(index, 2)
        x0 = margin_x + column * (cell_width + gap_x)
        y0 = top_y + row * (cell_height + gap_y)
        x1, y1 = x0 + cell_width, y0 + cell_height
        border_color = border_colors[item["category"]]
        draw.rounded_rectangle(
            (x0, y0, x1, y1),
            radius=12,
            outline=border_color,
            width=border_width,
            fill="white",
        )
        draw.text(
            (x0 + 18, y0 + 13),
            f'{item["id"]}  {item["title"]}',
            font=_font(20, bold=True),
            fill="#111820",
        )
        image_box = (
            x0 + 15,
            y0 + title_height,
            x1 - 15,
            y1 - 15,
        )
        source = Image.open(item["path"]).convert("RGB")
        contained = ImageOps.contain(
            source,
            (image_box[2] - image_box[0], image_box[3] - image_box[1]),
            method=Image.Resampling.LANCZOS,
        )
        paste_x = image_box[0] + (image_box[2] - image_box[0] - contained.width) // 2
        paste_y = image_box[1] + (image_box[3] - image_box[1] - contained.height) // 2
        canvas.paste(contained, (paste_x, paste_y))

    draw.text(
        (margin_x, canvas_height - 44),
        "Key: blue prior F5; orange new evidence; green annotation/QC; purple mixed candidate; grey overlaps tentative F4.",
        font=_font(19),
        fill="#3f4a54",
    )

    png_path = AUDIT_DIR / "Figure_5_candidate_source_audit_board_v1.png"
    pdf_path = AUDIT_DIR / "Figure_5_candidate_source_audit_board_v1.pdf"
    canvas.save(png_path, dpi=(220, 220), optimize=True)
    canvas.save(pdf_path, "PDF", resolution=220.0)
    return png_path, pdf_path


def main() -> None:
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    _check_inputs()
    manifest_path = _write_manifest()
    module_audit_path = _write_module_audit()
    png_path, pdf_path = _assemble_board()
    print(f"Audit board PNG: {png_path}")
    print(f"Audit board PDF: {pdf_path}")
    print(f"Source manifest: {manifest_path}")
    print(f"Module audit: {module_audit_path}")


if __name__ == "__main__":
    os.environ.setdefault("MPLBACKEND", "Agg")
    main()
