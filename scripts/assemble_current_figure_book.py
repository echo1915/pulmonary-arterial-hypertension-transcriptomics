"""Assemble the current PAH manuscript figures and legends into one PDF."""

from pathlib import Path

from PIL import Image
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from reportlab.platypus import Image as RLImage
from reportlab.platypus import PageBreak, Paragraph, SimpleDocTemplate, Spacer, KeepTogether

from project_paths import PROJECT_DATA_ROOT


OUT_ROOT = PROJECT_DATA_ROOT / "outputs"
BOOK_DIR = OUT_ROOT / "manuscript-figure-book"
BOOK_DIR.mkdir(parents=True, exist_ok=True)
PDF_PATH = BOOK_DIR / "PAH_current_figures_and_legends_v4_main_then_supplementary.pdf"


ITEMS = [
    (
        "Figure 1 | Overall study workflow and evidence architecture",
        OUT_ROOT / "figure1-publication-v4" / "Figure_1_overall_project_workflow_v5_draft.png",
        "Four subject-level bulk-lung cohorts were integrated by cross-cohort meta-analysis to define a 188-gene robust set, followed by independent raw-CEL validation in GSE53408 and nomination of PDE8B, PIEZO2 and SLC16A12. Two single-cell RNA-sequencing cohorts provide discovery/localization (GSE210248) and independent replication (GSE293580), while a separate PDE8B-centered bulk coexpression analysis evaluates program stability and pathway context. The convergent evidence supports an SMC-associated PDE8B/cAMP working hypothesis; it does not establish causality.",
    ),
    (
        "Table 1 | Overview of bulk lung cohorts",
        OUT_ROOT / "figure1-publication-v4" / "Table_1_bulk_lung_cohort_overview_v1_draft.png",
        "Cohort roles, platforms, tissue source and audited subject counts for the four discovery cohorts and the independent GSE53408 validation cohort. Counts represent biological subjects rather than technical columns or individual cells. GSE254617 technical columns were aggregated by PHBI subject; GSE117261 was not counted as an independent primary cohort because of subject overlap.",
    ),
    (
        "Figure 2 | Cross-cohort gene-set discovery and independent lung validation",
        OUT_ROOT / "manuscript-figures-v3" / "Figure_2_gene_set_and_independent_validation_v6_draft.png",
        "a, Heatmap of Hedges' g values for the 188 robust genes across GSE15197, GSE113439, GSE254617 and GSE208592. b, Comparison of four-cohort meta-analytic effects with effects from independently reprocessed GSE53408 raw CEL files; PDE8B, PIEZO2 and SLC16A12 are highlighted. c, Cohort-specific effects for the three candidates. d, Subject-level RNA-expression distributions in GSE53408 (control, n=11; PAH, n=12), with group means. Discovery and validation analyses use subjects as biological replicates. FDR denotes false-discovery-rate adjustment.",
    ),
    (
        "Figure 3 | Single-cell landscape, cell identity and PDE8B localization",
        OUT_ROOT / "figure3-localization-review-v1" / "Figure_3_cell_atlas_markers_PDE8B_SMC_review_v3.png",
        "a, Recomputed UMAP of GSE210248 cells colored by annotated cell type. b, PDE8B expression on the same embedding. c, Marker dot plot supporting cell-type annotation; color represents relative average expression and dot size represents the detected fraction. d,e, Patient-level pseudobulk PDE8B expression in SMC1 and SMC2 for donors (n=3) and PAH patients (n=3). UMAP and marker panels are descriptive cell-level views; disease contrasts use patients, not cells, as biological replicates. Exact P values are shown without implying significance in this small cohort.",
    ),
    (
        "Figure 4 | Smooth-muscle-cell state continuum and remodeling programs",
        OUT_ROOT / "figure4-smc-main-supp-review-v1" / "Figure_4_SMC_state_continuum_main_review_v1.png",
        "a, UMAP of SMC1 and SMC2 cells. b, Paired patient-pseudobulk scores for contractile, extracellular-matrix (ECM) and remodeling programs. c,d, MYH11 and FN1 expression along the ECM-minus-contractile SMC state score; smoothed lines summarize descriptive cell-level trends. e, Mean paired SMC2-minus-SMC1 program-score differences in PAH. f, Symmetric enrichment display contrasting SMC1-associated contractile/cytoskeletal pathways with SMC2-associated ECM/remodeling pathways. These data support state association and do not demonstrate lineage transition.",
    ),
    (
        "Figure 5 | PDE8B-centered vascular and cyclic-nucleotide program",
        OUT_ROOT / "figure5-program-review-v2" / "Figure_5_PDE8B_program_main_unified_v3.png",
        "a, Cross-cohort meta-correlation of gene expression with PDE8B, based on subject-level within-cohort correlations in four bulk-lung cohorts. b, Selected pathways enriched among positive and negative core genes. c, Meta-correlations and 95% confidence intervals for predefined PAH and cyclic-nucleotide genes. d, Subject-level PDE8B expression in the independent GSE293580 strict-SMC subset (donor n=2, IPAH n=3, SSc-PAH n=4); the two-donor comparison is descriptive. e, Leave-one-cohort-out stability of the positive and negative core program; genes are features, not biological replicates. f, Patient-level cell-type pseudobulk localization in GSE210248. The figure presents associative evidence only.",
    ),
    (
        "Figure 6 | PDE8B marks a cross-cohort cyclic-nucleotide program",
        OUT_ROOT / "figure6-cyclic-nucleotide-program" / "Figure_6_PDE8B_cyclic_nucleotide_program.png",
        "a, Network of the PDE8B-associated cyclic-nucleotide program. PDE8B hydrolyses cAMP, PDE3B provides a potential cAMP-cGMP crosstalk node, and PDE5A-PRKG1 represents the cGMP-PKG branch. b, Random-effects meta-correlations and 95% confidence intervals for prespecified cyclic-nucleotide genes across four bulk-lung cohorts. Coloured estimates pass Benjamini-Hochberg FDR < 0.05; grey estimates do not. c, Cohort-specific correlations for the same genes. d, Enrichment of the positive high-confidence PDE8B core module; point size reflects pathway recall. Subjects and cohorts define biological replication. These associative data do not establish direct PDE8B regulation of cGMP or causal pathway activity.",
    ),
    (
        "Supplementary Figure 1 | Cross-layer prioritization of PDE8B among three validated candidates",
        OUT_ROOT / "candidate-comparison" / "Supplementary_Figure_candidate_prioritization_bridge.png",
        "a, Subject-level effect sizes for PDE8B, PIEZO2 and SLC16A12 across the four discovery cohorts and independent GSE53408 raw-CEL validation. b, Patient-level single-cell evidence in GSE210248 SMCs and independent GSE293580 IPAH and SSc-PAH comparisons. c, Predefined evidence dimensions shown without a post-hoc composite score. All three genes retain bulk-lung validation support. PDE8B was prioritized for Figure 3 because its direction is concordant in SMCs and in both independent etiologic comparisons; PIEZO2 lacks concordant SMC localization in GSE210248, whereas SLC16A12 lacks concordant IPAH replication in GSE293580. This prioritization does not invalidate PIEZO2 or SLC16A12 and does not establish PDE8B causality.",
    ),
    (
        "Supplementary Figure 2 | PDE8B leave-one-cohort-out sensitivity",
        OUT_ROOT / "manuscript-figures-v3" / "Supplementary_Figure_PDE8B_LOCO_v4_draft.png",
        "Sensitivity analyses evaluating whether the cross-cohort PDE8B effect and its prioritization depend on any single discovery cohort. Each leave-one-cohort-out estimate is recalculated after omitting one cohort; cohorts, rather than genes or technical columns, define the independent meta-analytic units.",
    ),
    (
        "Supplementary Figure 3 | Supporting evidence for SMC state analyses",
        OUT_ROOT / "figure4-smc-main-supp-review-v1" / "Supplementary_Figure_SMC_state_support_review_v1.png",
        "Supporting patient-level and marker-level analyses for the SMC1/SMC2 state comparison, including the robustness of contractile and ECM/remodeling patterns. Inferential summaries use patient pseudobulk; cell-level displays are descriptive and are not treated as independent biological replicates.",
    ),
    (
        "Supplementary Figure 4 | PDE8B program support and reference-gene context",
        OUT_ROOT / "figure5-program-review-v1" / "Supplementary_Figure_PDE8B_program_support_review_v1.png",
        "Additional pathway, predefined-gene and reference-gene analyses supporting interpretation of the PDE8B-centered coexpression program. Associations are estimated within cohorts at the subject level and combined across cohorts; enrichment and gene-feature distributions do not constitute independent biological replication.",
    ),
    (
        "Supplementary Figure 5 | Robustness of the PDE8B-centered program",
        OUT_ROOT / "figure5-program-review-v1" / "Supplementary_Figure_PDE8B_program_robustness_review_v1.png",
        "Leave-one-cohort-out and related sensitivity analyses for the positive and negative PDE8B coexpression cores. The analyses assess cross-cohort stability and dependence on individual cohorts, but they remain associative and do not establish PDE8B-driven regulation.",
    ),
]


def register_font():
    candidates = [
        Path("C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibri.ttf"),
    ]
    for path in candidates:
        if path.exists():
            pdfmetrics.registerFont(TTFont("BookSans", str(path)))
            bold_path = Path("C:/Windows/Fonts/arialbd.ttf")
            if bold_path.exists():
                pdfmetrics.registerFont(TTFont("BookSans-Bold", str(bold_path)))
                pdfmetrics.registerFontFamily("BookSans", normal="BookSans", bold="BookSans-Bold")
            return "BookSans"
    return "Helvetica"


def scaled_image(path: Path, max_width: float, max_height: float):
    with Image.open(path) as im:
        width, height = im.size
    scale = min(max_width / width, max_height / height)
    return RLImage(str(path), width=width * scale, height=height * scale)


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont(FONT, 7)
    canvas.setFillColor(HexColor("#66727B"))
    canvas.drawRightString(283 * mm, 8 * mm, str(doc.page))
    canvas.restoreState()


FONT = register_font()
PAGE = landscape(A4)
doc = SimpleDocTemplate(
    str(PDF_PATH), pagesize=PAGE,
    leftMargin=14 * mm, rightMargin=14 * mm,
    topMargin=12 * mm, bottomMargin=14 * mm,
    title="PAH current figures and legends",
    author="project-001",
)

legend_style = ParagraphStyle(
    "Legend", fontName=FONT, fontSize=8.4, leading=10.8,
    textColor=HexColor("#20262D"), alignment=TA_JUSTIFY,
)

story = []
usable_width = PAGE[0] - doc.leftMargin - doc.rightMargin
for index, (title, path, legend) in enumerate(ITEMS):
    if not path.exists():
        raise FileNotFoundError(path)
    if "candidate_prioritization_bridge" in path.name:
        figure_height = 145 * mm
    elif "Figure_6_PDE8B" in path.name:
        figure_height = 148 * mm
    else:
        figure_height = 158 * mm
    figure = scaled_image(path, usable_width, figure_height)
    story.append(KeepTogether([
        figure,
        Spacer(1, 3 * mm),
        Paragraph(f"<b>{title}.</b> {legend}", legend_style),
    ]))
    if index != len(ITEMS) - 1:
        story.append(PageBreak())

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(PDF_PATH)
