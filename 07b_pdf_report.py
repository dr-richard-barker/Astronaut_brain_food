#!/usr/bin/env python3
"""
Script 07b: Generate the Phylo-branded PDF report.

Uses ReportLab per the pdf-report-generation skill guidelines.
Produces: /mnt/results/report_astronaut_opposite_forcing.pdf

Sections:
  1. Title page
  2. Executive Summary
  3. Methods
  4. Results (consensus signatures, LINCS reversers, enrichment, nutrient-gene map, vegan recipes)
  5. Discussion
  6. Conclusions & Next Steps
  7. References
"""

import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch, cm
from reportlab.lib.colors import HexColor, white, black
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, Image, KeepTogether, HRFlowable, ListFlowable, ListItem
)
from reportlab.platypus.flowables import Flowable
from reportlab.graphics.shapes import Drawing, Rect, String, Line
from reportlab.graphics.charts.barcharts import VerticalBarChart, HorizontalBarChart

# ── Phylo brand colors ──
PHYLO_BLACK     = HexColor("#000000")
PHYLO_WARM_GRAY = HexColor("#ECE9E2")
PHYLO_OFF_WHITE = HexColor("#FAF9F3")
PHYLO_LIME      = HexColor("#E9ED4C")
PHYLO_ORANGE    = HexColor("#FF9400")
PHYLO_GREEN     = HexColor("#75A025")
PHYLO_PINK      = HexColor("#FD9BED")
PHYLO_BLUE      = HexColor("#0279EE")
PHYLO_GOLD      = HexColor("#D4A04A")

HEADING_COLOR   = HexColor("#111111")
BODY_TEXT       = HexColor("#2C2A26")
MUTED_TEXT      = HexColor("#8A8378")
CAPTION_TEXT    = MUTED_TEXT
TABLE_HEADER_BG = PHYLO_GOLD
TABLE_HEADER_FG = HexColor("#FFFFFF")
TABLE_ALT_ROW   = HexColor("#F9F7F3")
TABLE_BORDER    = HexColor("#D5CFC5")
DIVIDER_COLOR   = PHYLO_GOLD
CALLOUT_BG      = PHYLO_OFF_WHITE
CALLOUT_BORDER  = PHYLO_GOLD
LINK_COLOR      = HexColor("#0563C1")

CHART_COLORS = [PHYLO_GOLD, PHYLO_BLUE, PHYLO_GREEN, PHYLO_ORANGE, PHYLO_PINK, PHYLO_BLACK]

FONT_HEADING = "Helvetica-Bold"
FONT_BODY    = "Helvetica"
FONT_ITALIC  = "Helvetica-Oblique"
FONT_MONO    = "Courier"

PROJ = "/workspace/astronaut-opposite-forcing"
FIG_DIR = os.path.join(PROJ, "results", "figures")
OUTPUT_PATH = "/mnt/results/report_astronaut_opposite_forcing.pdf"

REPORT_TITLE = "Astronaut Opposite Forcing: Spaceflight Transcriptomic Signatures and Vegan Nutritional Countermeasures"


# ──────────────────────────────────────────────────────────
# Styles
# ──────────────────────────────────────────────────────────
styles = getSampleStyleSheet()

styles.add(ParagraphStyle(name="ReportTitle", fontName=FONT_HEADING,
    fontSize=22, textColor=HEADING_COLOR, spaceBefore=0, spaceAfter=6, leading=28))

styles.add(ParagraphStyle(name="Subtitle", fontName=FONT_BODY,
    fontSize=12, textColor=PHYLO_GOLD, spaceAfter=4))

styles.add(ParagraphStyle(name="Attribution", fontName=FONT_ITALIC,
    fontSize=10, textColor=MUTED_TEXT, spaceAfter=8))

styles.add(ParagraphStyle(name="SectionHead", fontName=FONT_HEADING,
    fontSize=16, textColor=HEADING_COLOR, spaceBefore=24, spaceAfter=10, leading=20))

styles.add(ParagraphStyle(name="SubSectionHead", fontName=FONT_HEADING,
    fontSize=12, textColor=HEADING_COLOR, spaceBefore=14, spaceAfter=6, leading=16))

styles.add(ParagraphStyle(name="Body", fontName=FONT_BODY,
    fontSize=10.5, textColor=BODY_TEXT, alignment=TA_JUSTIFY,
    spaceAfter=8, leading=15))

styles.add(ParagraphStyle(name="BodyLeft", fontName=FONT_BODY,
    fontSize=10.5, textColor=BODY_TEXT, alignment=TA_LEFT,
    spaceAfter=8, leading=15))

styles.add(ParagraphStyle(name="Caption", fontName=FONT_ITALIC,
    fontSize=9, textColor=CAPTION_TEXT, alignment=TA_CENTER,
    spaceBefore=4, spaceAfter=14))

styles.add(ParagraphStyle(name="TableHeader", fontName=FONT_HEADING,
    fontSize=9, textColor=TABLE_HEADER_FG, alignment=TA_LEFT, leading=12))

styles.add(ParagraphStyle(name="TableCell", fontName=FONT_BODY,
    fontSize=8.5, textColor=BODY_TEXT, alignment=TA_LEFT, leading=11))

styles.add(ParagraphStyle(name="TableCellSmall", fontName=FONT_BODY,
    fontSize=7.5, textColor=BODY_TEXT, alignment=TA_LEFT, leading=10))

styles.add(ParagraphStyle(name="CalloutText", fontName=FONT_BODY,
    fontSize=10, textColor=BODY_TEXT, alignment=TA_LEFT, leading=14))

styles.add(ParagraphStyle(name="RefText", fontName=FONT_BODY,
    fontSize=9, textColor=BODY_TEXT, alignment=TA_LEFT, leading=12, spaceAfter=4))

styles.add(ParagraphStyle(name="BulletText", fontName=FONT_BODY,
    fontSize=10.5, textColor=BODY_TEXT, alignment=TA_LEFT, leading=15, spaceAfter=4))


# ──────────────────────────────────────────────────────────
# Page header/footer
# ──────────────────────────────────────────────────────────
def page_header_footer(canvas, doc):
    canvas.saveState()
    w, h = letter
    canvas.setFont(FONT_BODY, 9)
    canvas.setFillColor(MUTED_TEXT)
    canvas.drawString(60, h - 40, "Astronaut Opposite Forcing — Spaceflight Signatures & Nutritional Countermeasures")
    canvas.setStrokeColor(PHYLO_GOLD)
    canvas.setLineWidth(1)
    canvas.line(60, h - 48, w - 60, h - 48)
    canvas.setStrokeColor(TABLE_BORDER)
    canvas.setLineWidth(0.75)
    canvas.line(60, 40, w - 60, 40)
    canvas.setFont(FONT_BODY, 8)
    canvas.setFillColor(MUTED_TEXT)
    canvas.drawCentredString(w / 2, 26, f"Page {doc.page}")
    canvas.restoreState()


# ──────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────
def divider(width=480):
    return HRFlowable(width=width, thickness=1, color=DIVIDER_COLOR, spaceAfter=10, spaceBefore=4)


def callout_box(text, width=460):
    data = [[Paragraph(text, styles["CalloutText"])]]
    t = Table(data, colWidths=[width])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), CALLOUT_BG),
        ("BOX", (0, 0), (-1, -1), 0.5, TABLE_BORDER),
        ("LINEBEFORE", (0, 0), (0, -1), 3, PHYLO_GOLD),
        ("TOPPADDING", (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
        ("LEFTPADDING", (0, 0), (-1, -1), 14),
        ("RIGHTPADDING", (0, 0), (-1, -1), 14),
    ]))
    t.hAlign = "CENTER"
    return t


def make_table(headers, rows, col_widths, header_style=None, cell_style=None,
               max_rows=None, repeat_header=True):
    """Build a styled table with gold header and alternating rows."""
    if header_style is None:
        header_style = styles["TableHeader"]
    if cell_style is None:
        cell_style = styles["TableCell"]
    if max_rows and len(rows) > max_rows:
        rows = rows[:max_rows]

    data = [[Paragraph(f'<b>{h}</b>', header_style) for h in headers]]
    for row in rows:
        data.append([Paragraph(str(c), cell_style) for c in row])

    t = Table(data, colWidths=col_widths, repeatRows=1 if repeat_header else 0)
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), TABLE_HEADER_BG),
        ("TEXTCOLOR", (0, 0), (-1, 0), TABLE_HEADER_FG),
        ("GRID", (0, 0), (-1, -1), 0.5, TABLE_BORDER),
        ("BOX", (0, 0), (-1, -1), 0.75, TABLE_BORDER),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]
    for i in range(2, len(data), 2):
        style_cmds.append(("BACKGROUND", (0, i), (-1, i), TABLE_ALT_ROW))
    t.setStyle(TableStyle(style_cmds))
    t.hAlign = "CENTER"
    return t


def embed_figure(filename, caption_text, doc_width, max_height=350):
    """Embed a PNG figure with caption, preserving aspect ratio."""
    img_path = os.path.join(FIG_DIR, filename)
    from PIL import Image as PILImage
    pil_img = PILImage.open(img_path)
    iw, ih = pil_img.size
    aspect = ih / iw

    target_w = doc_width
    target_h = target_w * aspect
    if target_h > max_height:
        target_h = max_height
        target_w = target_h / aspect

    img = Image(img_path, width=target_w, height=target_h)
    img.hAlign = "CENTER"
    caption = Paragraph(caption_text, styles["Caption"])
    return KeepTogether([img, Spacer(1, 4), caption])


# ──────────────────────────────────────────────────────────
# Load data
# ──────────────────────────────────────────────────────────
def load_data():
    d = {}
    d["liver_sig"] = pd.read_csv(os.path.join(PROJ, "data/processed/consensus_signature_liver_immune.csv"))
    d["brain_sig"] = pd.read_csv(os.path.join(PROJ, "data/processed/consensus_signature_brain.csv"))
    d["liver_lincs"] = pd.read_csv(os.path.join(PROJ, "results/tables/lincs_liver_immune/tier1_ranking_annotated.csv"))
    d["brain_lincs"] = pd.read_csv(os.path.join(PROJ, "results/tables/lincs_brain/tier1_ranking_annotated.csv"))
    d["liver_gsea_h"] = pd.read_csv(os.path.join(PROJ, "results/tables/gsea_liver_immune_hallmark.csv"))
    d["brain_gsea_h"] = pd.read_csv(os.path.join(PROJ, "results/tables/gsea_brain_hallmark.csv"))
    d["nutrient_map"] = pd.read_csv(os.path.join(PROJ, "results/tables/nutrient_gene_map_combined.csv"))
    d["recipes"] = pd.read_csv(os.path.join(PROJ, "results/tables/vegan_recipes_scored.csv"))
    d["meal_plans"] = pd.read_csv(os.path.join(PROJ, "results/tables/vegan_meal_plans.csv"))
    d["countermeasure"] = pd.read_csv(os.path.join(PROJ, "results/tables/countermeasure_nutrient_summary.csv"))

    # Screening summary
    with open(os.path.join(PROJ, "results/tables/lincs_screening_summary.json")) as f:
        d["screening_summary"] = json.load(f)

    # Ortholog mapping stats
    d["ortholog"] = pd.read_csv(os.path.join(PROJ, "data/processed/ortholog_mapping.csv"))

    # Study catalog
    d["study_catalog"] = pd.read_csv(os.path.join(PROJ, "data/study_catalog.csv"))
    d["download_status"] = pd.read_csv(os.path.join(PROJ, "data/download_status.csv"))

    return d


# ──────────────────────────────────────────────────────────
# Build report
# ──────────────────────────────────────────────────────────
def build_report():
    data = load_data()
    date_str = datetime.now().strftime("%B %d, %Y")

    doc = SimpleDocTemplate(
        OUTPUT_PATH, pagesize=letter,
        topMargin=52, bottomMargin=52,
        leftMargin=60, rightMargin=60,
        title="Astronaut Opposite Forcing Report",
        author="Biomni / Phylo",
    )
    story = []
    dw = doc.width  # document content width

    # ═══════════════════════════════════════════════════════
    # TITLE PAGE
    # ═══════════════════════════════════════════════════════
    story.append(Spacer(1, 50))
    story.append(Paragraph("Astronaut Opposite Forcing", styles["ReportTitle"]))
    story.append(Paragraph("Spaceflight Transcriptomic Signatures and Vegan Nutritional Countermeasures", styles["Subtitle"]))
    story.append(Spacer(1, 8))
    story.append(Paragraph(f"<i>Generated by Biomni  |  {date_str}</i>", styles["Attribution"]))
    story.append(Spacer(1, 20))
    story.append(divider(dw))
    story.append(Spacer(1, 10))

    # Key metrics callout
    liver_up = (data["liver_sig"]["direction"] == "up").sum()
    liver_dn = (data["liver_sig"]["direction"] == "down").sum()
    liver_core = (data["liver_sig"]["core"] == True).sum()
    brain_up = (data["brain_sig"]["direction"] == "up").sum()
    brain_dn = (data["brain_sig"]["direction"] == "down").sum()
    brain_core = (data["brain_sig"]["core"] == True).sum()
    n_studies = len(data["study_catalog"])
    n_downloaded = (data["download_status"]["status"] == "downloaded").sum() if "status" in data["download_status"].columns else 31

    metrics_text = (
        f"<b>Studies screened:</b> {n_studies} NASA OSDR rodent spaceflight RNA-seq studies<br/>"
        f"<b>Studies meta-analyzed:</b> 15 (with valid DE signal)<br/>"
        f"<b>Liver/Immune signature:</b> {liver_up + liver_dn} consensus genes ({liver_up} up, {liver_dn} down), {liver_core} core<br/>"
        f"<b>Brain signature:</b> {brain_up + brain_dn} consensus genes ({brain_up} up, {brain_dn} down), {brain_core} core<br/>"
        f"<b>LINCS reversers:</b> 219 liver, 221 brain candidate compounds<br/>"
        f"<b>Vegan countermeasure nutrients:</b> 10 bioactive compounds<br/>"
        f"<b>Curated recipes:</b> 12 phase-optimized vegan meals"
    )
    story.append(callout_box(metrics_text, width=dw))
    story.append(Spacer(1, 20))

    # Abstract
    story.append(Paragraph("Abstract", styles["SubSectionHead"]))
    abstract = (
        "Prolonged spaceflight induces systemic transcriptomic changes in astronauts, including "
        "mitochondrial dysfunction, cell-cycle suppression, metabolic reprogramming, and immune dysregulation. "
        "We built a reproducible computational pipeline (<i>astronaut-opposite-forcing</i>) that (1) constructs "
        "consensus spaceflight transcriptomic signatures from NASA OSDR rodent ISS RNA-seq studies via "
        "random-effects meta-analysis, (2) screens the LINCS L1000 perturbation compendium for compounds that "
        "reverse these signatures (<i>opposite forcings</i>), (3) annotates hits with DrugBank, Broad Drug "
        "Repurposing Hub, and ChEMBL, (4) validates biological plausibility via functional enrichment "
        "(GSEA/ORA) and JensenLab DISEASES, and (5) translates therapeutic candidates into a vegan nutrition "
        "layer using FooDB and NutriChem 2.0. The liver/immune signature (525 genes) is dominated by "
        "downregulated E2F/G2M cell-cycle programs and upregulated fatty acid and xenobiotic metabolism. "
        "The brain signature (505 genes) shows profound downregulation of oxidative phosphorylation and "
        "mitochondrial complex I assembly. Top reversers include Imatinib and Metformin (liver), and "
        "Vemurafenib and Formoterol (brain). Vegan nutritional countermeasures featuring sulforaphane, "
        "curcumin, quercetin, and EGCG provide accessible, food-based adjuncts to pharmacological intervention."
    )
    story.append(Paragraph(abstract, styles["Body"]))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # EXECUTIVE SUMMARY
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("1. Executive Summary", styles["SectionHead"]))
    story.append(divider(dw))

    exec_summary = [
        "This report presents the complete results of the <b>astronaut-opposite-forcing</b> pipeline, "
        "a reproducible workflow for identifying compounds and nutrients that counteract spaceflight-induced "
        "transcriptomic changes. The pipeline integrates NASA Open Science Data Repository (OSDR) rodent "
        "spaceflight RNA-seq data, LINCS L1000 connectivity mapping, multi-database drug annotation, "
        "functional enrichment, and vegan nutritional translation.",

        "Two tissue-stratified consensus signatures were built from 15 rodent ISS studies using "
        "random-effects meta-analysis (metafor REML) with a consensus definition of nominal p < 0.05 "
        "and >=70% direction consistency across studies. The <b>liver/immune signature</b> (525 genes) "
        "captures cell-cycle suppression (E2F targets, G2M checkpoint downregulated) and metabolic "
        "reprogramming (fatty acid metabolism, oxidative phosphorylation, xenobiotic metabolism upregulated). "
        "The <b>brain signature</b> (505 genes) reveals mitochondrial dysfunction with oxidative "
        "phosphorylation and complex I assembly downregulated — a hallmark of spaceflight biology.",

        "LINCS L1000 connectivity mapping (local GMT fallback due to SigCom API unavailability) identified "
        "219 liver and 221 brain candidate reversers. Top liver reversers include <b>Imatinib</b> "
        "(Bcr-Abl kinase inhibitor), <b>Metformin</b> (insulin sensitizer/AMPK activator), and "
        "<b>Quercetin</b> (flavonoid antioxidant). Top brain reversers include <b>Vemurafenib</b> "
        "(RAF inhibitor), <b>Formoterol</b> (beta-2 adrenergic agonist), and <b>Sulforaphane</b> "
        "(Nrf2 activator).",

        "Functional enrichment confirmed biological plausibility: GSEA Hallmark identified 32 significant "
        "liver and 23 brain pathways (FDR < 0.25). The vegan nutrition layer mapped 10 countermeasure "
        "nutrients to FooDB food sources and curated 12 phase-optimized vegan recipes covering pre-flight, "
        "in-flight, and post-flight mission phases.",
    ]
    for para in exec_summary:
        story.append(Paragraph(para, styles["Body"]))

    story.append(Spacer(1, 10))
    story.append(callout_box(
        "<b>Key finding:</b> Spaceflight induces tissue-specific transcriptomic signatures with "
        "distinct metabolic and mitochondrial vulnerabilities. Pharmacological reversers (Imatinib, "
        "Metformin, Formoterol) and vegan nutritional countermeasures (sulforaphane, curcumin, "
        "quercetin, EGCG) target complementary pathways — suggesting a combined pharmacological + "
        "nutritional countermeasure strategy for long-duration spaceflight.",
        width=dw
    ))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # METHODS
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("2. Methods", styles["SectionHead"]))
    story.append(divider(dw))

    # 2.1 Data acquisition
    story.append(Paragraph("2.1 Data Acquisition", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Rodent spaceflight RNA-seq studies were identified via the NASA OSDR BDAPI "
        "(https://osdr.nasa.gov/bio/api/). The study catalog was filtered for Mus musculus ISS missions "
        "with RNA-seq count matrices available. 160 studies were screened; 31 valid STAR count matrices "
        "were downloaded. Human cross-check studies were also cataloged for validation.",
        styles["Body"]))

    # 2.2 Differential expression
    story.append(Paragraph("2.2 Differential Expression Analysis", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Per-study differential expression was computed with DESeq2 (design: ~condition, Ground Control "
        "as reference). Sample conditions were parsed from sample_metadata.json (factor value: spaceflight) "
        "with fallback parsing from sample names (FLT=flight, GC=ground, VIV=vivarium). ERCC spike-in "
        "controls were removed before meta-analysis. apeglm LFC shrinkage was applied; studies with "
        "over-shrunk LFCs (median SE < 0.01) were excluded, retaining 15 studies with valid DE signal.",
        styles["Body"]))

    # 2.3 Meta-analysis
    story.append(Paragraph("2.3 Random-Effects Meta-Analysis", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Gene-level random-effects meta-analysis was performed with metafor::rma (REML method, "
        "DerSimonian-Laird fallback). Mouse-to-human ortholog mapping used org.Mm.eg.db + homologene "
        "(biomaRt was unavailable). Consensus genes were defined as nominal p < 0.05 AND >=70% direction "
        "consistency across studies. Core genes additionally required |pooled_log2FC| >= 0.5. "
        "Technical gene filters removed ribosomal (Rpl/Rps), housekeeping (GAPDH, ACTB, B2M), and "
        "translation factors (Eef/Eif); metallothioneins (MT1/MT2) were retained.",
        styles["Body"]))

    # 2.4 LINCS screening
    story.append(Paragraph("2.4 LINCS L1000 Connectivity Mapping", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Consensus signatures (human HGNC symbols) were screened against the LINCS L1000 perturbation "
        "compendium via the SigCom API. Gene resolution succeeded at 95-99% coverage. The SigCom data API "
        "was unavailable (503) during analysis; a local GMT fallback (hypergeometric overlap against "
        "single_drug_perturbations-v1.0.gmt) was used. Compounds were annotated with the Broad Drug "
        "Repurposing Hub (clinical phase, MoA, target, indication). Compound name artifacts were filtered "
        "by regex validation.",
        styles["Body"]))

    # 2.5 Enrichment
    story.append(Paragraph("2.5 Functional Enrichment", styles["SubSectionHead"]))
    story.append(Paragraph(
        "GSEA was performed with clusterProfiler::GSEA on Hallmark and Reactome (msigdbr C2/CP:REACTOME) "
        "gene sets. ORA was performed with enricher() on GO:BP and Reactome. KEGG was excluded due to "
        "commercial license restrictions. Disease-gene enrichment used JensenLab DISEASES (knowledge and "
        "experiments channels). Significance thresholds: GSEA p.adjust < 0.25 (standard for exploratory "
        "GSEA), ORA p.adjust < 0.05.",
        styles["Body"]))

    # 2.6 Nutrition
    story.append(Paragraph("2.6 Nutrient-Gene Mapping and Vegan Recipe Optimization", styles["SubSectionHead"]))
    story.append(Paragraph(
        "LINCS reverser compounds were mapped to FooDB (foodb.ca, 2020-04-07 release) by compound name "
        "matching (exact then partial first-word match >= 4 chars). Vegan food sources were classified "
        "by keyword filtering. 10 countermeasure nutrients were defined based on LINCS hits, FooDB "
        "coverage, and known bioactivity. 12 curated vegan recipes were scored by nutrient coverage, "
        "LINCS validation, and mission-phase relevance, organized into pre-flight, in-flight, and "
        "post-flight meal plans.",
        styles["Body"]))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("3. Results", styles["SectionHead"]))
    story.append(divider(dw))

    # 3.1 Consensus signatures
    story.append(Paragraph("3.1 Consensus Spaceflight Signatures", styles["SubSectionHead"]))
    story.append(Paragraph(
        f"Random-effects meta-analysis of 15 rodent ISS studies yielded two tissue-stratified consensus "
        f"signatures. The <b>liver/immune signature</b> comprises {liver_up + liver_dn} consensus genes "
        f"({liver_up} upregulated, {liver_dn} downregulated) with {liver_core} core genes "
        f"(|log2FC| >= 0.5). The <b>brain signature</b> comprises {brain_up + brain_dn} consensus genes "
        f"({brain_up} upregulated, {brain_dn} downregulated) with {brain_core} core genes. "
        f"Ortholog mapping achieved 89% coverage for consensus genes (39.1% genome-wide).",
        styles["Body"]))

    story.append(embed_figure("fig1_signature_overview.png",
        "Figure 1. Consensus spaceflight transcriptomic signatures by tissue stratum. "
        "Hatched bars indicate core genes with |pooled log2FC| >= 0.5.", dw, max_height=280))

    # Top consensus genes table
    story.append(Paragraph("Top core consensus genes (liver/immune):", styles["BodyLeft"]))
    liver_core_genes = data["liver_sig"][data["liver_sig"]["core"] == True].sort_values(
        "pooled_log2FC", key=lambda x: x.abs(), ascending=False).head(10)
    liver_rows = []
    for _, r in liver_core_genes.iterrows():
        sym = r["human_symbol"] if pd.notna(r["human_symbol"]) else r["mouse_symbol"]
        liver_rows.append([sym, f"{r['pooled_log2FC']:.2f}", f"{r['pvalue']:.2e}",
                          r["direction"], str(r["k"])])
    story.append(make_table(
        ["Gene", "log2FC", "P-value", "Direction", "N studies"],
        liver_rows, [110, 70, 80, 70, 70]))
    story.append(Spacer(1, 8))

    story.append(Paragraph("Top core consensus genes (brain):", styles["BodyLeft"]))
    brain_core_genes = data["brain_sig"][data["brain_sig"]["core"] == True].sort_values(
        "pooled_log2FC", key=lambda x: x.abs(), ascending=False).head(10)
    brain_rows = []
    for _, r in brain_core_genes.iterrows():
        sym = r["human_symbol"] if pd.notna(r["human_symbol"]) else r["mouse_symbol"]
        brain_rows.append([sym, f"{r['pooled_log2FC']:.2f}", f"{r['pvalue']:.2e}",
                          r["direction"], str(r["k"])])
    story.append(make_table(
        ["Gene", "log2FC", "P-value", "Direction", "N studies"],
        brain_rows, [110, 70, 80, 70, 70]))

    story.append(PageBreak())

    # 3.2 LINCS reversers
    story.append(Paragraph("3.2 LINCS L1000 Signature Reversers", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Connectivity mapping identified compounds whose perturbation signatures reverse the spaceflight "
        "consensus signatures (negative z-sum = reversal). The local GMT fallback provides coarser "
        "results than the full SigCom API (no cell-line context or tiered scoring) but identifies real "
        "compounds with validated mechanisms.",
        styles["Body"]))

    story.append(embed_figure("fig2_top_reversers.png",
        "Figure 2. Top 10 LINCS L1000 signature reversers for liver/immune (left) and brain (right) "
        "signatures. More negative z-sum indicates stronger reversal.", dw, max_height=300))

    # Top reversers tables
    story.append(Paragraph("Top liver/immune reversers:", styles["BodyLeft"]))
    liver_lincs = data["liver_lincs"]
    liver_lincs_valid = liver_lincs[liver_lincs["compound"].apply(
        lambda x: len(str(x)) > 2 and not str(x).isdigit())].head(8)
    liver_lincs_rows = []
    for _, r in liver_lincs_valid.iterrows():
        phase = str(r.get("drh_clinical_phase", "")) if pd.notna(r.get("drh_clinical_phase")) else ""
        moa = str(r.get("drh_moa", "")) if pd.notna(r.get("drh_moa")) else ""
        moa_short = moa[:30] + "..." if len(moa) > 30 else moa
        liver_lincs_rows.append([
            str(r["compound"]), f"{r['best_z_sum']:.2f}",
            phase, moa_short
        ])
    story.append(make_table(
        ["Compound", "Best Z-Sum", "Clinical Phase", "Mechanism of Action"],
        liver_lincs_rows, [90, 70, 70, 200], cell_style=styles["TableCellSmall"]))
    story.append(Spacer(1, 8))

    story.append(Paragraph("Top brain reversers:", styles["BodyLeft"]))
    brain_lincs = data["brain_lincs"]
    brain_lincs_valid = brain_lincs[brain_lincs["compound"].apply(
        lambda x: len(str(x)) > 2 and not str(x).isdigit())].head(8)
    brain_lincs_rows = []
    for _, r in brain_lincs_valid.iterrows():
        phase = str(r.get("drh_clinical_phase", "")) if pd.notna(r.get("drh_clinical_phase")) else ""
        moa = str(r.get("drh_moa", "")) if pd.notna(r.get("drh_moa")) else ""
        moa_short = moa[:30] + "..." if len(moa) > 30 else moa
        brain_lincs_rows.append([
            str(r["compound"]), f"{r['best_z_sum']:.2f}",
            phase, moa_short
        ])
    story.append(make_table(
        ["Compound", "Best Z-Sum", "Clinical Phase", "Mechanism of Action"],
        brain_lincs_rows, [90, 70, 70, 200], cell_style=styles["TableCellSmall"]))

    story.append(PageBreak())

    # 3.3 Functional enrichment
    story.append(Paragraph("3.3 Functional Enrichment", styles["SubSectionHead"]))
    story.append(Paragraph(
        "GSEA Hallmark identified 32 significant liver and 23 significant brain pathways (p.adjust < 0.25). "
        "The liver/immune signature is characterized by downregulated cell-cycle programs (E2F targets, "
        "G2M checkpoint, MYC targets) and upregulated metabolic pathways (fatty acid metabolism, oxidative "
        "phosphorylation, xenobiotic metabolism, bile acid metabolism). The brain signature shows "
        "downregulated oxidative phosphorylation and mitochondrial complex I assembly — consistent with "
        "the mitochondrial dysfunction widely reported in spaceflight biology.",
        styles["Body"]))

    story.append(embed_figure("fig3_gsea_hallmark.png",
        "Figure 3. GSEA Hallmark pathway enrichment. Dot size represents -log10(p.adjust). "
        "Negative NES = downregulated in spaceflight; positive NES = upregulated.", dw, max_height=350))

    # Top GSEA pathways table
    story.append(Paragraph("Top GSEA Hallmark pathways (liver/immune):", styles["BodyLeft"]))
    liver_gsea = data["liver_gsea_h"].sort_values("p.adjust").head(8)
    liver_gsea_rows = []
    for _, r in liver_gsea.iterrows():
        pw = r["Description"].replace("HALLMARK_", "").replace("_", " ")
        liver_gsea_rows.append([pw, f"{r['NES']:.2f}", f"{r['p.adjust']:.2e}", str(r["setSize"])])
    story.append(make_table(
        ["Pathway", "NES", "p.adjust", "Set Size"],
        liver_gsea_rows, [200, 60, 80, 60], cell_style=styles["TableCellSmall"]))
    story.append(Spacer(1, 8))

    story.append(Paragraph("Top GSEA Hallmark pathways (brain):", styles["BodyLeft"]))
    brain_gsea = data["brain_gsea_h"].sort_values("p.adjust").head(8)
    brain_gsea_rows = []
    for _, r in brain_gsea.iterrows():
        pw = r["Description"].replace("HALLMARK_", "").replace("_", " ")
        brain_gsea_rows.append([pw, f"{r['NES']:.2f}", f"{r['p.adjust']:.2e}", str(r["setSize"])])
    story.append(make_table(
        ["Pathway", "NES", "p.adjust", "Set Size"],
        brain_gsea_rows, [200, 60, 80, 60], cell_style=styles["TableCellSmall"]))

    story.append(PageBreak())

    # 3.4 Nutrient-gene mapping
    story.append(Paragraph("3.4 Nutrient-Gene Mapping", styles["SubSectionHead"]))
    story.append(Paragraph(
        "LINCS reverser compounds were mapped to FooDB to identify vegan food sources. Of the top 50 "
        "compounds per stratum, 20 liver and 11 brain compounds matched FooDB entries. Known therapeutic "
        "nutrients with vegan food sources include Quercetin (212 vegan foods), Sulforaphane (3 foods), "
        "Ascorbic acid (327 foods), Luteolin (292 foods), and Curcumin (287 foods). Food contaminants "
        "(Cadmium, Arsenic) were correctly flagged as non-therapeutic.",
        styles["Body"]))

    # Nutrient map table
    nutrient = data["nutrient_map"]
    nutrient_with_food = nutrient[(nutrient["foodb_matched"] == True) & (nutrient["is_known_nutrient"] == True)].drop_duplicates("compound")
    nutrient_rows = []
    for _, r in nutrient_with_food.iterrows():
        foods = str(r.get("vegan_food_sources", ""))[:60]
        n_foods = str(r.get("n_vegan_foods", 0))
        nutrient_rows.append([str(r["compound"]), str(r["stratum"]), n_foods, foods])
    story.append(make_table(
        ["Compound", "Stratum", "N Vegan Foods", "Example Food Sources"],
        nutrient_rows, [80, 80, 60, 230], cell_style=styles["TableCellSmall"]))

    story.append(PageBreak())

    # 3.5 Vegan recipes
    story.append(Paragraph("3.5 Vegan Recipe Optimization", styles["SubSectionHead"]))
    story.append(Paragraph(
        "Twelve curated vegan recipes were scored by nutrient coverage, LINCS validation, and mission-phase "
        "relevance. Recipes are organized into three mission phases: pre-flight (senolytic and antioxidant "
        "preparation), in-flight (mitochondrial support and anti-inflammatory), and post-flight (recovery "
        "and reconditioning). The top-scoring recipe is the Broccoli Sulforaphane Power Bowl (score=143), "
        "targeting sulforaphane, curcumin, ascorbic acid, quercetin, and lutein.",
        styles["Body"]))

    story.append(embed_figure("fig4_nutrient_coverage.png",
        "Figure 4. Vegan recipe nutrient coverage heatmap. Green cells indicate the nutrient is "
        "targeted by the recipe. Mission phase shown at left.", dw, max_height=320))

    story.append(embed_figure("fig5_recipe_scores.png",
        "Figure 5. Recipe scores by mission phase. Higher scores indicate broader countermeasure "
        "nutrient coverage.", dw, max_height=300))

    # Recipe table
    recipes = data["recipes"].sort_values("nutrient_score", ascending=False).head(6)
    recipe_rows = []
    for _, r in recipes.iterrows():
        nutrients = str(r["target_nutrients"])[:50]
        recipe_rows.append([str(r["name"]), str(r["phase"]).replace("_", " "), str(int(r["nutrient_score"])), nutrients])
    story.append(make_table(
        ["Recipe", "Phase", "Score", "Target Nutrients"],
        recipe_rows, [130, 70, 50, 200], cell_style=styles["TableCellSmall"]))

    story.append(PageBreak())

    # Meal plan table
    story.append(Paragraph("Mission-Phase Meal Plans", styles["SubSectionHead"]))
    meal_plans = data["meal_plans"]
    for phase in ["pre_flight", "in_flight", "post_flight"]:
        phase_df = meal_plans[meal_plans["phase"] == phase]
        if len(phase_df) == 0:
            continue
        phase_label = phase.replace("_", " ").title()
        story.append(Paragraph(f"<b>{phase_label}</b>", styles["BodyLeft"]))
        mp_rows = []
        for _, r in phase_df.iterrows():
            mp_rows.append([str(r["meal_type"]).title(), str(r["recipe_name"]), str(r["target_nutrients"])[:60]])
        story.append(make_table(
            ["Meal", "Recipe", "Target Nutrients"],
            mp_rows, [70, 160, 220], cell_style=styles["TableCellSmall"]))
        story.append(Spacer(1, 8))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # DISCUSSION
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("4. Discussion", styles["SectionHead"]))
    story.append(divider(dw))

    discussion = [
        "The consensus signatures recapitulate known spaceflight biology with notable tissue specificity. "
        "The liver/immune signature's downregulation of E2F targets and G2M checkpoint genes reflects the "
        "cell-cycle suppression observed across multiple spaceflight studies — likely a consequence of "
        "altered mechanical loading, radiation exposure, and stress signaling. The concurrent upregulation "
        "of fatty acid metabolism, oxidative phosphorylation, and xenobiotic metabolism suggests a hepatic "
        "metabolic shift toward energy storage and detoxification, potentially compensating for altered "
        "nutrient availability and xenobiotic load in the spacecraft environment.",

        "The brain signature's downregulation of oxidative phosphorylation and mitochondrial complex I "
        "assembly is particularly significant. Mitochondrial dysfunction is one of the most consistently "
        "reported molecular effects of spaceflight, observed across rodent studies, NASA Twin Study data, "
        "and in vitro experiments. The downregulation of complex I (NADH dehydrogenase) specifically "
        "implicates impaired electron transport chain function, which could contribute to the cognitive "
        "changes and neuroinflammation reported in astronauts.",

        "The LINCS reversers span diverse mechanisms. <b>Imatinib</b> (top liver reverser) is a multi-kinase "
        "inhibitor targeting Bcr-Abl, KIT, and PDGFR — its reversal of the liver signature may relate to "
        "PDGFR signaling in hepatic stellate cells and immune modulation. <b>Metformin</b> activates AMPK, "
        "which counteracts the metabolic reprogramming seen in the liver signature and has documented "
        "mitochondrial protective effects. <b>Formoterol</b> (top brain reverser) is a beta-2 adrenergic "
        "agonist that may counteract muscle atrophy signaling and has neuroprotective properties. "
        "<b>Sulforaphane</b> activates Nrf2 antioxidant response, directly addressing the oxidative stress "
        "component of spaceflight mitochondrial dysfunction.",

        "The vegan nutritional countermeasures provide an accessible, low-risk adjunct to pharmacological "
        "intervention. Sulforaphane (from cruciferous vegetables), curcumin (from turmeric), quercetin "
        "(from berries and alliums), and EGCG (from green tea) all have documented anti-inflammatory, "
        "antioxidant, and mitochondrial protective effects. The curated recipes ensure practical, "
        "palatable delivery of these bioactives within mission-appropriate meal formats. The phase-based "
        "design (pre-flight senolytic preparation, in-flight mitochondrial support, post-flight recovery) "
        "aligns nutritional intervention with the temporal dynamics of spaceflight adaptation.",
    ]
    for para in discussion:
        story.append(Paragraph(para, styles["Body"]))

    story.append(Spacer(1, 10))

    # Limitations
    story.append(Paragraph("Limitations", styles["SubSectionHead"]))
    limitations = [
        "The LINCS L1000 screening used a local GMT fallback (hypergeometric overlap) rather than the "
        "full SigCom API, which was unavailable during analysis. The fallback provides coarser results "
        "without cell-line context or tiered scoring. Re-running with the full API when available would "
        "improve resolution.",
        "Meta-analysis heterogeneity (I-squared ~99%) is high, reflecting pooling across different tissues, "
        "missions, and durations. BH-FDR across ~44k genes with this heterogeneity yields 0 significant "
        "genes — the consensus definition (nominal p < 0.05 + >=70% direction consistency) is standard "
        "for signature building but less stringent than genome-wide FDR control.",
        "Ortholog mapping achieved 89% coverage for consensus genes but only 39.1% genome-wide, "
        "introducing potential bias toward well-conserved genes.",
        "JensenLab DISEASES enrichment yielded limited results because spaceflight is not a disease "
        "ontology term. Functional enrichment (GSEA/ORA) provides the stronger biological validation.",
        "RecipeNLG (2.2M recipes) was not programmatically accessible; curated recipes were used instead. "
        "While more targeted, this limits recipe diversity.",
        "FooDB matching by compound name may miss compounds with alternative nomenclature. InChIKey-based "
        "matching would improve coverage but requires additional identifier resolution.",
    ]
    for lim in limitations:
        story.append(Paragraph(f"• {lim}", styles["BulletText"]))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # CONCLUSIONS
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("5. Conclusions and Next Steps", styles["SectionHead"]))
    story.append(divider(dw))

    story.append(Paragraph(
        "This pipeline demonstrates a reproducible computational approach to identifying countermeasures "
        "for spaceflight-induced biological changes. The consensus signatures capture tissue-specific "
        "vulnerabilities — hepatic metabolic reprogramming and cerebral mitochondrial dysfunction — that "
        "are consistent with the broader spaceflight biology literature. The combination of pharmacological "
        "reversers (Imatinib, Metformin, Formoterol) and vegan nutritional countermeasures (sulforaphane, "
        "curcumin, quercetin, EGCG) suggests a multi-modal countermeasure strategy.",
        styles["Body"]))

    story.append(Paragraph("Recommended Next Steps", styles["SubSectionHead"]))
    next_steps = [
        "<b>Re-run LINCS screening with SigCom API</b> when the service recovers, to obtain cell-line-resolved, "
        "tiered connectivity scores and validate the GMT fallback results.",
        "<b>Validate top reversers in vitro</b> using spaceflight-analog models (clinostat, random "
        "positioning machine, or radiation exposure) with transcriptomic readout.",
        "<b>Expand the human cross-check</b> by meta-analyzing human spaceflight transcriptomic data "
        "(NASA Twin Study, JAXA studies) against the rodent-derived signatures.",
        "<b>Refine recipe optimization</b> by integrating RecipeNLG (when accessible) for broader recipe "
        "diversity, and add nutrient bioavailability and shelf-stability constraints for spaceflight logistics.",
        "<b>Conduct dose-response modeling</b> for the top nutritional countermeasures, estimating "
        "therapeutic intake levels achievable through diet alone vs. supplementation.",
        "<b>Integrate proteomics and metabolomics</b> data from OSDR to validate transcriptomic signatures "
        "at the protein and metabolite level.",
    ]
    for ns in next_steps:
        story.append(Paragraph(f"• {ns}", styles["BulletText"]))

    story.append(PageBreak())

    # ═══════════════════════════════════════════════════════
    # REFERENCES
    # ═══════════════════════════════════════════════════════
    story.append(Paragraph("6. References", styles["SectionHead"]))
    story.append(divider(dw))

    refs = [
        "NASA Open Science Data Repository (OSDR). https://osdr.nasa.gov. Accessed 2025.",
        "LINCS L1000 / Connectivity Map. Broad Institute. https://clue.io. SigCom API: https://maayanlab.cloud/sigcom-lincs.",
        "Subramanian A, et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. Cell. 2017;171(6):1437-1452.",
        "Dewey FE, et al. Druggable genome and systems biology. Nat Rev Drug Discov. 2011.",
        "Wishart DS, et al. FooDB: a food compound database. https://foodb.ca. 2020 release.",
        "Jensen LJ, et al. JensenLab DISEASES: text mining and data integration for disease-gene associations. Database. 2016.",
        "Liberzon A, et al. The Molecular Signatures Database (MSigDB) Hallmark Gene Set Collection. Cell Systems. 2015;1(6):417-425.",
        "Yu G, et al. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287.",
        "Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48.",
        "Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15:550.",
        "Zhu Q, et al. apeglm: approximate posterior estimation for RNA-seq log fold changes. Nat Methods. 2019.",
        "DrugBank 5.1.8. Wishart DS, et al. DrugBank 5.0: a major update. Nucleic Acids Res. 2018. Open-access XML, CC BY-NC 4.0.",
        "Broad Drug Repurposing Hub. https://clue.io/repurposing-app.",
        "ChEMBL. Gaulton A, et al. The ChEMBL database. Nucleic Acids Res. 2012.",
        "Beheshti A, et al. Multi-omics and systems biology approaches for understanding spaceflight effects. Cell. 2019.",
        "Garrett-Bakelman FE, et al. The NASA Twins Study: A multidimensional analysis of a year-long human spaceflight. Science. 2019;364(6436).",
    ]
    for i, ref in enumerate(refs, 1):
        story.append(Paragraph(f"[{i}] {ref}", styles["RefText"]))

    story.append(Spacer(1, 20))
    story.append(divider(dw))
    story.append(Paragraph(
        "<i>This report was generated by Biomni (Phylo). The pipeline code, data, and full results "
        "are available in the Zenodo package: astronaut-opposite-forcing.</i>",
        styles["Attribution"]))

    # Build
    print(f"Building PDF: {OUTPUT_PATH}")
    doc.build(story, onFirstPage=page_header_footer, onLaterPages=page_header_footer)
    print(f"PDF generated: {os.path.getsize(OUTPUT_PATH)} bytes")


# ──────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────
def validate_pdf():
    from pypdf import PdfReader
    reader = PdfReader(OUTPUT_PATH)
    page_count = len(reader.pages)
    file_size = os.path.getsize(OUTPUT_PATH)
    print(f"\nValidation:")
    print(f"  Pages: {page_count}")
    print(f"  File size: {file_size} bytes")
    assert page_count >= 5, f"Report has only {page_count} pages"
    assert file_size > 50000, f"Report is only {file_size} bytes"
    first_page_text = reader.pages[0].extract_text()
    assert len(first_page_text.strip()) > 0, "First page has no extractable text"
    print("  All checks passed.")


if __name__ == "__main__":
    print("=" * 60)
    print("Script 07b: PDF Report Generation")
    print("=" * 60)
    build_report()
    validate_pdf()
    print("\nDone.")
