#!/usr/bin/env python3
"""
Script 07a: Generate figures for the PDF report and manuscript.

Produces SVG + PNG figures:
  fig1_signature_overview  - up/down gene counts per stratum
  fig2_top_reversers        - top LINCS reversers (liver + brain)
  fig3_gsea_hallmark        - GSEA Hallmark NES dotplot (liver + brain)
  fig4_nutrient_coverage    - vegan recipe nutrient coverage heatmap
  fig5_recipe_scores        - recipe scores by phase

All figures use Phylo brand colors, Liberation Sans font, colorblind-friendly.
"""

import os
import json
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch

# ── Phylo brand colors ──
PHYLO_GOLD   = "#D4A04A"
PHYLO_BLUE   = "#0279EE"
PHYLO_GREEN  = "#75A025"
PHYLO_ORANGE = "#FF9400"
PHYLO_PINK   = "#FD9BED"
PHYLO_BLACK  = "#000000"
PHYLO_WARM   = "#ECE9E2"
PHYLO_OFF    = "#FAF9F3"

CHART_COLORS = [PHYLO_GOLD, PHYLO_BLUE, PHYLO_GREEN, PHYLO_ORANGE, PHYLO_PINK, PHYLO_BLACK]

# ── Matplotlib config ──
rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
rcParams['svg.fonttype'] = 'none'  # keep SVG text editable
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False

PROJ = "/workspace/astronaut-opposite-forcing"
FIG_DIR = os.path.join(PROJ, "results", "figures")
os.makedirs(FIG_DIR, exist_ok=True)


def save_fig(fig, name):
    """Save figure as both SVG and PNG."""
    svg_path = os.path.join(FIG_DIR, f"{name}.svg")
    png_path = os.path.join(FIG_DIR, f"{name}.png")
    fig.savefig(svg_path, format="svg", bbox_inches="tight", dpi=150)
    fig.savefig(png_path, format="png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  Saved: {name}.svg + {name}.png")


# ──────────────────────────────────────────────────────────
# Figure 1: Signature Overview — up/down gene counts per stratum
# ──────────────────────────────────────────────────────────
def fig1_signature_overview():
    print("Generating fig1_signature_overview...")
    liver = pd.read_csv(os.path.join(PROJ, "data/processed/consensus_signature_liver_immune.csv"))
    brain = pd.read_csv(os.path.join(PROJ, "data/processed/consensus_signature_brain.csv"))

    strata = ["Liver / Immune", "Brain"]
    up_counts = [
        (liver["direction"] == "up").sum(),
        (brain["direction"] == "up").sum(),
    ]
    dn_counts = [
        (liver["direction"] == "down").sum(),
        (brain["direction"] == "down").sum(),
    ]
    core_up = [
        ((liver["direction"] == "up") & (liver["core"] == True)).sum(),
        ((brain["direction"] == "up") & (brain["core"] == True)).sum(),
    ]
    core_dn = [
        ((liver["direction"] == "down") & (liver["core"] == True)).sum(),
        ((brain["direction"] == "down") & (brain["core"] == True)).sum(),
    ]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    x = np.arange(len(strata))
    w = 0.32

    bars1 = ax.bar(x - w/2, up_counts, w, label="Upregulated", color=PHYLO_ORANGE, edgecolor="white", linewidth=0.5)
    bars2 = ax.bar(x + w/2, dn_counts, w, label="Downregulated", color=PHYLO_BLUE, edgecolor="white", linewidth=0.5)

    # Core gene markers (hatched overlay)
    ax.bar(x - w/2, core_up, w, color=PHYLO_ORANGE, edgecolor=PHYLO_BLACK, linewidth=0.8, hatch="//", alpha=0.0)
    ax.bar(x - w/2, core_up, w, color="none", edgecolor=PHYLO_BLACK, linewidth=0.8, hatch="//")

    for bar, val in zip(bars1, up_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, str(val),
                ha="center", va="bottom", fontsize=11, fontweight="bold")
    for bar, val in zip(bars2, dn_counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, str(val),
                ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax.set_ylabel("Consensus Gene Count", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(strata, fontsize=12)
    ax.set_title("Spaceflight Consensus Transcriptomic Signatures", fontsize=14, fontweight="bold", pad=12)
    ax.legend(loc="upper right", fontsize=10, frameon=False)

    # Annotation for core genes
    ax.text(0.02, 0.95, "Hatched bars = core genes (|log2FC| >= 0.5)",
            transform=ax.transAxes, fontsize=8, color="#8A8378", va="top")

    ax.set_ylim(0, max(max(up_counts), max(dn_counts)) * 1.25)
    save_fig(fig, "fig1_signature_overview")


# ──────────────────────────────────────────────────────────
# Figure 2: Top LINCS Reversers (liver + brain side by side)
# ──────────────────────────────────────────────────────────
def fig2_top_reversers():
    print("Generating fig2_top_reversers...")
    liver = pd.read_csv(os.path.join(PROJ, "results/tables/lincs_liver_immune/tier1_ranking_annotated.csv"))
    brain = pd.read_csv(os.path.join(PROJ, "results/tables/lincs_brain/tier1_ranking_annotated.csv"))

    # Filter valid compounds, take top 10 by best_z_sum (most negative = best reversal)
    liver_valid = liver[liver["compound"].apply(lambda x: len(x) > 2 and not x.isdigit())].head(10)
    brain_valid = brain[brain["compound"].apply(lambda x: len(x) > 2 and not x.isdigit())].head(10)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

    for ax, df, title, color in [
        (ax1, liver_valid, "Liver / Immune Signature", PHYLO_GOLD),
        (ax2, brain_valid, "Brain Signature", PHYLO_BLUE),
    ]:
        compounds = df["compound"].values[::-1]
        z_sums = df["best_z_sum"].values[::-1]
        # Use clinical phase for color intensity
        phases = df.get("drh_clinical_phase", pd.Series([""] * len(df))).values[::-1]

        bars = ax.barh(range(len(compounds)), z_sums, color=color, edgecolor="white", linewidth=0.5, height=0.7)
        ax.set_yticks(range(len(compounds)))
        ax.set_yticklabels(compounds, fontsize=10)
        ax.set_xlabel("Best Z-Sum (negative = reversal)", fontsize=11)
        ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
        ax.axvline(x=0, color=PHYLO_BLACK, linewidth=0.5)

        # Annotate clinical phase
        for i, (bar, phase) in enumerate(zip(bars, phases)):
            if isinstance(phase, str) and phase:
                ax.text(bar.get_width() - 0.3, bar.get_y() + bar.get_height()/2,
                        phase[:8], ha="right", va="center", fontsize=7, color="white", fontweight="bold")

    fig.suptitle("Top LINCS L1000 Signature Reversers", fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    save_fig(fig, "fig2_top_reversers")


# ──────────────────────────────────────────────────────────
# Figure 3: GSEA Hallmark NES dotplot (liver + brain)
# ──────────────────────────────────────────────────────────
def fig3_gsea_hallmark():
    print("Generating fig3_gsea_hallmark...")
    liver = pd.read_csv(os.path.join(PROJ, "results/tables/gsea_liver_immune_hallmark.csv"))
    brain = pd.read_csv(os.path.join(PROJ, "results/tables/gsea_brain_hallmark.csv"))

    # Filter significant (p.adjust < 0.25 per GSEA convention) and take top 12 by |NES|
    liver_sig = liver[liver["p.adjust"] < 0.25].copy()
    brain_sig = brain[brain["p.adjust"] < 0.25].copy()
    liver_sig["abs_NES"] = liver_sig["NES"].abs()
    brain_sig["abs_NES"] = brain_sig["NES"].abs()
    liver_top = liver_sig.nlargest(12, "abs_NES")
    brain_top = brain_sig.nlargest(12, "abs_NES")

    # Combine for a shared dotplot
    all_pathways = list(set(liver_top["Description"].tolist() + brain_top["Description"].tolist()))
    all_pathways.sort()

    rows = []
    for pw in all_pathways:
        l_row = liver_sig[liver_sig["Description"] == pw]
        b_row = brain_sig[brain_sig["Description"] == pw]
        l_nes = l_row["NES"].values[0] if len(l_row) > 0 else 0
        b_nes = b_row["NES"].values[0] if len(b_row) > 0 else 0
        l_padj = l_row["p.adjust"].values[0] if len(l_row) > 0 else 1
        b_padj = b_row["p.adjust"].values[0] if len(b_row) > 0 else 1
        rows.append({"pathway": pw, "liver_NES": l_nes, "brain_NES": b_nes,
                      "liver_padj": l_padj, "brain_padj": b_padj})
    df = pd.DataFrame(rows).sort_values("liver_NES", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(6, len(df) * 0.35)))
    y = np.arange(len(df))
    w = 0.35

    # Size by -log10(p.adjust)
    l_sizes = -np.log10(df["liver_padj"].clip(lower=1e-10)) * 15 + 20
    b_sizes = -np.log10(df["brain_padj"].clip(lower=1e-10)) * 15 + 20

    ax.scatter(df["liver_NES"], y - w/2, s=l_sizes, c=PHYLO_GOLD, alpha=0.8, edgecolors="white", linewidth=0.5, label="Liver / Immune")
    ax.scatter(df["brain_NES"], y + w/2, s=b_sizes, c=PHYLO_BLUE, alpha=0.8, edgecolors="white", linewidth=0.5, label="Brain")

    ax.set_yticks(y)
    labels = [p.replace("HALLMARK_", "").replace("_", " ").title() for p in df["pathway"]]
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Normalized Enrichment Score (NES)", fontsize=11)
    ax.set_title("GSEA Hallmark Pathways — Spaceflight Signatures", fontsize=14, fontweight="bold", pad=12)
    ax.axvline(x=0, color=PHYLO_BLACK, linewidth=0.5, linestyle="--")
    ax.legend(loc="lower right", fontsize=10, frameon=False)

    # Size legend
    ax.text(0.02, 0.02, "Dot size = -log10(p.adjust)", transform=ax.transAxes,
            fontsize=8, color="#8A8378", va="bottom")

    save_fig(fig, "fig3_gsea_hallmark")


# ──────────────────────────────────────────────────────────
# Figure 4: Vegan Recipe Nutrient Coverage Heatmap
# ──────────────────────────────────────────────────────────
def fig4_nutrient_coverage():
    print("Generating fig4_nutrient_coverage...")
    recipes = pd.read_csv(os.path.join(PROJ, "results/tables/vegan_recipes_scored.csv"))

    # Parse target nutrients
    all_nutrients = set()
    for tn in recipes["target_nutrients"]:
        for n in str(tn).split(";"):
            n = n.strip()
            if n:
                all_nutrients.add(n)
    all_nutrients = sorted(all_nutrients)

    # Build coverage matrix
    matrix = []
    recipe_names = []
    for _, row in recipes.iterrows():
        recipe_names.append(row["name"])
        nutrients = [n.strip() for n in str(row["target_nutrients"]).split(";")]
        matrix.append([1 if n in nutrients else 0 for n in all_nutrients])

    mat = np.array(matrix)

    fig, ax = plt.subplots(figsize=(10, max(5, len(recipe_names) * 0.45)))
    cmap = matplotlib.colors.ListedColormap([PHYLO_OFF, PHYLO_GREEN])
    im = ax.imshow(mat, aspect="auto", cmap=cmap, interpolation="nearest")

    ax.set_xticks(range(len(all_nutrients)))
    ax.set_xticklabels([n.title() for n in all_nutrients], rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(len(recipe_names)))
    ax.set_yticklabels(recipe_names, fontsize=9)

    # Add phase annotation on y-axis
    phases = recipes["phase"].values
    for i, ph in enumerate(phases):
        ax.text(-0.5, i, ph.replace("_", " ").title(), ha="right", va="center", fontsize=7, color="#8A8378", fontweight="bold")

    ax.set_title("Vegan Recipe — Countermeasure Nutrient Coverage", fontsize=14, fontweight="bold", pad=12)

    # Grid lines
    for i in range(mat.shape[0] + 1):
        ax.axhline(i - 0.5, color="white", linewidth=1)
    for j in range(mat.shape[1] + 1):
        ax.axvline(j - 0.5, color="white", linewidth=1)

    save_fig(fig, "fig4_nutrient_coverage")


# ──────────────────────────────────────────────────────────
# Figure 5: Recipe Scores by Phase
# ──────────────────────────────────────────────────────────
def fig5_recipe_scores():
    print("Generating fig5_recipe_scores...")
    recipes = pd.read_csv(os.path.join(PROJ, "results/tables/vegan_recipes_scored.csv"))
    recipes = recipes.sort_values("nutrient_score", ascending=True)

    phase_colors = {"pre_flight": PHYLO_BLUE, "in_flight": PHYLO_GREEN, "post_flight": PHYLO_ORANGE}

    fig, ax = plt.subplots(figsize=(9, 5.5))
    colors = [phase_colors.get(p, PHYLO_GOLD) for p in recipes["phase"]]
    bars = ax.barh(range(len(recipes)), recipes["nutrient_score"], color=colors, edgecolor="white", linewidth=0.5, height=0.7)

    ax.set_yticks(range(len(recipes)))
    ax.set_yticklabels(recipes["name"], fontsize=9)
    ax.set_xlabel("Nutrient Coverage Score", fontsize=11)
    ax.set_title("Vegan Recipe Scores by Mission Phase", fontsize=14, fontweight="bold", pad=12)

    for bar, val in zip(bars, recipes["nutrient_score"]):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2, str(int(val)),
                ha="left", va="center", fontsize=9, fontweight="bold")

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=PHYLO_BLUE, label="Pre-flight"),
                       Patch(facecolor=PHYLO_GREEN, label="In-flight"),
                       Patch(facecolor=PHYLO_ORANGE, label="Post-flight")]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9, frameon=False)

    save_fig(fig, "fig5_recipe_scores")


# ──────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 60)
    print("Script 07a: Figure Generation")
    print("=" * 60)

    fig1_signature_overview()
    fig2_top_reversers()
    fig3_gsea_hallmark()
    fig4_nutrient_coverage()
    fig5_recipe_scores()

    print("\nAll figures generated in:", FIG_DIR)
    print("Done.")
