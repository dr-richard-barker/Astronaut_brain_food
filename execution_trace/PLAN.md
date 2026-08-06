# Plan: Neurological-Focused Reanalysis — Preserving Brain Function with Diet

## Summary

Refocus the astronaut-opposite-forcing pipeline on the **neurological narrative**: how spaceflight damages the brain at the transcriptomic level, and how vegan nutritional countermeasures can preserve neurological function. Re-use the existing brain consensus signature (505 genes), GSEA/ORA results, and LINCS brain reversers — no new data acquisition needed. Add a full pharmacokinetic layer grounded in published literature, produce a comprehensive set of primary + supplementary tables, mechanism + data figures with legends, and update all deliverables (PDF report, LaTeX manuscript, Zenodo package).

## Key Changes from Previous Run

1. **Brain-only focus**: Drop liver/immune stratum from the narrative. The brain signature (505 genes, 13 core) becomes the sole subject. Liver results remain in the Zenodo package as supplementary context but are not in the main report.

2. **New pharmacokinetic analysis**: Build a nutrient PK table with literature-grounded data for 7 brain-relevant nutrients (sulforaphane, curcumin, luteolin, EGCG, ascorbic acid, quercetin, resveratrol) covering: oral bioavailability, peak plasma concentration, time to peak, half-life, BBB penetration evidence, brain tissue concentration, recommended human dose, spaceflight-relevant notes. All values cited to specific papers.

3. **New figures (8 total)**: 4 data-driven + 4 mechanism-focused. All with formal figure legends.

4. **New tables (6 primary + 4 supplementary)**: All with formal table legends.

5. **Updated PDF report**: Restructured around the neurological narrative with 6 sections.

6. **Updated LaTeX manuscript**: Brain-focused, npj Microgravity style, with PK data.

7. **Updated Zenodo package**: All new tables, figures, scripts, updated README/CITATION.

---

## Section 1: New Analyses (Script 08)

### 1A. Brain pathway deep-dive (Python, from existing data)

Extract and categorize all significant brain GSEA Reactome pathways (527 total, p.adj < 0.25) into neurological theme groups:

- **Mitochondrial dysfunction** (Complex I biogenesis, respiratory electron transport, mitochondrial translation, Complex IV assembly, mitochondrial protein import, aerobic respiration)
- **Synaptic/neuronal compensation** (neuronal system, transmission across chemical synapses, NMDA receptors, neurotransmitter release cycle, long-term potentiation, protein-protein interactions at synapses)
- **DNA damage/repair** (DNA repair, base excision repair, nucleotide excision repair, telomere extension)
- **Protein homeostasis** (unfolded protein response, proteasome assembly, mRNA splicing, translation)
- **Circadian/cholesterol dysregulation** (BMAL1/CLOCK, RORA/NR1D1, SREBP cholesterol biosynthesis)
- **Cell cycle suppression** (E2F targets, G2M checkpoint, MYC targets, mitotic spindle)

Output: `results/tables/brain_pathway_themes.csv` (pathway, theme, NES, p.adj, direction)

### 1B. Nutrient-pathway mapping (Python, from existing data + literature)

For each of the 7 brain-relevant nutrients, map to:
- Affected neurological pathways (from GSEA results + literature mechanism)
- LINCS reversal evidence (z-sum, rank)
- FooDB vegan food sources (from existing nutrient_gene_map_brain.csv)
- Pharmacokinetic data (from literature, cited)

Output: `results/tables/brain_nutrient_pk_table.csv` — the core PK table with columns:
  nutrient, mechanism, target_pathways, bbb_penetration, oral_bioavailability, peak_plasma, t_peak, half_life, brain_concentration, human_dose, spaceflight_relevance, key_citations

### 1C. Brain gene-pathway-nutrient network (Python)

Build a tripartite network: brain DE genes → affected pathways → countermeasure nutrients. Used for the network figure and supplementary table.

Output: `results/tables/brain_gene_pathway_nutrient_network.csv`

---

## Section 2: Figures (Script 09)

All figures use Phylo brand colors, Liberation Sans, SVG + PNG, with formal legends saved to `results/figures/figure_legends.md`.

### Data-driven figures (matplotlib/seaborn)

**Figure 1: Brain signature volcano plot**
- X: pooled log2FC, Y: -log10(p-value)
- Color: up (orange) / down (blue) / non-significant (gray)
- Label top 15 genes with human symbols
- Data: `data/processed/meta_analysis_brain_full.csv`

**Figure 2: Brain GSEA Hallmark + Reactome combined dotplot**
- Y: pathway names (top 20 by |NES|, grouped by theme)
- X: NES (negative = downregulated = red, positive = upregulated = blue)
- Dot size: -log10(p.adj)
- Data: `results/tables/gsea_brain_hallmark.csv` + `gsea_brain_reactome.csv`

**Figure 3: Brain LINCS reversers — top 15 horizontal bar chart**
- Y: compound names
- X: best z-sum (negative = reversal)
- Color by clinical phase
- Annotate with MoA
- Data: `results/tables/lincs_brain/tier1_ranking_annotated.csv`

**Figure 4: Nutrient pharmacokinetic comparison radar/spider chart**
- 7 nutrients on 6 axes: BBB penetration, oral bioavailability, brain tissue accumulation, half-life, evidence strength, spaceflight relevance
- Each nutrient as a separate colored line
- Data: `results/tables/brain_nutrient_pk_table.csv` (normalized 0-1 scores)

### Mechanism-focused figures (GenerateImage)

**Figure 5: Spaceflight brain injury mechanism cascade (schematic)**
- Microgravity/radiation → oxidative stress → mitochondrial dysfunction (Complex I ↓) → ATP depletion → synaptic impairment → neuroinflammation → BBB disruption → cognitive decline
- Show where each nutrient intervenes (Nrf2 activation, mitochondrial biogenesis, anti-inflammatory, BBB protection)
- Generated via GenerateImage

**Figure 6: Nutrient → pathway → gene intervention network (schematic)**
- Central: "Spaceflight Brain Signature"
- Left: downregulated pathways (mitochondrial, synaptic)
- Right: 7 nutrients with their target mechanisms
- Arrows showing reversal direction
- Generated via GenerateImage

**Figure 7: Diet-neuroprotection strategy by mission phase (schematic)**
- Three columns: Pre-flight / In-flight / Post-flight
- Each with: neurological goal, key nutrients, representative recipe, PK consideration
- Generated via GenerateImage

**Figure 8: Vegan recipe nutrient coverage heatmap (data-driven)**
- Re-use existing fig4_nutrient_coverage but brain-focused (filter to brain-relevant nutrients only)
- Add PK bioavailability column annotation
- Data: `results/tables/vegan_recipes_scored.csv` + PK table

---

## Section 3: Tables (Script 08/09)

### Primary Tables (in PDF report + CSV)

**Table 1: Brain consensus signature summary**
- Stratum, total genes, up, down, core, n studies, ortholog coverage
- Legend: "Consensus spaceflight brain transcriptomic signature derived from random-effects meta-analysis of 4 rodent ISS studies."

**Table 2: Top 20 brain LINCS reversers**
- Compound, z-sum, clinical phase, MoA, target, disease area, FooDB match, is_nutrient
- Legend: "Top compounds reversing the brain spaceflight signature in LINCS L1000 connectivity mapping."

**Table 3: Brain GSEA pathway themes**
- Theme, pathway, NES, p.adj, direction, key genes
- Legend: "Significant GSEA pathways grouped by neurological theme (p.adj < 0.25)."

**Table 4: Brain-relevant nutrient pharmacokinetic profile**
- Nutrient, mechanism, BBB penetration, oral bioavailability, peak plasma, t_peak, half-life, brain concentration, human dose, spaceflight relevance, citations
- Legend: "Pharmacokinetic and neuroprotective profile of brain-relevant vegan countermeasure nutrients, grounded in published literature."

**Table 5: Vegan recipe nutrient coverage by mission phase**
- Recipe, phase, target nutrients, nutrient score, PK notes
- Legend: "Curated vegan recipes organized by mission phase with nutrient coverage and pharmacokinetic considerations."

**Table 6: Nutrient → pathway → gene intervention map**
- Nutrient, target pathway, target genes (from brain signature), mechanism, evidence type
- Legend: "Mapping of countermeasure nutrients to affected neurological pathways and brain signature genes."

### Supplementary Tables (CSV only, in Zenodo package)

**Supplementary Table S1: Full brain consensus signature (505 genes)**
- gene_id, mouse_symbol, human_symbol, pooled_log2FC, pvalue, padj, direction, I2, k, n_studies_up, n_studies_down, core

**Supplementary Table S2: Full brain GSEA Hallmark results (all pathways)**
- ID, Description, setSize, ES, NES, pvalue, p.adjust, leading_edge, core_enrichment

**Supplementary Table S3: Full brain GSEA Reactome results (all significant)**
- Same columns as S2, filtered to p.adj < 0.25

**Supplementary Table S4: Full brain LINCS reverser ranking (all compounds)**
- rank, compound, n_reversing_sigs, median_z_sum, best_z_sum, best_p, method, DRH annotation fields

---

## Section 4: PDF Report (Script 10)

Phylo-branded PDF using ReportLab per pdf-report-generation skill. ~20-25 pages.

### Structure:
1. **Title page** — "Preserving Neurological Function in Spaceflight: A Diet-Based Countermeasure Strategy"
2. **Executive Summary** — 3 paragraphs + key metrics callout
3. **Introduction** — Spaceflight brain biology, mitochondrial dysfunction, BBB damage
4. **Methods** — Meta-analysis, LINCS screening, enrichment, nutrient PK literature review
5. **Results**
   - 5.1 Brain consensus signature (Fig 1, Table 1)
   - 5.2 Pathway analysis (Fig 2, Table 3)
   - 5.3 LINCS reversers (Fig 3, Table 2)
   - 5.4 Nutrient pharmacokinetics (Fig 4, Table 4)
   - 5.5 Nutrient-pathway-gene network (Fig 6, Table 6)
   - 5.6 Vegan recipe optimization (Fig 8, Table 5)
   - 5.7 Mechanism cascade (Fig 5, Fig 7)
6. **Discussion** — Biological interpretation, PK limitations (especially curcumin), spaceflight-specific considerations
7. **Conclusions & Next Steps**
8. **References** — All literature citations

Output: `/mnt/results/report_neurological_countermeasures.pdf`

---

## Section 5: LaTeX Manuscript (Script 10)

Updated npj Microgravity-style manuscript, brain-focused.

### Structure:
- Title: "Dietary preservation of neurological function during spaceflight: a computational pipeline linking brain transcriptomic signatures to vegan nutritional countermeasures"
- Abstract (brain-focused)
- Introduction (spaceflight neurobiology)
- Methods (meta-analysis, LINCS, enrichment, PK literature review)
- Results (signature, pathways, reversers, nutrient PK, recipes)
- Discussion (mitochondrial dysfunction, BBB, PK limitations)
- Conclusions
- References (expanded with PK papers)

Output: `/mnt/results/manuscript_neurological_countermeasures.tex`

---

## Section 6: Zenodo Package (Script 11)

Updated package with all new artifacts.

### Contents:
- All scripts (01-07 + 08-11)
- All result tables (existing + new brain-focused)
- All figures (existing + 8 new with legends)
- PDF report (neurological version)
- LaTeX manuscript (neurological version)
- README.md (updated, brain-focused)
- CITATION.cff (updated)
- Figure legends file

Output: `/mnt/results/astronaut_neurological_zenodo.zip`

---

## Execution Plan

| Step | Script | Runtime | Output |
|------|--------|---------|--------|
| 1 | `08_brain_neuro_analysis.py` | ~5 min | New tables (pathway themes, PK, network) |
| 2 | `09_neuro_figures.py` | ~10 min | 8 figures (SVG+PNG) + figure legends |
| 3 | `10_neuro_report.py` | ~5 min | PDF report + LaTeX manuscript |
| 4 | `11_neuro_zenodo.py` | ~5 min | Zenodo package zip |

**Total estimated runtime: ~25 min** (all foreground, default machine)

---

## Key Assumptions

1. **Curcumin honesty**: The 2025 critical reappraisal [12] shows unconjugated curcumin plasma levels remain ~1000-fold below in vitro effective concentrations, and piperine does NOT improve bioavailability. The report will present curcumin's neuroprotective mechanism evidence while honestly flagging this PK limitation. This is a critical scientific integrity point.

2. **Sulforaphane as lead nutrient**: Best BBB penetration evidence (crosses BBB, accumulates in ventral midbrain/striatum within 15 min), Nrf2 activation directly addresses oxidative stress/mitochondrial dysfunction, and is the top brain-specific nutrient from LINCS. Will be positioned as the primary dietary countermeasure.

3. **Vitamin C brain gradient**: Plasma ~50 µM → CSF ~200 µM → higher intracellular brain. SVCT2 transporter actively concentrates vitamin C in brain. This is the best-documented brain PK profile among the nutrients.

4. **No new data acquisition**: All analyses re-use existing brain consensus signature, GSEA/ORA, and LINCS results. New work is synthesis, visualization, PK literature integration, and report generation.

5. **GenerateImage for mechanism figures**: Figures 5, 6, 7 are schematic/conceptual (mechanism cascades, network diagrams, strategy infographics) — these require GenerateImage per visualization guidelines, not matplotlib.

6. **Figure legends**: All 8 figures get formal legends in a separate `figure_legends.md` file and embedded in the PDF report.

7. **Citation format**: Inline [N] references matching the literature search document indices, with full reference list in the PDF and LaTeX.
