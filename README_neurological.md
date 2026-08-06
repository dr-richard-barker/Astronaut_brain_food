# Neurological Countermeasures for Spaceflight

## A Transcriptomic-Driven, Diet-Based Neuroprotection Strategy

This package presents a neurological-focused reanalysis of rodent ISS spaceflight brain transcriptomic data, combined with LINCS L1000 drug connectivity mapping, nutrient pharmacokinetic profiling, and vegan dietary translation.

## Contents

### Reports
- `report_neurological_countermeasures.pdf` — Full PDF report (~20 pages)
- `manuscript_neurological_countermeasures.tex` — LaTeX manuscript

### Figures (results/figures/)
- `fig1_brain_volcano.svg/png` — Brain transcriptomic volcano plot
- `fig2_gsea_theme_dotplot.svg/png` — GSEA pathways by neurological theme
- `fig3_lincs_reversers_bar.svg/png` — Top 15 LINCS reversers
- `fig4_nutrient_pk_radar.svg/png` — Nutrient PK comparison radar
- `fig5_spaceflight_brain_mechanism.png` — Spaceflight brain injury cascade
- `fig6_nutrient_intervention_network.png` — Nutrient-pathway-gene network
- `fig7_diet_mission_phase.png` — Diet strategy by mission phase
- `fig8_recipe_nutrient_heatmap.svg/png` — Recipe nutrient coverage heatmap
- `figure_legends.md` — Detailed figure legends

### Tables (results/tables/)
- `brain_pathway_themes.csv` — 550 pathways in 10 neurological themes
- `brain_nutrient_pk_table.csv` — 7 nutrients × 16 PK columns
- `brain_gene_pathway_nutrient_network.csv` — 27-edge intervention network
- `brain_supplementary_tables/` — Tables S1-S6 (full signatures, GSEA, LINCS, ORA)

### Scripts (scripts/)
- `08_brain_neuro_analysis.py` — Pathway themes, PK table, network
- `09_neuro_figures.py` — Data-driven figures (1-4, 8)
- `10_neuro_report.py` — PDF report + LaTeX manuscript
- `11_neuro_zenodo.py` — Zenodo package assembly

## Key Findings

1. **Brain consensus signature**: 505 genes (208 up, 297 down, 13 core) from 4 rodent ISS studies
2. **Dominant pathology**: Mitochondrial dysfunction, protein homeostasis failure, DNA repair suppression
3. **Top nutrient reversers**: Sulforaphane (Nrf2 activator), Ascorbic acid (BBB transport), Quercetin (BBB integrity)
4. **Best PK profiles**: Ascorbic acid and Sulforaphane (confirmed BBB penetration, therapeutic brain concentrations)
5. **Critical caveat**: Curcumin bioavailability ~1000x below therapeutic (Kroon et al. 2025)

## Data Sources

- NASA OSDR: OSD-525, OSD-564, OSD-612, OSD-613 (rodent ISS brain RNA-seq)
- LINCS L1000: NIH Common Fund (~1.8M perturbational profiles)
- Broad Drug Repurposing Hub: Clinical phase and MOA annotations
- FooDB: Vegan food source mapping
- Published literature: Nutrient pharmacokinetic data

## Reproducibility

All analyses are reproducible from the provided scripts. Run in order:
```
python scripts/08_brain_neuro_analysis.py
python scripts/09_neuro_figures.py
python scripts/10_neuro_report.py
python scripts/11_neuro_zenodo.py
```

## Citation

See `CITATION_neurological.cff` for citation information.
