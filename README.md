# Astronaut Opposite Forcing

## Spaceflight Transcriptomic Signatures and Vegan Nutritional Countermeasures

A reproducible computational pipeline that:
1. Builds consensus spaceflight transcriptomic signatures from NASA OSDR rodent ISS RNA-seq studies
2. Screens LINCS L1000 perturbation compendium for compounds that reverse the signatures ("opposite forcings")
3. Annotates hits with DrugBank 5.1.8 (open-access), Broad Drug Repurposing Hub, and ChEMBL
4. Validates biological plausibility via functional enrichment (GSEA/ORA) and JensenLab DISEASES
5. Translates therapeutic candidates into a vegan nutrition layer (FooDB + NutriChem 2.0)
6. Produces curated vegan recipes organized by mission phase (pre/in/post-flight)

## Key Results

| Metric | Value |
|--------|-------|
| Studies screened | 160 NASA OSDR rodent spaceflight RNA-seq |
| Studies meta-analyzed | 15 (with valid DE signal) |
| Liver/Immune signature | 525 consensus genes (283 up, 242 down), 25 core |
| Brain signature | 505 consensus genes (208 up, 297 down), 13 core |
| LINCS reversers | 219 liver, 221 brain candidate compounds |
| Vegan countermeasure nutrients | 10 bioactive compounds |
| Curated recipes | 12 phase-optimized vegan meals |

## Pipeline Structure

```
astronaut-opposite-forcing/
├── scripts/
│   ├── 01_OSD_query.R              # NASA OSDR study enumeration & data download
│   ├── 02_biomarker_identification.R  # Per-study DESeq2 differential expression
│   ├── 02b_meta_analysis_resume.R  # Random-effects meta-analysis + consensus signature
│   ├── 03_drug_screening.py        # LINCS L1000 connectivity mapping + drug annotation
│   ├── 04_enrichment.R             # GSEA/ORA + JensenLab DISEASES enrichment
│   ├── 05_nutrient_gene_mapping.py # FooDB compound-food mapping
│   ├── 06_recipe_optimization.py   # Vegan recipe scoring & meal plan generation
│   ├── 07_generate_figures.py      # Figure generation (SVG + PNG)
│   ├── 07b_pdf_report.py           # Phylo-branded PDF report
│   └── 07c_manuscript_readme.py    # LaTeX manuscript, README, CITATION.cff
├── data/
│   ├── study_catalog.csv           # Study metadata
│   ├── download_status.csv         # Download tracking
│   ├── human_study_catalog.csv     # Human studies for cross-check
│   ├── raw/
│   │   ├── foodb/                  # FooDB 2020-04-07 CSV release
│   │   └── jensenlab_diseases_*.tsv
│   └── processed/
│       ├── consensus_signature_{liver_immune,brain}.csv
│       ├── meta_analysis_{liver_immune,brain}_full.csv
│       ├── ortholog_mapping.csv
│       ├── per_study_de/           # Per-study DE results
│       └── signature_{stratum}_{up,down}_genes.txt
├── results/
│   ├── tables/
│   │   ├── lincs_{liver_immune,brain}/
│   │   │   ├── tier1_ranking.csv
│   │   │   ├── tier1_ranking_annotated.csv
│   │   │   ├── reversers_raw.csv
│   │   │   └── robustness_summary.json
│   │   ├── gsea_{stratum}_{hallmark,reactome}.csv
│   │   ├── ora_{stratum}_{up,dn}_{gobp,reactome}.csv
│   │   ├── jensen_diseases_{brain,liver_immune}.csv
│   │   ├── nutrient_gene_map_{liver_immune,brain}.csv
│   │   ├── nutrient_gene_map_combined.csv
│   │   ├── vegan_recipes_scored.csv
│   │   ├── vegan_meal_plans.csv
│   │   └── countermeasure_nutrient_summary.csv
│   └── figures/
│       ├── fig1_signature_overview.{svg,png}
│       ├── fig2_top_reversers.{svg,png}
│       ├── fig3_gsea_hallmark.{svg,png}
│       ├── fig4_nutrient_coverage.{svg,png}
│       └── fig5_recipe_scores.{svg,png}
└── report_astronaut_opposite_forcing.pdf
```

## Dependencies

### R (>= 4.3)
- DESeq2, apeglm
- metafor
- clusterProfiler, msigdbr
- org.Mm.eg.db, homologene
- AnnotationDbi (note: use `dplyr::select()` explicitly to avoid namespace conflict)

### Python (>= 3.10)
- pandas, numpy
- matplotlib, seaborn
- reportlab, pypdf
- PIL/Pillow

## Data Sources

| Source | URL | License |
|--------|-----|---------|
| NASA OSDR | https://osdr.nasa.gov | Public domain |
| LINCS L1000 / SigCom | https://maayanlab.cloud/sigcom-lincs | Public |
| Broad Drug Repurposing Hub | https://clue.io/repurposing-app | CC BY 4.0 |
| DrugBank 5.1.8 | https://www.drugbank.ca | CC BY-NC 4.0 (academic) |
| ChEMBL | https://www.ebi.ac.uk/chembl | CC BY-SA 3.0 |
| FooDB | https://foodb.ca | CC BY-NC 4.0 |
| JensenLab DISEASES | https://diseases.jensenlab.org | CC BY 4.0 |
| MSigDB Hallmark | https://www.gsea-msigdb.org | Free for academic |

## Reproduction

```bash
# R scripts
Rscript scripts/01_OSD_query.R
Rscript scripts/02_biomarker_identification.R
Rscript scripts/02b_meta_analysis_resume.R
Rscript scripts/04_enrichment.R

# Python scripts
python scripts/03_drug_screening.py
python scripts/05_nutrient_gene_mapping.py
python scripts/06_recipe_optimization.py
python scripts/07_generate_figures.py
python scripts/07b_pdf_report.py
python scripts/07c_manuscript_readme.py
```

## Citation

See CITATION.cff

## License

- Code: MIT License
- Data: Respective source licenses (see Data Sources table)
- Report: CC BY 4.0

## Contact

Generated by Biomni (Phylo, Inc.)
