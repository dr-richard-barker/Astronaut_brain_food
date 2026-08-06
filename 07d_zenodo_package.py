#!/usr/bin/env python3
"""
Script 07d: Assemble the Zenodo-ready package.

Creates a zip archive containing:
  - All scripts (01-07)
  - Results tables (CSV/JSON)
  - Figures (SVG + PNG)
  - README.md, CITATION.cff
  - PDF report
  - LaTeX manuscript
  - Data catalogs (not raw data — too large / license restrictions)

Output: /mnt/results/astronaut_opposite_forcing_zenodo.zip
"""

import os
import zipfile
import shutil
from datetime import datetime

PROJ = "/workspace/astronaut-opposite-forcing"
OUTPUT_ZIP = "/mnt/results/astronaut_opposite_forcing_zenodo.zip"

# Files to include (relative to project root)
INCLUDE_FILES = [
    # Scripts
    "scripts/01_OSD_query.R",
    "scripts/02_biomarker_identification.R",
    "scripts/02b_meta_analysis_resume.R",
    "scripts/03_drug_screening.py",
    "scripts/04_enrichment.R",
    "scripts/05_nutrient_gene_mapping.py",
    "scripts/06_recipe_optimization.py",
    "scripts/07_generate_figures.py",
    "scripts/07b_pdf_report.py",
    "scripts/07c_manuscript_readme.py",
    "scripts/07d_zenodo_package.py",
    # Data catalogs (small, license-clear)
    "data/study_catalog.csv",
    "data/download_status.csv",
    "data/human_study_catalog.csv",
    # Processed data (consensus signatures, meta-analysis results)
    "data/processed/consensus_signature_liver_immune.csv",
    "data/processed/consensus_signature_brain.csv",
    "data/processed/meta_analysis_liver_immune_full.csv",
    "data/processed/meta_analysis_brain_full.csv",
    "data/processed/ortholog_mapping.csv",
    "data/processed/study_sample_info.csv",
    "data/processed/signature_liver_immune_up_genes.txt",
    "data/processed/signature_liver_immune_down_genes.txt",
    "data/processed/signature_brain_up_genes.txt",
    "data/processed/signature_brain_down_genes.txt",
]

# Directories to include recursively
INCLUDE_DIRS = [
    "results/tables",
    "results/figures",
    "data/processed/per_study_de",
]

# Files from /mnt/results to include
RESULTS_FILES = [
    "/mnt/results/report_astronaut_opposite_forcing.pdf",
    "/mnt/results/manuscript_astronaut_opposite_forcing.tex",
    "/mnt/results/README.md",
    "/mnt/results/CITATION.cff",
]


def build_zenodo_package():
    print("=" * 60)
    print("Script 07d: Zenodo Package Assembly")
    print("=" * 60)

    # Write to /workspace first (zip is sequential, but let's be safe)
    temp_zip = os.path.join(PROJ, "astronaut_opposite_forcing_zenodo.zip")

    with zipfile.ZipFile(temp_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        # Individual files from project
        for fpath in INCLUDE_FILES:
            full = os.path.join(PROJ, fpath)
            if os.path.exists(full):
                arcname = f"astronaut-opposite-forcing/{fpath}"
                zf.write(full, arcname)
                print(f"  + {fpath}")
            else:
                print(f"  ! MISSING: {fpath}")

        # Directories
        for dpath in INCLUDE_DIRS:
            full_dir = os.path.join(PROJ, dpath)
            if not os.path.isdir(full_dir):
                print(f"  ! MISSING DIR: {dpath}")
                continue
            for root, dirs, files in os.walk(full_dir):
                for fname in files:
                    ffull = os.path.join(root, fname)
                    relpath = os.path.relpath(ffull, PROJ)
                    arcname = f"astronaut-opposite-forcing/{relpath}"
                    zf.write(ffull, arcname)
                    print(f"  + {relpath}")

        # Results files from /mnt/results
        for rpath in RESULTS_FILES:
            if os.path.exists(rpath):
                fname = os.path.basename(rpath)
                arcname = f"astronaut-opposite-forcing/{fname}"
                zf.write(rpath, arcname)
                print(f"  + {fname}")
            else:
                print(f"  ! MISSING: {rpath}")

    # Copy to /mnt/results (use shell cp to avoid S3 FUSE copystat issues)
    import subprocess
    subprocess.run(["cp", temp_zip, OUTPUT_ZIP], check=True)
    os.remove(temp_zip)

    size_mb = os.path.getsize(OUTPUT_ZIP) / (1024 * 1024)
    print(f"\nZenodo package: {OUTPUT_ZIP}")
    print(f"Size: {size_mb:.1f} MB")

    # Verify
    with zipfile.ZipFile(OUTPUT_ZIP, "r") as zf:
        names = zf.namelist()
        print(f"Files in archive: {len(names)}")
        # Check key files
        for key in ["README.md", "CITATION.cff", "report_astronaut_opposite_forcing.pdf"]:
            found = any(key in n for n in names)
            print(f"  {key}: {'OK' if found else 'MISSING'}")

    print("\nDone.")


if __name__ == "__main__":
    build_zenodo_package()
