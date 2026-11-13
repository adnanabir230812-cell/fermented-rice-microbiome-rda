"""
Master End-to-End Reproducibility Pipeline (Phases 1 to 4)
Manuscript: Microbial Diversity Dynamics and Nutritional Characterisation of Selected Fermented Sour Rice Varieties
Executes:
  - Phase 1: Alpha Diversity Non-Parametric Hypothesis Testing (Table S1)
  - Phase 2: Beta Diversity PERMANOVA, PERMDISP, ANOSIM & PCoA (Table S2, Figure 4b)
  - Phase 3: PICRUSt2 Predictive Functional Metagenomics (Table S3, Figure 8a, Figure 8b)
  - Phase 4: Spearman Rank Correlation & Constrained Multivariate RDA (Table S4, Figure 9, Figure 10)
"""

import os
import sys
import subprocess

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_DIR = os.path.join(REPO_ROOT, "scripts")

phases = [
    ("Phase 1: Alpha Diversity Statistical Significance (Table S1)", "phase1_alpha_diversity.py"),
    ("Phase 2: Beta Diversity PERMANOVA/PERMDISP & PCoA (Table S2, Fig 4b)", "phase2_beta_diversity_pcoa.py"),
    ("Phase 3: PICRUSt2 Predictive Functional Metagenomics (Table S3, Fig 8a, Fig 8b)", "phase3_picrust2_functional.py"),
    ("Phase 4: Microbiome-Nutritional Spearman & RDA Integration (Table S4, Fig 9, Fig 10)", "phase4_spearman_rda_integration.py"),
]

print("=" * 78)
print("  EXECUTING COMPLETE 4-PHASE BIOINFORMATICS & MULTIVARIATE PIPELINE")
print("=" * 78)

for title, script_name in phases:
    print(f"\n>>> Running {title}...")
    script_path = os.path.join(SCRIPTS_DIR, script_name)
    res = subprocess.run([sys.executable, script_path], cwd=REPO_ROOT)
    if res.returncode != 0:
        print(f"[ERROR] {script_name} exited with code {res.returncode}")
        sys.exit(res.returncode)

print("\n" + "=" * 78)
print("  ALL 4 PHASES COMPLETED SUCCESSFULLY!")
print("  All Tables (S1-S4) -> ./results/tables/")
print("  All 600 DPI Figures (4b, 8a, 8b, 9, 10) -> ./results/figures/")
print("=" * 78)
