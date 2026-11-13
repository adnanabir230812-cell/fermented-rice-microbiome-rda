# Microbial Diversity Dynamics & Nutritional Characterisation of Fermented Sour Rice
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-276DC3.svg)](https://www.r-project.org/)
[![PICRUSt2](https://img.shields.io/badge/PICRUSt2-v2.5.2-green.svg)](https://github.com/picrust/picrust2)

This repository contains the complete, end-to-end reproducible bioinformatic, functional metagenomic (PICRUSt2), and multivariate ecological analysis pipeline (Phases 1 to 4) accompanying the research manuscript:

> **"Microbial Diversity Dynamics and Nutritional Characterisation of Selected Fermented Sour Rice Varieties: A 16S rRNA Amplicon Sequencing Study"**  
> *Target Varieties:* ADT36, CO43, Mappillai Samba (MS), Improved White Ponni (IWP), White Ponni Unpolished (WPU), and CR1009 (Ration Rice) at 12 h and 16 h fermentation intervals.

---

## 🔬 Complete 4-Phase Analytical Workflow

### Phase 1: Alpha Diversity Formal Statistical Hypothesis Testing (`scripts/phase1_alpha_diversity.py`)
* **Metrics Evaluated:** Observed OTUs, Chao1, Shannon ($H'$), Gini-Simpson ($1-D$), Pielou's Evenness ($J'$), Fisher's $\alpha$, and Faith's Phylogenetic Diversity (PD).
* **Statistical Tests:** Non-parametric paired Wilcoxon signed-rank test (12 h vs. 16 h temporal shifts) and omnibus Kruskal–Wallis test across the 6 cultivars followed by Dunn's post-hoc pairwise tests with Benjamini–Hochberg False Discovery Rate (FDR $q$-values).
* **Output:** `results/tables/Table_S1_Alpha_Diversity_Statistical_Significance.xlsx`

### Phase 2: Beta Diversity, Two-Way PERMANOVA, PERMDISP & PCoA (`scripts/phase2_beta_diversity_pcoa.py` & `scripts/phase2_permanova_vegan.R`)
* **Distance Matrices (`data/distance_matrices/`):** Weighted UniFrac, Unweighted UniFrac, Bray–Curtis, and Jaccard.
* **Multivariate Tests:** Two-way sequential PERMANOVA (`vegan::adonis2`, 9,999 permutations: $\text{Distance} \sim \text{Variety} + \text{Time}$), multivariate homogeneity of group dispersions (`vegan::betadisper` / PERMDISP), and ANOSIM.
* **Outputs:**
  * `results/figures/Figure_4b_PCoA_Weighted_UniFrac` (600 DPI TIFF, PDF, PNG)
  * `results/tables/Table_S2_Beta_Diversity_PERMANOVA_and_PERMDISP.xlsx`

### Phase 3: PICRUSt2 Predictive Functional Metagenomics (`scripts/phase3_picrust2_functional.py`)
* **Inputs:** 16S rRNA amplicon sequences (`data/dna-sequences.fasta`), BIOM table (`data/feature-table.biom`), and taxonomic profiles (`data/taxasums.xlsx`) mapped to KEGG Orthology (KO) and MetaCyc pathways across 23 core metabolic pathways (benchmarked against *Shruthi et al., 2025* and *Shangpliang et al., 2023*).
* **Statistical Profiling:** Standardized functional abundance ($Z$-scores) and temporal $\log_2$ fold-change ($\log_2\text{FC}$, 16 h vs. 12 h) with paired Wilcoxon and FDR adjustments.
* **Outputs:**
  * `results/figures/Figure_8a_Functional_Pathway_Heatmap` (600 DPI TIFF, PDF, PNG)
  * `results/figures/Figure_8b_Differential_Pathway_Enrichment` (600 DPI TIFF, PDF, PNG)
  * `results/tables/Table_S3_Predicted_Metabolic_Pathways.xlsx`

### Phase 4: Microbiome–Nutritional Integration via Spearman Correlation & RDA (`scripts/phase4_spearman_rda_integration.py` & `scripts/phase4_dca_rda_vegan.R`)
* **Spearman Rank Correlation:** Two-tailed Spearman correlation ($r$) and Benjamini–Hochberg FDR ($q$-values) coupling 12 dominant bacterial genera with 13 nutritional, mineral, and antioxidant parameters (`data/physicochemical_data.xlsx`).
* **DCA & Redundancy Analysis (RDA):** Detrended Correspondence Analysis (`vegan::decorana`) confirmed a linear response gradient ($\text{DCA1} = 2.665\text{ SD} < 3.0\text{ SD}$), validating linear Redundancy Analysis (`vegan::rda`), which captured **75.58% of total constrained variance** (RDA1 = 52.68%, RDA2 = 22.90%; benchmarked against *Das et al., AMB Express, 2026*).
* **Outputs:**
  * `results/figures/Figure_9_Spearman_Microbiome_Nutrition_Heatmap` (600 DPI TIFF, PDF, PNG)
  * `results/figures/Figure_10_Multivariate_RDA_Biplot` (600 DPI TIFF, PDF, PNG)
  * `results/tables/Table_S4_Microbiome_Nutritional_Correlations_and_RDA.xlsx`

---

## 📁 Repository Structure

```text
fermented-rice-microbiome-rda/
├── data/                                               # Raw 16S QIIME2/BIOM & nutritional input files
│   ├── distance_matrices/                              # Beta diversity distance matrices (.tsv)
│   │   ├── weighted_unifrac.tsv
│   │   ├── unweighted_unifrac.tsv
│   │   ├── bray_curtis.tsv
│   │   └── jaccard.tsv
│   ├── dna-sequences.fasta                             # 16S rRNA representative amplicon sequences
│   ├── feature-table.biom                              # BIOM v2.1 observation matrix (for PICRUSt2)
│   ├── feature-table.tsv                               # Tab-separated OTU/ASV count table
│   ├── metadata.tsv                                    # Sample metadata (Cultivar & Fermentation Time)
│   ├── taxonomy.qza                                    # QIIME2 taxonomy artifact
│   ├── taxa.xlsx                                       # Multi-level taxonomic classifications
│   ├── taxasums.xlsx                                   # Genus-level abundance summary
│   └── physicochemical_data.xlsx                       # Triplicate nutritional, mineral & antioxidant data
├── notebooks/
│   └── reproducibility_analysis.ipynb                 # Pre-rendered Jupyter Notebook with inline figures
├── results/
│   ├── figures/                                        # All 600 DPI publication figures (TIFF, PDF, PNG)
│   │   ├── Figure_4b_PCoA_Weighted_UniFrac.*
│   │   ├── Figure_8a_Functional_Pathway_Heatmap.*
│   │   ├── Figure_8b_Differential_Pathway_Enrichment.*
│   │   ├── Figure_9_Spearman_Microbiome_Nutrition_Heatmap.*
│   │   └── Figure_10_Multivariate_RDA_Biplot.*
│   └── tables/                                         # All Supplementary Excel Tables (S1 to S4)
│       ├── Table_S1_Alpha_Diversity_Statistical_Significance.xlsx
│       ├── Table_S2_Beta_Diversity_PERMANOVA_and_PERMDISP.xlsx
│       ├── Table_S3_Predicted_Metabolic_Pathways.xlsx
│       └── Table_S4_Microbiome_Nutritional_Correlations_and_RDA.xlsx
├── scripts/                                            # Modular Phase-specific Python & R scripts
│   ├── phase1_alpha_diversity.py                       # Phase 1: Alpha diversity non-parametric stats
│   ├── phase2_beta_diversity_pcoa.py                   # Phase 2: Beta diversity PERMANOVA, PERMDISP & Fig 4b
│   ├── phase2_permanova_vegan.R                        # Phase 2: R vegan::adonis2 & betadisper verification
│   ├── phase3_picrust2_functional.py                   # Phase 3: PICRUSt2 KEGG/MetaCyc profiling & Fig 8a/8b
│   ├── phase4_spearman_rda_integration.py              # Phase 4: Spearman heatmap & RDA biplot (Fig 9 & 10)
│   └── phase4_dca_rda_vegan.R                          # Phase 4: R vegan::decorana (DCA) & rda verification
├── run_analysis.py                                     # Master all-in-one runner executing Phases 1 to 4
├── run_windows.bat                                     # One-click Windows batch launcher
├── run_linux_mac.sh                                    # One-click Linux/macOS shell launcher
├── requirements.txt                                    # Python package dependencies
└── LICENSE                                             # MIT Open-Source License
```

---

## ⚡ Quickstart: One-Click Execution (All 4 Phases)

### Option A: Windows (Zero-Configuration)
Double-click `run_windows.bat`. It will set up the environment, install dependencies, and execute **all 4 phases** (`run_analysis.py`), regenerating Tables S1–S4 and Figures 4b, 8a, 8b, 9, and 10 inside `./results/`.

### Option B: Command Line (Linux / macOS / Windows)
```bash
git clone https://github.com/adnanabir230812-cell/fermented-rice-microbiome-rda.git
cd fermented-rice-microbiome-rda
pip install -r requirements.txt

# Run the entire 4-phase pipeline at once:
python run_analysis.py

# Or run any individual phase script:
python scripts/phase1_alpha_diversity.py
python scripts/phase2_beta_diversity_pcoa.py
python scripts/phase3_picrust2_functional.py
python scripts/phase4_spearman_rda_integration.py
```

### Option C: R Verification Scripts (`vegan` package)
```bash
Rscript scripts/phase2_permanova_vegan.R
Rscript scripts/phase4_dca_rda_vegan.R
```
