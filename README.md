# Microbial Diversity Dynamics & Nutritional Characterisation of Fermented Sour Rice
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-276DC3.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)

This repository contains the standalone, fully reproducible bioinformatic and numerical ecology analysis pipeline accompanying the research manuscript:

> **"Microbial Diversity Dynamics and Nutritional Characterisation of Selected Fermented Sour Rice Varieties: A 16S rRNA Amplicon Sequencing Study"**  
> *Target Varieties:* ADT36, CO43, Mappillai Samba (MS), Improved White Ponni (IWP), White Ponni Unpolished (WPU), and CR1009 (Ration Rice) at 12 h and 16 h fermentation intervals.

---

## 🔬 Scientific Highlights & Workflow

1. **Compositional Data Harmonization:**
   * Genus-level relative abundance parsed across 16S amplicon libraries (`taxasums.xlsx`).
   * Triplicate nutritional and antioxidant profiling (`physicochemical_data.xlsx`).
2. **Spearman Rank Correlation & Multiple Testing Correction:**
   * Two-tailed non-parametric Spearman rank correlation ($r$) evaluating associations between 12 predominant bacterial genera and 13 nutritional/antioxidant parameters ($n = 6, df = 4$).
   * Multiplicity adjustments via **Benjamini–Hochberg False Discovery Rate (FDR $q$-values)**.
3. **Model Selection via Detrended Correspondence Analysis (DCA):**
   * DCA gradient length: **$2.6653\text{ SD} < 3.0\text{ SD}$**, strictly justifying a linear **Redundancy Analysis (RDA)** model over unimodal CCA (*Legendre & Legendre, 2012*).
4. **Constrained Multivariate Ordination (RDA):**
   * The canonical model captures **$75.58\%$ of total constrained variance** across the first two ordination axes (RDA1 = $52.68\%$, RDA2 = $22.90\%$), benchmarked against *Das et al. (AMB Express, 2026)*.
   * Proves that traditional pigmented *Mapillai Samba* clusters with protein, fiber, iron ($Fe$), magnesium ($Mg$), and DPPH scavenging activity via probiotic Firmicutes (*Bacillus*, *Lactococcus*) and enzymatic phytate dephosphorylation.

---

## 📁 Repository Structure

```text
fermented-rice-microbiome-rda/
├── data/                                 # Raw input matrices
│   ├── taxasums.xlsx                     # 16S amplicon genus relative abundance
│   └── physicochemical_data.xlsx         # Proximate, mineral, and antioxidant assays
├── notebooks/
│   └── reproducibility_analysis.ipynb   # Interactive Jupyter Notebook with inline plots
├── results/                              # Generated figures & statistical tables
│   ├── figures/                          # 600 DPI publication-grade figures (TIFF, PDF, PNG)
│   │   ├── Figure_9_Spearman_Microbiome_Nutrition_Heatmap.png
│   │   └── Figure_10_Multivariate_RDA_Biplot.png
│   └── tables/                           # Supplementary tables
│       └── Table_S4_Microbiome_Nutritional_Correlations_and_RDA.xlsx
├── scripts/
│   └── run_dca_rda.R                     # Standalone R script (vegan::decorana & vegan::rda)
├── run_analysis.py                       # Master one-click Python execution script
├── run_windows.bat                       # Double-click launcher for Windows
├── run_linux_mac.sh                      # Shell script for macOS / Linux
├── requirements.txt                      # Python dependencies
├── LICENSE                               # MIT License
└── README.md                             # Documentation & user guide
```

---

## ⚡ Quickstart: One-Click Reproduction

### Option A: Windows (Zero-Configuration)
Simply **double-click** `run_windows.bat` in the repository root directory. It will automatically initialize a virtual environment, install dependencies, execute the analysis, and save the 600 DPI figures and Excel tables into `./results/`.

### Option B: Command Line (Linux / macOS / Windows)
```bash
# 1. Clone the repository
git clone https://github.com/adnanabir230812-cell/fermented-rice-microbiome-rda.git
cd fermented-rice-microbiome-rda

# 2. Install dependencies
pip install -r requirements.txt

# 3. Execute the pipeline
python run_analysis.py
```

### Option C: R Users (vegan package)
```bash
Rscript scripts/run_dca_rda.R
```

---

## 📊 Summary of Deliverables

* **Figure 9:** Spearman rank correlation heatmap displaying significant positive correlations for *Idiomarina* ($r = +0.93^{**}$ with fat) and *Weissella* ($r = +0.75$ with DPPH), alongside inverse correlations for opportunistic *Ralstonia* ($r = -0.85^*$ with Mg and DPPH).
* **Figure 10:** 2D multivariate RDA biplot illustrating niche partitioning among fermented rice varieties under nutritional constraints.
* **Supplementary Table S4:** Multi-sheet Excel workbook containing Spearman correlation coefficients ($r$), exact $P$-values, FDR $q$-values, RDA site scores, species loadings, and environmental biplot vectors.

---

## 📜 Citation & Code Availability

If using this code or dataset, please cite:
```bibtex
@article{sour_rice_microbiome_2026,
  title={Microbial Diversity Dynamics and Nutritional Characterisation of Selected Fermented Sour Rice Varieties: A 16S rRNA Amplicon Sequencing Study},
  author={Lead Authors},
  journal={AMB Express},
  year={2026},
  publisher={Springer Nature}
}
```

*Permanent Archive:* Zenodo DOI: `10.5281/zenodo.1234567` (or corresponding release DOI).
