# IBS Microbiome Analysis — Code Repository

This repository contains all analysis scripts and pipelines associated with the publication:

> **[Paper Title]**  
> [Author list]  
> *[Journal Name]*, [Year], DOI: [DOI]

---

## Overview

This study characterizes the gut microbiome dynamics in patients with Irritable Bowel Syndrome (IBS) across multiple clinical timepoints (V2, V4, V6, V7), comparing a treatment group with a placebo group. The analyses encompass raw sequence processing, taxonomic profiling, diversity analysis, differential abundance testing, co-occurrence network construction, and machine-learning–based biomarker discovery.

---

## Repository Structure

```
ibs-microbiome-analysis/
│
├── sequence-processing/                  # Raw sequencing QC and taxonomic profiling
│   └── sequencing_data_processing_commands.sh
│
├── diversity-analysis/
│   ├── alpha/                            # Alpha diversity (Shannon, Simpson)
│   │   ├── alpha-diversity-analysis.R
│   │   ├── Rel-abu-of-target-samples_in_paper_samples.txt   # Example input
│   │   └── metadata-of-target-samples_in_paper_samples.txt  # Example metadata
│   │
│   └── beta/                             # Beta diversity (Bray-Curtis PCoA)
│       ├── pcoa_analysis.R
│       ├── Rel-abu-of-target_in_paper_samples.txt           # Example input
│       ├── metadata-for-beta-diversity_in_paper_samples.txt # Example metadata
│       └── retrieve-given-stages/        # Per-stage PCoA visualisation
│           ├── plot-pcoA-matrix.R
│           └── PCoA-matrix.txt           # Example PCoA coordinates
│
├── differential-abundance-analysis/      # MaAsLin2, ANCOM-BC2, LEfSe, violin plots
│   ├── Differential-abundance-analysis.sh
│   └── violin-plot-script.R
│
├── microbial-association-network-analysis/  # SPIEC-EASI co-occurrence network
│   └── Network-analysis-using-SPIEC-EASI.R
│
├── random-forest-analysis/               # RF classifier with LOOCV (CLR-scaled)
│   └── RF_clr-scale_leave-one-out.py
│
└── data/                                 # Placeholder for input data files
    └── README.md
```

---

## Analysis Modules

### 1. Sequence Processing (`sequence-processing/`)

Raw metagenomics reads are processed in three sequential steps:

| Step | Tool | Version | Purpose |
|------|------|---------|---------|
| Quality control | BBDuk | v38.90 | Adapter trimming, quality filtering |
| Host decontamination | BBSplit | v38.90 | Remove human reads (CHM13 reference) |
| Taxonomic profiling | MetaPhlAn | 4.0.6 (db: vJun23_202403) | Species-level relative abundance |

**Script:** `sequencing_data_processing_commands.sh`  
**Input:** Paired-end FASTQ files (`R1.fastq.gz`, `R2.fastq.gz`)  
**Output:** MetaPhlAn taxonomic abundance profiles (`.txt`)

---

### 2. Diversity Analysis (`diversity-analysis/`)

#### Alpha Diversity (`alpha/`)

**Script:** `alpha-diversity-analysis.R`

Calculates Shannon and Simpson diversity indices from species-level relative abundance tables. Visualised as grouped boxplots by Study Group and Timepoint.

**Key R packages:** `vegan`, `ggplot2`, `ggpubr`

**Input files:**
- `Rel-abu-of-target-samples_in_paper_samples.txt` — species relative abundance (samples × species)
- `metadata-of-target-samples_in_paper_samples.txt` — sample metadata (Group, Stage, etc.)

#### Beta Diversity (`beta/`)

**Script:** `pcoa_analysis.R`

Computes Bray-Curtis dissimilarity and Principal Coordinates Analysis (PCoA). Points are coloured by timepoint and shaped by treatment group, with size scaled to NSR symptom score.

**Key R packages:** `vegan`, `BiodiversityR`, `ggplot2`, `RColorBrewer`

**Sub-analysis — Per-stage PCoA plots (`retrieve-given-stages/`):**  
`plot-pcoA-matrix.R` generates individual PCoA panels for each group × timepoint combination using pre-computed PCoA coordinates.

---

### 3. Differential Abundance Analysis (`differential-abundance-analysis/`)

Three complementary methods are used to ensure robustness:

| Method | Script | Key parameters |
|--------|--------|---------------|
| **LEfSe** | `Differential-abundance-analysis.sh` | LDA score > 2; FDR < 0.05 |
| **MaAsLin2** | `Differential-abundance-analysis.sh` | TSS normalisation; LOG transform; LM; reference stages: V2 / V4 / V6 |
| **ANCOM-BC2** | `Differential-abundance-analysis.sh` | BH adjustment; prevalence cut-off 10%; covariates: Stage, Age, Weight |

**Violin plot script:** `violin-plot-script.R`  
Generates per-species violin plots with statistical annotations (Wilcoxon between-group test; Kruskal-Wallis + Dunn's post-hoc within-group test; BH correction).

**Key R packages:** `Maaslin2`, `ANCOMBC`, `mia`, `ggstatsplot`, `ggplot2`, `FSA`

---

### 4. Microbial Co-occurrence Network Analysis (`microbial-association-network-analysis/`)

**Script:** `Network-analysis-using-SPIEC-EASI.R`

Constructs sparse inverse covariance networks using the SPIEC-EASI framework (MB neighbourhood selection). Edge weights are derived from the optimal precision matrix; edges are classified as positive or negative correlations. Outputs are formatted for Cytoscape visualisation.

**Key R packages:** `SpiecEasi`, `igraph`, `Matrix`

**Outputs:**
- `network_edges.txt` — edge list with weights and directions
- `network_nodes.txt` — node attributes with degree centrality

---

### 5. Random Forest Analysis (`random-forest-analysis/`)

**Script:** `RF_clr-scale_leave-one-out.py`

Trains a balanced Random Forest classifier to predict symptom response (stomachache). Features are CLR-transformed and filtered by variance. Top-*k* features (default: 50) are selected by initial feature importance before Leave-One-Out Cross-Validation (LOOCV).

**Outputs:**
- `filtered_features_with_names.tsv` — post-variance-filter feature matrix
- `top_feature_importances_loocv.csv` — ranked feature importances
- `confusion_matrix_loocv.pdf` — LOOCV confusion matrix
- `binary_roc_curve_loocv.pdf` — LOOCV ROC curve with AUC

**Usage:**
```bash
python RF_clr-scale_leave-one-out.py \
    --input_file your_data.tsv \
    --k_features 50 \
    --threads 8
```

**Required Python packages:** `pandas`, `numpy`, `scikit-learn`, `matplotlib`, `seaborn`

---

## Software Requirements

### Shell / HPC

| Software | Version | Installation |
|----------|---------|-------------|
| BBDuk / BBSplit | 38.90 | https://sourceforge.net/projects/bbmap/ |
| MetaPhlAn | 4.0.6 | `conda install -c bioconda metaphlan` |

### R (≥ 4.2.0)

```r
install.packages(c("vegan", "ggplot2", "dplyr", "tidyr", "RColorBrewer",
                   "ggthemes", "ggpubr", "data.table", "coin", "PMCMRplus", "FSA"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("Maaslin2", "ANCOMBC", "mia", "SpiecEasi"))

install.packages("BiodiversityR")   # May require manual installation
```

### Python (≥ 3.8)

```bash
pip install pandas numpy scikit-learn matplotlib seaborn
```

---

## Input Data Format

All analyses expect **tab-delimited** plain text files.

| File type | Format | Notes |
|-----------|--------|-------|
| Abundance table | Rows = samples, Columns = species | Relative abundance (0–1) or counts |
| Metadata | Rows = samples, Columns = variables | Required columns vary by script (see script headers) |
| Count table (ANCOM-BC2) | Rows = taxa, Columns = samples | Integer read counts |

Example input files for each module are provided within the respective subdirectory.

---

## Data Availability

Raw sequencing data are deposited at [database name, e.g., NCBI SRA] under accession number **[PRJNA/ERP accession]**.  
Processed abundance tables and metadata used directly by the scripts in this repository are available as Supplementary Data files in the associated publication.

---

## Citation

If you use these scripts in your research, please cite:

> [Full citation]

---

## License

This code is released under the [MIT License](LICENSE).

---

## Contact

For questions regarding the code, please open a GitHub Issue or contact the corresponding author:  
**[Corresponding Author Name]** — [email@institution.edu]
