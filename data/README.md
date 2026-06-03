# Data Directory

This directory is a placeholder for input data files required by the analysis scripts.

## Data Availability

Raw sequencing data are deposited at [NCBI SRA / ENA / DDBJ] under accession number **[accession]**.

Processed abundance tables and metadata files are available as **Supplementary Data** in the associated publication. Download the relevant files and place them in the corresponding analysis subdirectory before running the scripts.

## File Placement Guide

| File(s) | Place in directory |
|---------|-------------------|
| Paired-end FASTQ files | `sequence-processing/` |
| Species relative abundance table (alpha) | `diversity-analysis/alpha/` |
| Sample metadata (alpha) | `diversity-analysis/alpha/` |
| Species relative abundance table (beta) | `diversity-analysis/beta/` |
| Sample metadata (beta) | `diversity-analysis/beta/` |
| Abundance table (MaAsLin2) | `differential-abundance-analysis/` |
| Count-based abundance tables (ANCOM-BC2) | `differential-abundance-analysis/` |
| Species-level count abundance (SPIEC-EASI) | `microbial-association-network-analysis/` |
| RF input table (samples × features, with metadata columns) | `random-forest-analysis/` |

See individual script headers or the top-level `README.md` for detailed column/format requirements.
