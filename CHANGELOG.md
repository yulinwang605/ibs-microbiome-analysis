# Changelog

All notable changes to this repository will be documented in this file.

## [1.0.0] - 2026-06-03

### Added
- Initial public release accompanying manuscript submission
- `sequence-processing/sequencing_data_processing_commands.sh`: BBDuk QC, BBSplit host removal, MetaPhlAn4 taxonomic profiling pipeline
- `diversity-analysis/alpha/alpha-diversity-analysis.R`: Shannon and Simpson alpha diversity calculation and visualisation
- `diversity-analysis/beta/pcoa_analysis.R`: Bray-Curtis PCoA computation and plotting
- `diversity-analysis/beta/retrieve-given-stages/plot-pcoA-matrix.R`: Per-stage PCoA panel visualisation
- `differential-abundance-analysis/Differential-abundance-analysis.sh`: LEfSe parameters, MaAsLin2, and ANCOM-BC2 pipelines
- `differential-abundance-analysis/violin-plot-script.R`: Violin plots with Wilcoxon/Kruskal-Wallis/Dunn's test annotations
- `microbial-association-network-analysis/Network-analysis-using-SPIEC-EASI.R`: SPIEC-EASI co-occurrence network construction and Cytoscape export
- `random-forest-analysis/RF_clr-scale_leave-one-out.py`: CLR-scaled Random Forest with LOOCV and ROC/confusion matrix outputs
