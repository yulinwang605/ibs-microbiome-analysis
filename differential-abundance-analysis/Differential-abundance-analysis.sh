#1 LEfSe analsyis parameters
#For LEfSe analysis, taxa with a Log10 LDA score >2 were considered significant, while FDR <0.05 was considered significant

#2 #########MaAsLin2 analsysi used R script

# ==============================================================================
# MaAsLin2 Statistical Analysis Pipeline
# ==============================================================================
# This script performs multivariable association testing between microbial 
# profiles and metadata using MaAsLin2.

# Load required libraries
library(Maaslin2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. ENVIRONMENT SETTINGS & DATA LOADING
# ------------------------------------------------------------------------------
# Set the working directory to the path containing your data files
setwd("/path/to/taxonomic-profile")

# Load the microbial relative abundance dataset (rows as features, columns as samples)
abundance_data <- read.delim("abundance_table.txt", 
                             header = TRUE, sep = "\t", row.names = 1)

# Load the companion metadata file
metadata <- read.table("metadata_table.txt", 
                       header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# Check the structure of the metadata to ensure key study groups are present
head(metadata)

# Convert the grouping variables to categorical factors
metadata$Twogroups <- as.factor(metadata$Twogroups)
metadata$Threegroups <- as.factor(metadata$Threegroups)

# ------------------------------------------------------------------------------
# 2. RUN MAASLIN2 MODELS WITH DIFFERENT REFERENCE STAGES
# ------------------------------------------------------------------------------

# --- Run 1: Comparison using 'V2' as the reference baseline ---
results_ref_v2 <- Maaslin2(
  input_data = abundance_data, 
  input_metadata = metadata, 
  output = "Maaslin2_output_vs_V2",  # Directory for saving outputs
  fixed_effects = c("Stage", "NSR.score", "Response", "diarrhea_response", "stomach_ache_response"),
  random_effects = NULL,             # Specify random effects here if using a mixed-effects model
  normalization = "TSS",             # Total Sum Scaling normalization
  transform = "LOG",                 # Logarithmic transformation
  analysis_method = "LM",            # General linear modeling (LM)
  min_prevalence = 0.1,              # Filter out features present in less than 10% of samples
  reference = c("Stage,V2")          # Set V2 as the reference level
)

# --- Run 2: Comparison using 'V4' as the reference baseline ---
results_ref_v4 <- Maaslin2(
  input_data = abundance_data, 
  input_metadata = metadata, 
  output = "Maaslin2_output_vs_V4", 
  fixed_effects = c("Stage", "NSR.score", "Response", "diarrhea_response", "stomach_ache_response"),
  random_effects = NULL, 
  normalization = "TSS", 
  transform = "LOG", 
  analysis_method = "LM", 
  min_prevalence = 0.1, 
  reference = c("Stage,V4")          # Set V4 as the reference level
)

# --- Run 3: Comparison using 'V6' as the reference baseline ---
results_ref_v6 <- Maaslin2(
  input_data = abundance_data, 
  input_metadata = metadata, 
  output = "Maaslin2_output_vs_V6", 
  fixed_effects = c("Stage", "NSR.score", "Response", "diarrhea_response", "stomach_ache_response"),
  random_effects = NULL, 
  normalization = "TSS", 
  transform = "LOG", 
  analysis_method = "LM", 
  min_prevalence = 0.1, 
  reference = c("Stage,V6")          # Set V6 as the reference level
)


#3 #########ANCOM-BC2 analsysi used R script
library(mia)
library(ANCOMBC)

###################################ANCOM_BC2 for control group samples, samples with v2v4v6v7
setwd("/path/to/taxonomic-profile-count-based")

#using V2-vs-V4 as an example
# microbial count table

otu_mat = read.table("Counts-of-V2-vs-V4-control -samples.txt",sep = "\t",header = TRUE,
                     row.names = 1,check.names = FALSE) #counts-based abundance of V2 and V4 samples from control can be retrieved from the source data (Figure S5 related source data).

otu_mat = as.matrix(otu_mat)
assays = SimpleList(counts = otu_mat)

# sample metadata
smd = read.table("metadata-of-control-v2-vs-V4.txt",sep = "\t",header = TRUE,row.names = 1) #metadata of V2 and V4 samples from control can be retrieved from the source data (Tables 1-3 and Figure 1 related source data).


smd = DataFrame(smd)

# taxonomy table
tax_tab = read.table("tax-info.txt",sep = "\t",header = TRUE,row.names = 1) #tax-info can be retrieved from the source data (Figure S5 related source data)
tax_tab = DataFrame(tax_tab)

# create TSE
tse = TreeSummarizedExperiment(assays = assays,
                               colData = smd,
                               rowData = tax_tab)

# convert TSE to phyloseq
pseq = makePhyloseqFromTreeSummarizedExperiment(tse)

#####ancombc2
out = ancombc2(data = pseq, assay_name = "counts",
               tax_level = "Species", 
               p_adj_method = "BH", prv_cut = 0.10, fix_formula="Stage",
               group = "Stage", struc_zero = FALSE, neg_lb = FALSE,
               em_control = list(tol = 1e-5, max_iter = 100), 
               alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
res_prim = out$res
res_global = out$res_global
write.table(res_prim,file="ancombc2_V2-vs-V4.txt",sep = "\t",quote = FALSE,row.names = FALSE) 
#####ancombc2
out = ancombc2(data = pseq, assay_name = "counts",
               tax_level = "Species", 
               p_adj_method = "BH", prv_cut = 0.10, fix_formula="Stage+Age+Weight",
               group = "Stage", struc_zero = FALSE, neg_lb = FALSE,
               em_control = list(tol = 1e-5, max_iter = 100), 
               alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
res_prim = out$res
res_global = out$res_global
write.table(res_prim,file="ancombc2_V2-vs-V4-adjusted.txt",sep = "\t",quote = FALSE,row.names = FALSE)
