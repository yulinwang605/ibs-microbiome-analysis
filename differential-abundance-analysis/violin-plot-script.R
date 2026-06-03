# ==============================================================================
# 1. LOAD LIBRARIES
# ==============================================================================
library(ggstatsplot)
library(ggplot2)
library(palmerpenguins)
library(tidyverse)
library(data.table)
library(dplyr)
library(coin)       # For independent-samples Wilcoxon test
library(PMCMRplus)  # For post-hoc pairwise comparisons
library(FSA)        # For Dunn's test

# ==============================================================================
# 2. DEFINE REUSABLE STATISTICAL FUNCTION
# ==============================================================================
# Function to perform full statistical pipeline (Between-groups and Within-groups)
run_complete_statistical_analysis <- function(df, target_species) {
  cat("\n", paste(rep("#", 60), collapse = ""), "\n")
  cat("  STATISTICAL ANALYSIS FOR:", target_species, "\n")
  cat(paste(rep("#", 60), collapse = ""), "\n")
  
  # Ensure the target formula can handle dynamic column names
  formula_between <- as.formula(paste(target_species, "~ Group"))
  formula_within  <- as.formula(paste(target_species, "~ Stage"))
  
  # --- Part A: Between-Group Comparison (Placebo vs. Treatment at each Stage) ---
  group_comparison <- df %>%
    group_by(Stage) %>%
    do(w_test = wilcox.test(formula_between, data = ., exact = FALSE)) %>%
    summarise(Stage = Stage, p_value = w_test$p.value)
  
  # FDR adjustment for between-group comparison
  group_comparison$p_adj <- p.adjust(group_comparison$p_value, method = "fdr")
  
  print("--- Between-Group Results (Placebo vs. Treatment at each Stage) ---")
  print(group_comparison)
  
  # --- Part B: Within-Group Pairwise Stage Comparison (Dunn's Test) ---
  analyze_pairwise_stages <- function(group_name, sub_df) {
    sub_data <- sub_df %>% filter(Group == group_name)
    
    cat("\n", paste(rep("=", 20), collapse = ""), " Analyzing Group:", group_name, paste(rep("=", 20), collapse = ""), "\n")
    
    # Global Kruskal-Wallis Test
    kw_test <- kruskal.test(formula_within, data = sub_data)
    cat("Kruskal-Wallis Global P-value:", round(kw_test$p.value, 5), "\n")
    
    # Post-hoc Pairwise Comparisons using Dunn's Test (with Benjamini-Hochberg adjustment)
    dunn_res <- dunnTest(formula_within, data = sub_data, method = "bh")
    pairwise_results <- dunn_res$res
    
    # Filter for significant pairs (P.adj < 0.05)
    sig_pairs <- pairwise_results %>% filter(P.adj < 0.05)
    
    return(list(all_pairs = pairwise_results, sig_pairs = sig_pairs))
  }
  
  placebo_res   <- analyze_pairwise_stages("Placebo", df)
  treatment_res <- analyze_pairwise_stages("Treatment", df)
  
  cat("\n--- Placebo Group: All Pairwise Stage Comparisons ---\n")
  print(placebo_res$all_pairs)
  
  cat("\n--- Treatment Group: All Pairwise Stage Comparisons ---\n")
  print(treatment_res$all_pairs)
}


# ==============================================================================
# 3. ANALYSIS FOR: Agathobacter rectalis abundace was used as an example
# ==============================================================================

# Data loading and preprocessing
data_agathobacter <- fread("Agathobacter_rectalis.txt", header = TRUE, sep = "\t")
data_agathobacter$Stage <- factor(data_agathobacter$Stage, levels = c("V2", "V4", "V6", "V7"))
data_agathobacter$Group <- as.factor(data_agathobacter$Group)

# Log10 transformation for visualization
data_agathobacter$Agathobacter_rectalis_log10 <- log10(data_agathobacter$Agathobacter_rectalis + 1)  

# Violin plot using ggbetweenstats
grouped_ggbetweenstats(
  data = data_agathobacter,                         
  x = Stage,                           
  y = Agathobacter_rectalis_log10, 
  grouping.var = Group,               
  plot.type = "violin",               
  pairwise.comparisons = TRUE,        
  pairwise.display = "significant",   
  xlab = "Stage",                     
  ylab = "log10(Abundance of Agathobacter rectalis)", 
  ggtheme = ggplot2::theme_minimal(), 
  results.subtitle = TRUE,            
  centrality.point.args = list(size = 2, color = "darkred"),
  ggplot.component = scale_y_continuous(limits = c(0, 1.7))
)

# Execute statistical pipeline
run_complete_statistical_analysis(data_agathobacter, "Agathobacter_rectalis")