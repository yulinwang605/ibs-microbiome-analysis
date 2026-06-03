# Install and load the required packages
install.packages("vegan")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("tidyr")

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(ggpubr)

# Load input data
# Read the microbial relative abundance data (tab-separated)
microbiome_data <- read.table("Rel-abu-of-target-samples_in_paper_samples.txt", header = TRUE, row.names = 1, sep = "\t")
# Read the sample metadata (tab-separated)
metadata <- read.table("metadata-of-target-samples_in_paper_samples.txt", header = TRUE, row.names = 1, sep = "\t")
metadata$Sample <- rownames(metadata)

# Filter out low-abundance rare taxa
table <- microbiome_data
table[table>0.001] <- 1
table[table<0.001] <- 0
table.generalist <- microbiome_data[, which(colSums(table)>=50)]
microbiome_data <- table.generalist


# Calculate Alpha Diversity Indices
# 1. Shannon diversity index
shannon_index <- diversity(microbiome_data, index = "shannon")

# 2. Simpson diversity index
simpson_index <- diversity(microbiome_data, index = "simpson")

# 3. Chao1 index
#chao1_index <- estimateR(microbiome_data)["Chao1", ]

# 4. ACE index
#ace_index <- estimateR(microbiome_data)["ACE", ]

# Combine the calculated diversity indices into a single data frame
diversity_df <- data.frame(
  Sample = rownames(metadata),
  Shannon = shannon_index,
  Simpson = simpson_index
)

# Merge metadata with alpha diversity indices
diversity_df <- diversity_df %>%
  left_join(metadata, by = "Sample")
write.table(diversity_df, file = "Shannon-simpson-index.tsv", quote = FALSE, sep="\t")


# Grouped boxplot for Shannon index
ggplot(diversity_df, aes(x=Stage, y=Shannon, fill=Group)) + 
  geom_boxplot()+
  #geom_jitter(color="gray32",size=2,alpha=0.5)+
  geom_signif(comparisons = list(c("Treated","Control")), map_signif_level = TRUE, step_increase = 0.1)+
  stat_compare_means(method = "wilcox.test")+
  theme_calc()

ggplot(diversity_df, aes(x=Group, y=Shannon, fill=Stage)) + 
  geom_boxplot()+
  #geom_jitter(color="gray32",size=2,alpha=0.5)+
  #geom_signif(comparisons = list(c("Treated","Control")), map_signif_level = TRUE, step_increase = 0.1)+
  #stat_compare_means(aes(group = Stage), method = "wilcox.test")+
  theme_calc()