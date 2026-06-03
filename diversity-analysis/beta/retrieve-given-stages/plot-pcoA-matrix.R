library("ggplot2")
library("RColorBrewer")
library(ggthemes)

#all-samples
pcoa_all <-read.table("PCoA-matrix.txt",header = TRUE, sep = "\t")
#control-v2
pcoa_control_v2 <- read.table("control-v2.txt",header = TRUE, sep = "\t")
#control-v4
pcoa_control_v4 <- read.table("control-v4.txt", header = TRUE, sep = "\t")
#control-v6
pcoa_control_v6 <- read.table("control-v6.txt", header = TRUE, sep = "\t")
#control-v7
pcoa_control_v7 <- read.table("control-v7.txt", header = TRUE, sep = "\t")

#treatment-v2
pcoa_treatment_v2 <- read.table("treatment-v2.txt",header = TRUE, sep = "\t")
#treatment-v4
pcoa_treatment_v4 <- read.table("treatment-v4.txt", header = TRUE, sep = "\t")
#treatment-v6
pcoa_treatment_v6 <- read.table("treatment-v6.txt", header = TRUE, sep = "\t")
#treatment-v7
pcoa_treatment_v7 <- read.table("treatment-v7.txt", header = TRUE, sep = "\t")

# Determine the global range for NRS.score across both datasets
nrs_min <- min(c(pcoa_all$NRS.score, pcoa_control_v2$NRS.score, pcoa_control_v4$NRS.score, pcoa_control_v6$NRS.score, pcoa_control_v7$NRS.score,
                 pcoa_treatment_v2$NRS.score,pcoa_treatment_v4$NRS.score,pcoa_treatment_v6$NRS.score,pcoa_treatment_v7$NRS.score), na.rm = TRUE)
nrs_max <- max(c(pcoa_all$NRS.score, pcoa_control_v2$NRS.score, pcoa_control_v4$NRS.score, pcoa_control_v6$NRS.score, pcoa_control_v7$NRS.score,
                 pcoa_treatment_v2$NRS.score,pcoa_treatment_v4$NRS.score,pcoa_treatment_v6$NRS.score,pcoa_treatment_v7$NRS.score), na.rm = TRUE)

# Set color for the Stage
colourCount <- length(unique(pcoa_all$Stage))
getPalette <- colorRampPalette(brewer.pal(7, "Set1"))

# Plot for all-samples data
ggplot(pcoa_all, aes(x = PC1, y = PC2, color = Stage, shape = Group)) + 
  geom_point(aes(size = NRS.score)) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - All Samples")

# Plot for control-v2 data
ggplot(pcoa_control_v2, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 16) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Placebo V2 Samples")


# Plot for control-v4 data
ggplot(pcoa_control_v4, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 16) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Placebo V4 Samples")


# Plot for control-v6 data
ggplot(pcoa_control_v6, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 16) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Placebo V6 Samples")



# Plot for control-v7 data
ggplot(pcoa_control_v7, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 16) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Placebo V7 Samples")



# Plot for treatment-v2 data
ggplot(pcoa_treatment_v2, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 17) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Treatment V2 Samples")



# Plot for treatment-v4 data
ggplot(pcoa_treatment_v4, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 17) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Treatment V4 Samples")


# Plot for treatment-v6 data
ggplot(pcoa_treatment_v6, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 17) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Treatment V6 Samples")


# Plot for treatment-v7 data
ggplot(pcoa_treatment_v7, aes(x = PC1, y = PC2, color = Stage)) + 
  geom_point(aes(size = NRS.score), shape = 17) + 
  scale_size_continuous(limits = c(nrs_min, nrs_max), range = c(1, 5)) + 
  scale_color_manual(values = getPalette(colourCount)) + 
  theme_calc() +
  ggtitle("PCoA Plot - Treatment V7 Samples")
