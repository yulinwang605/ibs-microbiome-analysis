###############################################################################
########## PCoA statistical analyses microbial abundance ######################
###############################################################################

library(ggplot2)
library(vegan)
library(BiodiversityR)
library(RColorBrewer)
library(ggthemes)
library(dplyr)

####################################################################################################
#PCoA of samples including v2 v4 v6 and v7
####################################################################################################
setwd("/path/to/abundance-and-metadata-files")
#use rel abu
#subtype.cell<-read.delim(file="metaphlan_rel_ab_species-remove-given-samples-treated.txt",header=T,sep="\t", row.names = 1,stringsAsFactors = F)
#use rel abu with unclassified
subtype.cell<-read.delim(file="Rel-abu-of-target_in_paper_samples.txt",header=T,sep="\t", row.names = 1,stringsAsFactors = F)

#subtype.cell<-subtype.cell[-1,]#remove first column
#filter relative abundance
table <- subtype.cell
table[table>0.001]<-1
table[table<1]<-0
table.generalist <- subtype.cell[,which(colSums(table)>=80)]
subtype.cell <- table.generalist

#subtype.cell<-table[,!colSums(table)==0]  ## 388
season<-read.table(file="metadata-for-beta-diversity_in_paper_samples.txt",header=T,sep="\t",stringsAsFactors = F)

scaleFUN1 <- function(x) sprintf("%.1f", x)
scaleFUN2 <- function(x) sprintf("%.2f", x)

###############################################################################
############################# 3D PCoA statistical analyses (Bray-curtis) ######
dist_subtype.cell<-vegdist(subtype.cell,method="bray") #Bray-Curtis distance calculation 
cmdscale_subtype.cell<-cmdscale(dist_subtype.cell, k=nrow(subtype.cell)-1,eig=T, add=F) #PCoA calculation
cmdscale_subtype.cell<- add.spec.scores(cmdscale_subtype.cell,subtype.cell,method="pcoa.scores", Rscale=T, scaling=1, multi=1) #Calculate scores (coordinates) for ARG-MRG types
tmp1<-paste("PC1 (",round(cmdscale_subtype.cell[[6]][1],2),"%)",sep="") 
tmp2<-paste("PC2 (",round(cmdscale_subtype.cell[[6]][2],2),"%)",sep="") 
tmp3<-paste("PC3 (",round(cmdscale_subtype.cell[[6]][3],2),"%)",sep="") 

#### prepare for plotting
pcoa_subtype.cell<-as.data.frame(cmdscale_subtype.cell[[1]][,1:3]) #extract data frame for PCoA plotting
colnames(pcoa_subtype.cell)<-c("PC1","PC2","PC3")
pcoa_subtype.cell$Sample<-rownames(pcoa_subtype.cell)

pcoa_season<-merge(pcoa_subtype.cell,season,by="Sample",all.x=T) #get season information of sample

#####set color
colourCount = length(unique(pcoa_season$Stage))
getPalette = colorRampPalette(brewer.pal(7, "Set1"))

###using ggplot##colored by taxonomy.
ggplot(pcoa_season, aes(x=PC1 , y=PC2,color = Stage,shape = Group)) + 
  geom_point(aes(size=NSR.score))+scale_color_manual(values = getPalette(colourCount))+
  scale_size_continuous(range = c(0.5,7))+
  scale_color_manual(values = c("V2" = "#C4BFBF", "V4" = "#4DAF4A", "V6" = "#4DA8F4", "V7" = "#F4442C"))+
  theme_calc()


#output pcoa matrix
write.table(pcoa_season,file = "PCoA-matrix.txt",sep = "\t",quote = FALSE)


