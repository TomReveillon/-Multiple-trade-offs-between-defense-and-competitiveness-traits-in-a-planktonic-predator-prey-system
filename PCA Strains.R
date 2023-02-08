setwd("~/LIMNO 2019-2022/Experiments/Prey Phenotyping")

rm(list=ls())

library(ade4)
library(caret)
library(corrplot)
library(cowplot)
library(data.table)
library(deSolve)
library(dplyr)
library(dunn.test)
library(factoextra)
library(foreach)
library(ggforce)
library(ggplot2)
library(ggrepel)
library(gmodels)
library(grid)
library(gridExtra)
library(knitr)
library(lme4)
library(lmtest)
library(magrittr)
library(nlme)
library(plotly)
library(plyr)
library(randomForest)
library(propagate)
library(reshape2)
library(scales)

###################################################################################
###################################################################################
##### ANALYSIS OF MORPHOLOGICAL FEATURES IN STRAINS WITH DIFFERENT MORPHOTYPES ####
###################################################################################
###################################################################################

# Load the datasets
Folder="Collection Strains"
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
DataInter=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")

# Verify dataset classes
lapply(DataInter, class)

# Mutate dataset classes to numeric
DataInter %<>% mutate_if(is.factor,as.numeric); DataInter %<>% mutate_if(is.integer,as.numeric)

# Create strain identities
IDCR1=c(rep("CR1",5000))
IDCR2=c(rep("CR2",5000))
IDCR3=c(rep("CR3",5000))
IDCR4=c(rep("CR4",5000))
IDCR5=c(rep("CR5",5000))
IDCR6=c(rep("CR6",5000))

# Add strain identities and change column names
DataInter$Strain=c(IDCR1,IDCR2,IDCR3,IDCR4,IDCR5,IDCR6)
colnames(DataInter)=gsub("\\_", "\\.", colnames(DataInter))

# Create strain identity vector
DataInter$Strain=factor(DataInter$Strain)
Strain=DataInter$Strain


#####################################
#### Principal component analysis ###
#####################################

# Preserve the order of strains
DataInter$Strain=factor(DataInter$Strain, levels=unique(DataInter$Strain))

# Remove strain identities and rename columns
DataPCA=subset(DataInter, select=-c(Strain))
DataPCA=DataPCA %>% select_if(grepl("M05", names(.)))
colnames(DataPCA)=gsub("\\.M05", "", colnames(DataPCA))
colnames(DataPCA)=gsub("\\.", " ", colnames(DataPCA))

# Calculate the feature volume 
DataPCA$Volume=(4*3.1416*(DataPCA$Diameter/2)^3)/3

# Select columns
DataPCA=sapply(DataPCA[,c(1:40)], as.numeric)
DataPCA=as.data.frame(DataPCA)

# Principal component analysis
PCA=dudi.pca(DataPCA, scannf=F, nf=2)

# Create a dataset of features
DataFeat=data.frame(Feature=colnames(DataPCA),Category=colnames(DataPCA))

# Create categories of features
DataFeat$Category %<>% 
  gsub(c("^Area$|^Diameter$|^Height$|^Length$|^Major Axis$|^Maximum Thickness$|^Minimum Thickness$|^Minor Axis$|^Spot Area$|^Perimeter$|^Width$|^Volume$"), "Size", .) %>%
  gsub(c("^Aspect Ratio$|^Circularity$|^Compactness$|^Elongatedness$|^Lobe Count$|^Shape Ratio$|^Symmetry 2$|^Symmetry 3$|^Symmetry 4$"), "Shape", .) %>%
  gsub(c("^Aspect Ratio Intensity$|^Bright Detail Intensity$|^Contrast$|^Intensity$|^Major Axis Intensity$|^Maximum Pixel$|^Maximum Spot Intensity$|^Mean Pixel$|^Median Pixel$|^Minimum Pixel$|^Minimum Spot Intensity$|^Minor Axis Intensity$|^Modulation$|^Raw Intensity$|^Raw Maximum Pixel$|^Raw Mean Pixel$|^Raw Median Pixel$|^Raw Minimum Pixel$|^Shape Intensity$|^Spot Count$"), "Signal", .)

# Extract features coordinates
DataFeat$Coord1=PCA$co[,1]
DataFeat$Coord2=PCA$co[,2]
DataFeat$Contr1=PCA$c1[,1]
DataFeat$Contr2=PCA$c1[,2]

# Subset per category of features
DataFeat1=subset(DataFeat, Category=="Size")
DataFeat1[,c(5,6)]=abs(DataFeat1[,c(5,6)])
head(DataFeat1[rev(order(DataFeat1$Contr1)),],5)
head(DataFeat1[rev(order(DataFeat1$Contr2)),],5)
DataFeat1=DataFeat1[c(1,6,12),]

DataFeat2=subset(DataFeat, Category=="Shape")
DataFeat2[,c(5,6)]=abs(DataFeat2[,c(5,6)])
head(DataFeat2[rev(order(DataFeat2$Contr1)),],5)
head(DataFeat2[rev(order(DataFeat2$Contr2)),],5)
DataFeat2=DataFeat2[c(2,6,9),]

DataFeat3=subset(DataFeat, Category=="Signal")
DataFeat3[,c(5,6)]=abs(DataFeat3[,c(5,6)])
head(DataFeat3[rev(order(DataFeat3$Contr1)),],5)
head(DataFeat3[rev(order(DataFeat3$Contr2)),],5)
DataFeat3=DataFeat3[c(1,4,7),]

# Create a dataset and rename features
DataFeat=rbind(DataFeat1,DataFeat2,DataFeat3)
DataFeat$Feature=c("Area","Thickness","Volume","Circularity","Sphericity","Asymmetry","Roundness","Chlorophyll","Aggregation")
DataFeat$Trait=c("Size","Size","Size","Shape","Shape","Shape","Shape","Size","Size")

# Position of features
tiff('PCA Features.tiff', units="in", width=8, height=8, res=1000)
ggplot(DataFeat, aes(Coord1, Coord2, group=Trait)) +
  geom_hline(yintercept=0, color="grey50", linetype="dotted", size=1) + 
  geom_vline(xintercept=0, color="grey50", linetype="dotted", size=1) +
  geom_segment(aes(x=0, y=0, xend=Coord1, yend=Coord2), color="black", linetype="solid", size=0.7) +
  geom_point(aes(color=Trait), pch=16, size=3) +
  geom_label_repel(aes(label=Feature), fill="white", color="black", segment.color=NA, size=6) +
  ylab(expression('Dimension 2 (19.9 %)')) + xlab(expression('Dimension 1 (46.5 %)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-1.0,1.0,by=0.5), limits=c(-1.1,1.1)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-1.0,1.0,by=0.5), limits=c(-1.1,1.1)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("Size"="firebrick3","Shape"="royalblue3")) +
  scale_color_manual(values=c("Size"="firebrick3","Shape"="royalblue3")) +
  theme(legend.position="none")
dev.off()


#############################
#### Interstrain features ###
#############################

# Add identities to the dataset
DataInter=DataPCA
DataInter$Strain=Strain

# Extract individuals coordinates
DataInter$Coord1=PCA$li[,1]
DataInter$Coord2=PCA$li[,2]
DataInter$Contr1=PCA$l1[,1]
DataInter$Contr2=PCA$l1[,2]

# Select representative individuals 
DataInter$Contr=(DataInter$Contr1+DataInter$Contr2)/2

# Split the dataset
SplitIndiv=split(DataInter, list(DataInter$Strain))

# Position of individuals
PlotFunc=function(x) {
  ggplot(x, aes(Coord1, Coord2, group=Strain)) +
    geom_hline(yintercept=0, color="grey50", linetype="dotted", size=1) + 
    geom_vline(xintercept=0, color="grey50", linetype="dotted", size=1) +
    geom_point(aes(color=Strain), size=0.7, pch=16, alpha=0.5) +
    stat_ellipse(aes(Coord1, Coord2, fill=Strain, color=Strain), level=0.90, alpha=0.5, geom="polygon") +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-10,10,by=5), limits=c(-11,11)) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-30,30,by=15), limits=c(-30,30)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Strain, ncol=3, nrow=2) +
    theme(legend.position="none")
}

tiff('PCA Individuals Interstrain.tiff', units="in", width=15, height=8, res=1000)
Panel1=lapply(SplitIndiv, PlotFunc)
Yaxis=textGrob(expression('Dimension 2 (19.9 %)'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Dimension 1 (46.5 %)'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(grobs=Panel1, left=Yaxis, bottom=Xaxis, ncol=3, nrow=2)
dev.off()

# Select relevant features
FeatInter=data.frame(DataPCA[,c(1,5)],Strain)
FeatInter[,1]=FeatInter[,1]/10^3
FeatInter[,2]=FeatInter[,2]/10^1
colnames(FeatInter)=c("Area","Circularity","Strain")
  
# Sort and split the dataset
MeltInter=melt(FeatInter, id.vars=c("Strain"))
colnames(MeltInter)[2:3]=c("Feature","Value")

# Split the dataset
SplitInter=split(MeltInter, list(MeltInter$Feature))

# Violin plots of features
PlotFunc=function(x) {
  ggplot(x, aes(Feature, Value, group=Strain)) +
    geom_violin(aes(fill=Strain, color=Strain)) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_blank()) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
    theme(axis.ticks.x=element_blank()) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.5)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    facet_wrap(~Strain, ncol=3, nrow=2) +
    theme(legend.position="none")
}

tiff('PCA Features Interstrain 1.tiff', units="in", width=15, height=8, res=1000)
Panel2=lapply(SplitInter, PlotFunc)
Panel2[[1]]=Panel2[[1]] + scale_y_continuous(expression('Cell area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,3,by=1), limits=c(0,3))
Panel2[[2]]=Panel2[[2]] + scale_y_continuous(expression('Cell circularity'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,6,by=2), limits=c(0,6))
Panel2[[1]]=Panel2[[1]] + scale_x_discrete(expression('Strain'))
Panel2[[2]]=Panel2[[2]] + scale_x_discrete(expression('Strain'))
grid.arrange(grobs=Panel2, ncol=2, nrow=1)
dev.off()

# Density plots of features
PlotFunc=function(x) {
  ggplot(x, aes(Value, group=Strain)) +
    geom_density(aes(fill=Strain, color=Strain)) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.3)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    theme(legend.position="none")
}

tiff('PCA Features Interstrain 2.tiff', units="in", width=15, height=8, res=1000)
Panel3=lapply(SplitInter, PlotFunc)
Panel3[[1]]=Panel3[[1]] + scale_y_continuous(expression('Density'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,25,by=5), limits=c(0,25))
Panel3[[1]]=Panel3[[1]] + scale_x_continuous(expression('Cell area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,1,by=0.2), limits=c(0,1))
Panel3[[2]]=Panel3[[2]] + scale_y_continuous(expression('Density'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,1,by=0.2), limits=c(0,1))
Panel3[[2]]=Panel3[[2]] + scale_x_continuous(expression('Cell circularity'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,4,by=1), limits=c(0,4))
grid.arrange(grobs=Panel3, ncol=2, nrow=1)
dev.off()


################################################
### Morphotypes identification with features ###
################################################

# Add identities to the dataset
DataMorpho=DataPCA
DataMorpho$Strain=Strain

# Calculate confidence intervals
AreaCI=ci(subset(DataMorpho, Strain=="CR4")$Area, confidence=0.975)

# Calculate mean area of single cells
AreaM=mean(subset(DataMorpho, !Area < 107.49 & !Area > 108.59)$Area)

# Distribution of area of single cells
ggplot(subset(DataMorpho, Strain=="CR4"), aes(Area)) +
  geom_density(fill="gold2", color="gold2", linetype="solid", alpha=0.3,size=1) +
  geom_vline(xintercept=AreaM, color="gold2", linetype="dashed", size=1) +
  ylab(expression('Frequency')) + xlab(expression('Cell area'~'('*µm^2*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
  scale_y_continuous(labels=function(x) sprintf("%.3f", x), breaks=seq(0,0.020,by=0.005), limits=c(0,0.021)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,600,by=200), limits=c(0,600)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  theme(legend.position="none")

# Calculate cell numbers per image
DataMorpho$Cells=round(DataMorpho$Area/AreaM,0)
DataMorpho$Cells[DataMorpho$Cells == 0]=1

# Assign clump categories per image
Morpho=data.frame(DataMorpho %>% mutate(Clumps=case_when(Cells <= 2 ~ "S", Cells <= 6 ~ "SC", Cells <= 10 ~ "MC", Cells > 10 ~ "LC")))[,43]


#############################
#### Intrastrain features ###
#############################

# Add identities to the dataset
DataIntra=DataPCA
DataIntra$Strain=Strain
DataIntra$Morpho=Morpho

# Extract individuals coordinates
DataIntra$Coord1=PCA$li[,1]
DataIntra$Coord2=PCA$li[,2]
DataIntra$Contr1=PCA$l1[,1]
DataIntra$Contr2=PCA$l1[,2]

# Split the dataset
SplitIndiv=split(DataIntra, list(DataIntra$Strain))

# Position of individuals
PlotFunc=function(x) {
  ggplot(x, aes(Coord1, Coord2, group=interaction(Strain,Morpho))) +
    geom_hline(yintercept=0, color="grey50", linetype="dotted", size=1) + 
    geom_vline(xintercept=0, color="grey50", linetype="dotted", size=1) +
    geom_point(aes(color=Strain), size=0.7, pch=16, alpha=0.5) +
    stat_ellipse(aes(Coord1, Coord2, fill=Strain, color=Strain), level=0.90, alpha=0.3, geom="polygon") +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_blank()) + theme(axis.title.x=element_blank()) +
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-10,10,by=5), limits=c(-11,11)) +
    scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-30,30,by=15), limits=c(-30,30)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Strain, ncol=3, nrow=2) +
    theme(legend.position="none")
}

tiff('PCA Individuals Intrastrain.tiff', units="in", width=15, height=8, res=1000)
Panel4=lapply(SplitIndiv, PlotFunc)
Yaxis=textGrob(expression('Dimension 2 (19.9 %)'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Dimension 1 (46.5 %)'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(grobs=Panel4, left=Yaxis, bottom=Xaxis, ncol=3, nrow=2)
dev.off()

# Select relevant features
FeatIntra=data.frame(DataPCA[,c(1,5)],Strain,Morpho)
FeatIntra[,1]=FeatIntra[,1]/10^3
FeatIntra[,2]=FeatIntra[,2]/10^1
colnames(FeatIntra)=c("Area","Circularity","Strain","Morpho")

# Sort and split the dataset
MeltIntra=melt(FeatIntra, id.vars=c("Strain","Morpho"))
colnames(MeltIntra)[3:4]=c("Feature","Value")
MeltIntra$Morpho=factor(MeltIntra$Morpho, levels=c("S","SC","MC","LC"))

# Create codes binding strains and morphotypes
MeltIntra$Code=factor(paste(MeltIntra$Strain, MeltIntra$Morpho, sep=""))

# Replace some morphotype categories
MeltIntra$Code=gsub("CR2LC", "CR2MC", MeltIntra$Code) 
MeltIntra$Code=gsub("CR3LC", "CR3MC", MeltIntra$Code) 
MeltIntra$Code=gsub("CR4MC", "CR4SC", MeltIntra$Code) 
MeltIntra$Code=gsub("CR4LC", "CR4SC", MeltIntra$Code) 
MeltIntra$Code=gsub("CR5MC", "CR5SC", MeltIntra$Code)
MeltIntra$Code=gsub("CR5LC", "CR5SC", MeltIntra$Code)
MeltIntra$Code=gsub("CR6MC", "CR6SC", MeltIntra$Code)
MeltIntra$Code=gsub("CR6LC", "CR6SC", MeltIntra$Code)

MeltIntra$Morpho[MeltIntra$Code %in% c("CR2MC")]="MC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR4SC")]="SC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR5SC")]="SC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR6SC")]="SC"
MeltIntra$Code=factor(MeltIntra$Code, levels=c("CR1S","CR1SC","CR1MC","CR1LC","CR2S","CR2SC","CR2MC","CR3S","CR3SC","CR3MC","CR4S","CR4SC","CR5S","CR5SC","CR6S","CR6SC"))

# Split the dataset
SplitIntra=split(MeltIntra, list(MeltIntra$Feature))

# Violin plots of features
PlotFunc=function(x) {
  ggplot(x, aes(Feature, Value, group=Code)) +
    geom_violin(aes(fill=Strain, color=Strain)) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_blank()) +
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.ticks.x=element_blank()) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.5)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    facet_wrap(~Code, ncol=4, nrow=4) +
    theme(legend.position="none")
}

tiff('PCA Features Intrastrain 1.tiff', units="in", width=15, height=8, res=1000)
Panel5=lapply(SplitIntra, PlotFunc)
Panel5[[1]]=Panel5[[1]] + scale_y_continuous(expression('Cell area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,3,by=1), limits=c(0,3))
Panel5[[2]]=Panel5[[2]] + scale_y_continuous(expression('Cell circularity'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,6,by=2), limits=c(0,6))
Panel5[[1]]=Panel5[[1]] + scale_x_discrete(expression('Morphotype'))
Panel5[[2]]=Panel5[[2]] + scale_x_discrete(expression('Morphotype'))
grid.arrange(grobs=Panel5, ncol=2, nrow=1)
dev.off()

# Density plots of features
PlotFunc=function(x) {
  ggplot(x, aes(Value, group=Morpho)) +
    geom_density(aes(fill=Strain, color=Strain)) + ylab("Density") +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.3)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    facet_wrap(~Strain, ncol=3, nrow=2) +
    theme(legend.position="none")
}

tiff('PCA Features Intrastrain 2.tiff', units="in", width=15, height=8, res=1000)
Panel6=lapply(SplitIntra, PlotFunc)
Panel6[[1]]=Panel6[[1]] + scale_y_continuous(expression('Density'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,25,by=5), limits=c(0,25))
Panel6[[1]]=Panel6[[1]] + scale_x_continuous(expression('Cell area'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,1,by=0.2), limits=c(0,1))
Panel6[[2]]=Panel6[[2]] + scale_y_continuous(expression('Density'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,4,by=1), limits=c(0,4))
Panel6[[2]]=Panel6[[2]] + scale_x_continuous(expression('Cell circularity'), labels=function(x) sprintf("%.0f", x), breaks=seq(0,4,by=1), limits=c(0,4))
grid.arrange(grobs=Panel6, ncol=2, nrow=1)
dev.off()


#############################################
### Coefficients of variation in features ###
#############################################

# Mutate dataset classes to numeric
DataCV=as.data.frame(sapply(DataPCA, as.numeric))
Feature=colnames(DataCV)
DataCV$Strain=Strain
DataCV$Morpho=Morpho
  
# Calculate interstrain coefficients
SplitCV=split(DataCV, list(DataCV$Strain))
SplitCV=SplitCV[c("CR1","CR2","CR3","CR4","CR5","CR6")]
CVInter=abs(unlist(lapply(SplitCV, function(x) {apply(x[,c(1:40)], 2, function(x) sd(x)/mean(x)*100)})))
CVInter=data.frame(CVInter)
CVInter$Feature=rep(Feature,6)
rownames(CVInter)=c()

# Create a dataset per category
FeatCV=data.frame(CV=CVInter[,1])
FeatCV$Feature=rep(Feature,6)
FeatCV$Strain=rep(c("CR1","CR2","CR3","CR4","CR5","CR6"), each=40)

# Statistical tests
shapiro.test(FeatCV$CV)
bartlett.test(FeatCV$CV~FeatCV$Strain)
kruskal.test(FeatCV$CV, FeatCV$Strain)
dunn.test(FeatCV$CV, FeatCV$Strain)

# Plot the coefficients of variation
tiff('CV Features Interstrain.tiff', units="in", width=8, height=8, res=1000)
ggplot(FeatCV, aes(Strain, CV, group=Strain)) +
  geom_boxplot(aes(fill=Strain, color=Strain), size=0.5, width=0.7, outlier.shape=NA) + 
  geom_point(aes(color=Strain), size=1.5, pch=16, alpha=0.7, position=position_jitter(0.3)) +
  ylab("Coefficient of variation (%)") + xlab("Strain") +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_blank()) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.ticks.x=element_blank()) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,120,by=40), limits=c(0,125)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.3)) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()

# Calculate intrastrain coefficients
SplitCV2=split(DataCV, list(DataCV$Strain,DataCV$Morpho))
SplitCV2=SplitCV2[c("CR1.S","CR1.SC","CR1.MC","CR1.LC","CR2.S","CR2.SC","CR2.MC","CR3.S","CR3.SC","CR3.MC","CR4.S","CR4.SC","CR5.S","CR5.SC","CR6.S","CR6.SC")]
CVIntra=abs(unlist(lapply(SplitCV2, function(x) {apply(x[,c(1:40)], 2, function(x) sd(x)/mean(x)*100)})))
CVIntra=data.frame(CVIntra)
CVIntra$Feature=rep(Feature,8)
rownames(CVIntra)=c()

# Create a dataset per category
FeatCV2=data.frame(CV=CVIntra[,1])
FeatCV2$Variable=rep("CV",8*40)
FeatCV2$Feature=rep(Feature,8)
FeatCV2$Strain=rep(c(rep("CR1",4),rep("CR2",3),rep("CR3",3),rep("CR4",2),rep("CR5",2),rep("CR6",2)), each=40)
FeatCV2$Morpho=rep(c("CR1S","CR1SC","CR1MC","CR1LC","CR2S","CR2SC","CR2MC","CR3S","CR3SC","CR3MC","CR4S","CR4SC","CR5S","CR5SC","CR6S","CR6SC"), each=40)

# Statistical test
SplitFeatCV2=split(FeatCV2, list(FeatCV2$Strain))
lapply(SplitFeatCV2, function(x) shapiro.test(x$CV))
lapply(SplitFeatCV2, function(x) bartlett.test(x$CV~x$Morpho))
lapply(SplitFeatCV2, function(x) kruskal.test(x$CV, x$Morpho))
lapply(SplitFeatCV2, function(x) dunn.test(x$CV, x$Morpho))

# Preserve the order of variables
FeatCV2$Strain=factor(FeatCV2$Strain, levels=unique(FeatCV2$Strain))
FeatCV2$Morpho=factor(FeatCV2$Morpho, levels=unique(FeatCV2$Morpho))

# Split the dataset
SplitFeatCV=list(FeatCV2)

# Plot the coefficients of variation
PlotFunc=function(x) {
  ggplot(x, aes(Variable, CV, group=Strain)) +
    geom_boxplot(aes(fill=Strain, color=Strain), size=0.5, width=0.5, outlier.shape=NA) + 
    geom_point(aes(color=Strain), size=1.5, pch=16, alpha=0.7, position=position_jitter(0.2)) +
    ylab("Coefficient of variation (%)") + xlab("Morphotype") +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_blank()) +
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.ticks.x=element_blank()) +
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,120,by=40), limits=c(0,125)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.5)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
    facet_wrap(~Morpho, ncol=4, nrow=4) +
    theme(legend.position="none")
}

tiff('CV Features Intrastrain.tiff', units="in", width=15, height=8, res=1000)
Panel7=lapply(SplitFeatCV, PlotFunc)
grid.arrange(grobs=Panel7, ncol=1, nrow=1)
dev.off()


############################################
### Interstrain allometric relationships ###
############################################

# Correlation matrix for features
CorPCA=cor(DataPCA[,c(1:40)], method="spearman")

# Correlation plot for features
corrplot(cor(DataPCA[,c(1:40)]), method="circle", order="FPC", type="lower", diag=F, 
col=colorRampPalette(c("darkred","white","darkgreen"))(100), addgrid.col="grey90",
cl.lim=c(-1.0, 1.0), cl.pos="b", cl.cex=1.0, tl.col="black", tl.cex=1.0, tl.srt=45)

# Regression test
RegFunc=function(x) {summary(lm(x$Area~x$Circularity))}
OutReg=lapply(split(FeatIntra, list(FeatIntra$Strain)), RegFunc)

# Correlation test
CorFunc=function(x) {cor.test(x$Area, x$Circularity, method="pearson", exact=F)}
OutCor=lapply(split(FeatIntra, list(FeatIntra$Strain)), CorFunc)

# Create a labeling dataset
R2=c(as.character(expression(italic(R^2) *~'= -0.66')),as.character(expression(italic(R^2) *~'= -0.49')),as.character(expression(italic(R^2) *~'= -0.35')),as.character(expression(italic(R^2) *~'= 0.22')),as.character(expression(italic(R^2) *~'= 0.39')),as.character(expression(italic(R^2) *~'= 0.03')))
P=c(as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'= 0.521')))
DataL=data.frame(Strain=c("CR1","CR2","CR3","CR4","CR5","CR6"), R2, P, stringsAsFactors=F)

# Replace some morphotype categories
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR2" & FeatIntra$Morpho=="LC", "MC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR3" & FeatIntra$Morpho=="LC", "MC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR4" & FeatIntra$Morpho=="MC", "SC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR4" & FeatIntra$Morpho=="LC", "SC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR5" & FeatIntra$Morpho=="MC", "SC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR5" & FeatIntra$Morpho=="LC", "SC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR6" & FeatIntra$Morpho=="MC", "SC", FeatIntra$Morpho)
FeatIntra$Morpho=ifelse(FeatIntra$Strain=="CR6" & FeatIntra$Morpho=="LC", "SC", FeatIntra$Morpho)

tiff('Correlation Features Interstrain.tiff', units="in", width=15, height=8, res=1000)
ggplot(FeatIntra, aes(Circularity, Area, group=Strain)) + 
  geom_point(aes(color=Strain), size=0.3, pch=16, alpha=0.5) + coord_cartesian(ylim=c(0,1.5), xlim=c(0,5)) +
  geom_smooth(method='lm', formula=y~x, color="black", linetype="11", size=0.8, se=F, fullrange=T) +
  stat_ellipse(aes(Circularity, Area, group=Strain, fill=Strain, color=Strain), level=0.90, alpha=0.5, geom="polygon") +
  ylab(expression('Cell area')) + xlab(expression("Cell circularity")) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.5,by=0.5), limits=c(-0.1,1.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,5.0,by=1.0), limits=c(-2.5,5.0)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  facet_wrap(~Strain, ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()

tiff('Correlation Features Intrastrain.tiff', units="in", width=15, height=8, res=1000)
ggplot(FeatIntra, aes(Circularity, Area, group=Strain)) + 
  geom_point(aes(color=Strain), size=0.3, pch=16, alpha=0.5) + coord_cartesian(ylim=c(0,1.5), xlim=c(0,5)) +
  geom_smooth(method='lm', formula=y~x, color="black", linetype="11", size=0.8, se=F, fullrange=T) +
  stat_ellipse(aes(Circularity, Area, group=interaction(Strain,Morpho), fill=Strain, color=Strain), level=0.90, alpha=0.5, geom="polygon") +
  ylab(expression('Cell area')) + xlab(expression("Cell circularity")) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.5,by=0.5), limits=c(-0.1,1.5)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,5.0,by=1.0), limits=c(-2.5,5.0)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.background=element_blank(), strip.text.x=element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  facet_wrap(~Strain, ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
