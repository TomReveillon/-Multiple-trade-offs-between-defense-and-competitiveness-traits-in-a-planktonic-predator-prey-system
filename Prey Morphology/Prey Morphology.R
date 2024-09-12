setwd("~/LIMNO 2019-2023/Experiments/Prey Morphology")

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

#################################################################################
#################################################################################
##### ANALYSIS OF MORPHOLOGICAL TRAITS IN STRAINS WITH DIFFERENT MORPHOTYPES ####
#################################################################################
#################################################################################

# Load the datasets
Folder="Collection Strains"
Files=list.files(path=Folder, pattern="*.txt", full.names=T)
DataInter=ldply(Files, read.table, sep="\t", fill=T, header=T, dec=",")

# Verify dataset classes
lapply(DataInter, class)

# Specify the variables as numeric numbers
DataInter %<>% mutate_if(is.integer,as.numeric)
DataInter %<>% mutate_if(is.factor,as.numeric)

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

# Select columns
DataPCA=DataPCA[,c(1,3,5,8,9,12,13,27,39)]
DataPCA=sapply(DataPCA[,c(1:9)], as.numeric)
DataPCA=as.data.frame(DataPCA)

# Principal component analysis
PCA=dudi.pca(DataPCA, scannf=F, nf=2)

# Create a dataset of traits
DataFeat=data.frame(Trait=colnames(DataPCA),Category=colnames(DataPCA))

# Create categories of traits
DataFeat$Category %<>% gsub(c("^Area$|^Diameter$|^Height$|^Length$|^Perimeter$|^Width$"), "Size", .) %>% gsub(c("^Aspect Ratio$|^Circularity$|^Elongatedness$|^Lobe Count$"), "Shape", .)
DataFeat$Category=gsub("Circularity", "Roundness", DataFeat$Category)
DataFeat$Category=gsub("Elongatedness", "Elongation", DataFeat$Category)

# Extract trait coordinates
DataFeat$Coord1=PCA$co[,1]
DataFeat$Coord2=PCA$co[,2]
DataFeat$Contr1=PCA$c1[,1]
DataFeat$Contr2=PCA$c1[,2]

# Position of traits
tiff('PCA Traits.tiff', units="in", width=8, height=8, res=1000)
ggplot(DataFeat, aes(Coord1, Coord2, group=Trait)) +
  geom_hline(yintercept=0, color="grey50", linetype="dotted", size=1) + 
  geom_vline(xintercept=0, color="grey50", linetype="dotted", size=1) +
  geom_segment(aes(x=0, y=0, xend=Coord1, yend=Coord2), color="black", linetype="solid", size=0.7) +
  geom_point(aes(color=Category), pch=16, size=3) +
  geom_label_repel(aes(label=Trait), fill="white", color="black", segment.color=NA, size=6) +
  ylab(expression('Dimension 2 (18.67 %)')) + xlab(expression('Dimension 1 (69.43 %)')) +
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


###########################
#### Interstrain traits ###
###########################

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
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-10,10,by=5), limits=c(-10,10)) +
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
Panel=lapply(SplitIndiv, PlotFunc)
Yaxis=textGrob(expression('Dimension 2 (16.94 %)'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Dimension 1 (72.15 %)'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=3, nrow=2)
dev.off()

# Select relevant traits
FeatInter=data.frame(DataPCA[,c(1,3)],Strain)
FeatInter[,1]=FeatInter[,1]/10^3
FeatInter[,2]=FeatInter[,2]/10^1
colnames(FeatInter)=c("Area","Roundness","Strain")

# Sort and split the dataset
MeltInter=melt(FeatInter, id.vars=c("Strain"))
colnames(MeltInter)[2:3]=c("Trait","Value")

# Calculate mean trait values
MeltInter=setDT(MeltInter)[, Mean := mean(Value), by=list(Strain,Trait)] 
MeltInter=as.data.frame(MeltInter)

# Split the dataset
SplitInter=split(MeltInter, list(MeltInter$Trait))

# Violin plots of traits
PlotFunc=function(x) {
  ggplot(x, aes(Trait, Value, group=Strain)) +
    geom_violin(aes(fill=Strain, color=Strain)) +
    geom_point(aes(Trait, Mean, color=Strain), size=3, pch=16) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_blank()) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_blank()) +
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

tiff('PCA Traits Interstrain.tiff', units="in", width=15, height=8, res=1000)
Panel=lapply(SplitInter, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Particle area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.6,by=0.4), limits=c(0,1.6))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Particle roundness'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,4.8,by=1.2), limits=c(0,4.8))
Xaxis=textGrob(expression(italic('C. reinhardtii')~'strain'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=1)
dev.off()

# Split the dataset
SplitInter=split(MeltInter, list(MeltInter$Trait))

# Box plots of traits
PlotFunc=function(x) {
  ggplot(x, aes(Strain, Value, group=Strain)) +
    geom_boxplot(aes(fill=Strain, color=Strain), size=0.5, width=0.7, fatten=NA, outlier.shape=NA) +
    stat_boxplot(aes(color=Strain), geom="errorbar", width=0.2) +
    geom_point(aes(Strain, Mean, color=Strain), size=3, pch=16) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
    scale_x_discrete(labels=c("CR1"=expression(C[R1]),"CR2"=expression(C[R2]),"CR3"=expression(C[R3]),"CR4"=expression(C[R4]),"CR5"=expression(C[R5]),"CR6"=expression(C[R6]))) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.5)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(legend.position="none")
}

tiff('PCA Traits Interstrain.tiff', units="in", width=15, height=8, res=1000)
Panel=lapply(SplitInter, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Particle area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.6,by=0.4), limits=c(0,1.6))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Particle roundness'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,4.8,by=1.2), limits=c(0,4.8))
Xaxis=textGrob(expression(italic('C. reinhardtii')~'strain'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=1)
dev.off()


##############################################
### Morphotypes identification with traits ###
##############################################

# Add identities to the dataset
DataMorpho=DataPCA
DataMorpho$Strain=Strain

# Calculate confidence intervals
AreaCI=ci(subset(DataMorpho, Strain=="CR4")$Area, confidence=0.975)

# Calculate mean area of single cells
AreaM=mean(subset(DataMorpho, !Area < 107.49 & !Area > 108.59)$Area)

# Calculate cell numbers per image
DataMorpho$Cells=round(DataMorpho$Area/AreaM,0)
DataMorpho$Cells[DataMorpho$Cells == 0]=1

# Assign clump categories per image
Morpho=data.frame(DataMorpho %>% mutate(Clumps=case_when(Cells <= 2 ~ "S", Cells <= 6 ~ "SC", Cells <= 10 ~ "MC", Cells > 10 ~ "LC")))[,12]
Morpho=as.factor(Morpho)


###########################
#### Intrastrain traits ###
###########################

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
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(-10,10,by=5), limits=c(-10,10)) +
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
Panel=lapply(SplitIndiv, PlotFunc)
Yaxis=textGrob(expression('Dimension 2 (16.94 %)'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Xaxis=textGrob(expression('Dimension 1 (72.15 %)'), gp=gpar(fontface="bold", fontsize=18))
grid.arrange(grobs=Panel, left=Yaxis, bottom=Xaxis, ncol=3, nrow=2)
dev.off()

# Select relevant traits
FeatIntra=data.frame(DataPCA[,c(1,3)],Strain,Morpho)
FeatIntra[,1]=FeatIntra[,1]/10^3
FeatIntra[,2]=FeatIntra[,2]/10^1
colnames(FeatIntra)=c("Area","Roundness","Strain","Morpho")

# Sort and split the dataset
MeltIntra=melt(FeatIntra, id.vars=c("Strain","Morpho"))
colnames(MeltIntra)[3:4]=c("Trait","Value")

# Create codes binding strains and morphotypes
MeltIntra$Code=paste(MeltIntra$Strain, MeltIntra$Morpho, sep="")

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
MeltIntra$Morpho[MeltIntra$Code %in% c("CR3MC")]="MC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR4SC")]="SC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR5SC")]="SC"
MeltIntra$Morpho[MeltIntra$Code %in% c("CR6SC")]="SC"

# Specify the variables as factors
MeltIntra$Code=factor(MeltIntra$Code, levels=c("CR1S","CR1SC","CR1MC","CR1LC","CR2S","CR2SC","CR2MC","CR3S","CR3SC","CR3MC","CR4S","CR4SC","CR5S","CR5SC","CR6S","CR6SC"))

# Calculate mean trait values
MeltIntra=setDT(MeltIntra)[, Mean := mean(Value), by=list(Code,Trait)]
MeltIntra=as.data.frame(MeltIntra)

# Split the dataset
SplitIntra=split(MeltIntra, list(MeltIntra$Trait))

# Violin plots of traits
PlotFunc=function(x) {
  ggplot(x, aes(Trait, Value, group=Trait)) +
    geom_violin(aes(fill=Strain, color=Strain)) +
    geom_point(aes(Trait, Mean, color=Strain), size=3, pch=16) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_blank()) +
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_blank()) + 
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

tiff('PCA Traits Intrastrain.tiff', units="in", width=15, height=8, res=1000)
Panel=lapply(SplitIntra, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Particle area'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.6,by=0.4), limits=c(0,1.6))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Particle roundness'), labels=function(x) sprintf("%.1f", x), breaks=seq(0,4.8,by=1.2), limits=c(0,4.8))
Xaxis=textGrob(expression(italic('C. reinhardtii')~'morphotype'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=1)
dev.off()

# Split the dataset
SplitIntra=split(MeltIntra, list(MeltIntra$Code,MeltIntra$Trait))

# Box plots of traits
PlotFunc=function(x) {
  ggplot(x, aes(Code, Value, group=Code)) +
    geom_boxplot(aes(fill=Code, color=Code), size=0.5, width=0.7, fatten=NA, outlier.shape=NA) +
    stat_boxplot(aes(color=Code), geom="errorbar", width=0.2) +
    geom_point(aes(Code, Mean, color=Code), size=3, pch=16) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_discrete(labels=c("CR1S"=expression(C[R1S]),"CR1SC"=expression(C[R1SC]),"CR1MC"=expression(C[R1MC]),"CR1LC"=expression(C[R1LC]),"CR2S"=expression(C[R2S]),"CR2SC"=expression(C[R2SC]),"CR2MC"=expression(C[R2MC]),"CR3S"=expression(C[R3S]),"CR3SC"=expression(C[R3SC]),"CR3MC"=expression(C[R3MC]),"CR4S"=expression(C[R4S]),"CR4SC"=expression(C[R4SC]),"CR5S"=expression(C[R5S]),"CR5SC"=expression(C[R5SC]),"CR6S"=expression(C[R6S]),"CR6SC"=expression(C[R6SC]))) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1S"="mediumpurple3","CR1SC"="mediumpurple3","CR1MC"="mediumpurple3","CR1LC"="mediumpurple3","CR2S"="cornflowerblue","CR2SC"="cornflowerblue","CR2MC"="cornflowerblue","CR3S"="chartreuse3","CR3SC"="chartreuse3","CR3MC"="chartreuse3","CR4S"="gold2","CR4SC"="gold2","CR5S"="darkorange1","CR5SC"="darkorange1","CR6S"="tomato2","CR6SC"="tomato2"),0.5)) +
    scale_color_manual(values=c("CR1S"="mediumpurple3","CR1SC"="mediumpurple3","CR1MC"="mediumpurple3","CR1LC"="mediumpurple3","CR2S"="cornflowerblue","CR2SC"="cornflowerblue","CR2MC"="cornflowerblue","CR3S"="chartreuse3","CR3SC"="chartreuse3","CR3MC"="chartreuse3","CR4S"="gold2","CR4SC"="gold2","CR5S"="darkorange1","CR5SC"="darkorange1","CR6S"="tomato2","CR6SC"="tomato2")) +
    theme(legend.position="none")
}

tiff('PCA Traits Intrastrain.tiff', units="in", width=15, height=8, res=1000)
Panel1=lapply(SplitIntra[c(1:16)], PlotFunc)
Panel1=lapply(Panel1, function(x) {x=x + scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.6,by=0.4), limits=c(0,1.6))})
Yaxis=textGrob(expression('Particle area'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Grid1=grid.arrange(grobs=Panel1, left=Yaxis, ncol=4, nrow=4)
Panel2=lapply(SplitIntra[c(17:32)], PlotFunc)
Panel2=lapply(Panel2, function(x) {x=x + scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,4.8,by=1.2), limits=c(0,4.8))})
Yaxis=textGrob(expression('Particle roundness'), gp=gpar(fontface="bold", fontsize=18), rot=90)
Grid2=grid.arrange(grobs=Panel2, left=Yaxis, ncol=4, nrow=4)
Xaxis=textGrob(expression(italic('C. reinhardtii')~'morphotype'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(Grid1, Grid2, bottom=Xaxis, ncol=2, nrow=1)
dev.off()


###########################################
### Coefficients of variation in traits ###
###########################################

# Mutate dataset classes to numeric
DataCV=as.data.frame(sapply(DataPCA, as.numeric))
DataCV$Strain=Strain
DataCV$Morpho=Morpho

# Split the dataset
DataCV=DataCV %>% arrange(factor(Strain, levels=c("CR1","CR2","CR3","CR4","CR5","CR6")), factor(Morpho, levels=c("S","SC","MC","LC")))
SplitCV=split(DataCV, list(DataCV$Strain))
SplitCV=SplitCV[sapply(SplitCV, function(x) dim(x)[1]) > 0]

# Calculate interstrain coefficients
CVInter=abs(unlist(lapply(SplitCV, function(x) {apply(x[,c(1:9)], 2, function(x) sd(x)/mean(x)*100)})))
CVInter=data.frame(CVInter,Trait=rep(DataFeat$Trait,6))
rownames(CVInter)=c()

# Create a dataset per category
FeatCV=data.frame(CV=CVInter[,1],Trait=rep(DataFeat$Trait,6))
FeatCV$Strain=rep(unique(DataCV[c("Strain")])[,1], each=9)

# Statistical tests
shapiro.test(FeatCV$CV)
bartlett.test(FeatCV$CV~FeatCV$Strain)
kruskal.test(FeatCV$CV, FeatCV$Strain)
dunn.test(FeatCV$CV, FeatCV$Strain)

# Split the dataset
SplitFeatCV=list(FeatCV)

PlotFunc=function(x) {
  ggplot(x, aes(Strain, CV, group=Strain)) +
    geom_boxplot(aes(fill=Strain, color=Strain), size=0.5, width=0.7, outlier.shape=NA) + 
    geom_point(aes(color=Strain), size=1.5, pch=16, alpha=0.7, position=position_jitter(0.3)) +
    ylab("Coefficient of variation (%)") + xlab(expression(italic('C. reinhardtii')~'strain')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.x=element_text(face="plain", colour="black", size=18)) + 
    scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,100,by=25), limits=c(0,100)) +
    scale_x_discrete(labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R6]),expression(C[R7]))) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.3)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(legend.position="none")
}

tiff('CV Traits Interstrain.tiff', units="in", width=8, height=8, res=1000)
Panel=lapply(SplitFeatCV, PlotFunc)
grid.arrange(grobs=Panel, ncol=1, nrow=1)
dev.off()


############################################
### Interstrain allometric relationships ###
############################################

# Correlation matrix for traits
CorPCA=cor(DataPCA[,c(1:9)], method="spearman")

# Correlation plot for traits
corrplot(cor(DataPCA[,c(1:9)]), method="circle", order="FPC", type="lower", diag=F, 
col=colorRampPalette(c("darkred","white","darkgreen"))(100), addgrid.col="grey90",
cl.lim=c(-1.0, 1.0), cl.pos="b", cl.cex=1.0, tl.col="black", tl.cex=1.0, tl.srt=45)

# Correlation tests
CorFunc=function(x) {cor.test(x$Area, x$Roundness, method="pearson", exact=F)}
OutCor=lapply(split(FeatIntra, list(FeatIntra$Strain)), CorFunc)

# Create a labeling dataset
R2=c(as.character(expression(italic(r) *~'= -0.53')),as.character(expression(italic(r) *~'= -0.43')),as.character(expression(italic(r) *~'= -0.45')),as.character(expression(italic(r) *~'= 0.07')),as.character(expression(italic(r) *~'= 0.12')),as.character(expression(italic(r) *~'= -0.08')))
P=c(as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')),as.character(expression(italic(P) *~'< 0.001')))
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

tiff('Correlation Traits Interstrain.tiff', units="in", width=15, height=8, res=1000)
ggplot(FeatIntra, aes(Roundness, Area, group=Strain)) + 
  geom_point(aes(color=Strain), size=0.3, pch=16, alpha=0.5) +
  geom_smooth(method='lm', formula=y~x, aes(color=Strain), linetype="solid", size=0.8, se=F, fullrange=T) +
  ylab(expression('Particle area')) + xlab(expression("Particle roundness")) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.5,by=0.3), limits=c(-0.1,1.5)) +
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
