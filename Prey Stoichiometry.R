setwd("~/LIMNO 2019-2023/Experiments/Prey Stoichiometry")

rm(list=ls())

library(cowplot)
library(data.table)
library(deSolve)
library(directlabels)
library(dplyr)
library(dynlm)
library(foreach)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(lme4)
library(lmtest)
library(magrittr)
library(nlme)
library(plotly)
library(plyr)
library(propagate)
library(reshape2)
library(scales)

##############################################################################
##############################################################################
##### PREY STOICHIOMETRY FOR PREDATOR-PREY SYSTEM WITH DIFFERENT SPECIES #####
##############################################################################
##############################################################################

# Import the datasets
Data=read.table("Data_CN.txt", h=T, dec=",")
names(Data)
summary(Data)

# Specify the variables as numeric or factor
Data[,c(3:12)] %<>% mutate_if(is.character,as.numeric)

# Calculate mean controls
Data$WaterC=rep(ddply(Data, .(Strain), summarize, WaterC=round(mean(WaterC),6))[,2], each=3)
Data$WaterN=rep(ddply(Data, .(Strain), summarize, WaterN=round(mean(WaterN),6))[,2], each=3)
Data$MediaC=rep(ddply(Data, .(Strain), summarize, MediaC=round(mean(MediaC),6))[,2], each=3)
Data$MediaN=rep(ddply(Data, .(Strain), summarize, MediaN=round(mean(MediaN),6))[,2], each=3)

# Correct element masses
Data$MassC=Data$MassC-Data$MediaC
Data$MassN=Data$MassN-Data$MediaN

# Convert carbon and nitrogen masses
Data$MassC=Data$MassC*10^6
Data$MassN=Data$MassN*10^6

# Calculate cell densities
Data$Dens=round(Data$Cells*Data$Volu*Data$Site*Data$Dilu*Data$Cove,0)

# Calculate cell element masses
Data$CellC=round(Data$MassC/Data$Dens,6)
Data$CellN=round(Data$MassN/Data$Dens,6)

# Calculate C:N ratios
Data$CellCN=round(Data$CellC/Data$CellN,6)

# Melt the dataset
Data2=melt(Data[,c(1,6,15:17)], id.vars=c("Strain","Trial"))
colnames(Data2)[3:4]=c("Element","Mass")


############################################
### Plotting stoichiometric compositions ###
############################################

# Create a dataset
Data2=Data2[,c(1:4)]
Data2[,c(4)]=round(Data2[,c(4)],6)
Data2[,c(4)][Data2[,c(4)]<0]=0
Data2=Data2[order(Data2$Strain,Data2$Element),]

# Export the dataset
Data2[,c(4)]=replace(Data2[,c(4)],Data2[,c(4)]<0,0)
write.table(Data2, file="Data_SP.txt", sep="\t", row.names=F)

# Calculate mean trait values
Data2=setDT(Data2)[, Mean := mean(Mass), by=list(Strain,Element)]
Data2=as.data.frame(Data2)

# Split the dataset
SplitData2=split(Data2, list(Data2$Element))

PlotFunc=function(x) {
  ggplot(x, aes(Strain, Mass, group=Strain)) +
    geom_boxplot(aes(fill=Strain, color=Strain), size=0.5, width=0.7, fatten=NA, outlier.shape=NA) +
    stat_boxplot(aes(color=Strain), geom="errorbar", width=0.2) +
    geom_point(aes(Strain, Mean, color=Strain), size=3, pch=16, alpha=0.7) +
    ylab(expression('Carbon:nitrogen ratio ratio'~'('*ng~C~ng~N^-1*')')) + xlab(expression(italic('C. reinhardtii')~'strain')) +
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.x=element_blank()) +
    scale_y_continuous(labels=sprintf(seq(12,20,by=2), fmt="%.0f"), breaks=seq(12,20,by=2), limits=c(12,20)) +
    scale_x_discrete(labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R6]),expression(C[R7]))) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=alpha(c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2"),0.3)) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.y=element_blank()) +
    theme(strip.background=element_blank(), strip.text.x=element_text(face="plain", colour="black", size=18, angle=0, vjust=2, hjust=0.5)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1.2) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1.2) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Element, scales="free") +
    theme(legend.position="none")
}

tiff('Prey Stoichiometry.tiff', units="in", width=8, height=8, res=1000)
Panel=lapply(SplitData2[c(3)], PlotFunc)
grid.arrange(grobs=Panel, ncol=1, nrow=1)
dev.off()
