setwd("~/LIMNO 2019-2022/Experiments")

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

###########################################################################
###########################################################################
##### TRAIT SPACE FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
###########################################################################
###########################################################################

MeltData$Strain=gsub("CR6", "CR7", MeltData$Strain)
MeltData$Strain=gsub("CR5", "CR6", MeltData$Strain)
MeltData=droplevels(subset(MeltData, !Trait=="PredG"))
Labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R6]),expression(C[R7]))
SplitData=split(MeltData, list(MeltData$Trait))

# Simulate regression solutions
DataA=data.frame(
  XPreyG=seq(0.09,0.36,length.out=1000),
  XAffin=seq(0.00,0.60,length.out=1000),
  XSatur=seq(0.40,2.80,length.out=1000),
  YPreyG=SplitData[[1]]$InterPreyG + SplitData[[1]]$SlopePreyG*seq(0.09,0.36,length.out=1000),
  YAffin=SplitData[[1]]$InterAffin + SplitData[[1]]$SlopeAffin*seq(0.00,0.60,length.out=1000),
  YSatur=SplitData[[1]]$InterSatur + SplitData[[1]]$SlopeSatur*seq(0.40,2.80,length.out=1000))

DataH=data.frame(
  XPreyG=seq(0.09,0.36,length.out=1000),
  XAffin=seq(0.00,0.60,length.out=1000),
  XSatur=seq(0.40,2.80,length.out=1000),
  YPreyG=SplitData[[2]]$InterPreyG + SplitData[[2]]$SlopePreyG*seq(0.09,0.36,length.out=1000),
  YAffin=SplitData[[2]]$InterAffin + SplitData[[2]]$SlopeAffin*seq(0.00,0.60,length.out=1000),
  YSatur=SplitData[[2]]$InterSatur + SplitData[[2]]$SlopeSatur*seq(0.40,2.80,length.out=1000))

DataS=data.frame(
  XPreyG=seq(0.09,0.36,length.out=1000),
  XAffin=seq(0.00,0.60,length.out=1000),
  XSatur=seq(0.40,2.80,length.out=1000),
  YPreyG=SplitData[[3]]$InterPreyG + SplitData[[3]]$SlopePreyG*seq(0.09,0.36,length.out=1000),
  YAffin=SplitData[[3]]$InterAffin + SplitData[[3]]$SlopeAffin*seq(0.00,0.60,length.out=1000),
  YSatur=SplitData[[3]]$InterSatur + SplitData[[3]]$SlopeSatur*seq(0.40,2.80,length.out=1000))

DataCN=data.frame(
  XPreyG=seq(0.09,0.36,length.out=1000),
  XAffin=seq(0.00,0.60,length.out=1000),
  XSatur=seq(0.40,2.80,length.out=1000),
  YPreyG=SplitData[[4]]$InterPreyG + SplitData[[4]]$SlopePreyG*seq(0.09,0.36,length.out=1000),
  YAffin=SplitData[[4]]$InterAffin + SplitData[[4]]$SlopeAffin*seq(0.00,0.60,length.out=1000),
  YSatur=SplitData[[4]]$InterSatur + SplitData[[4]]$SlopeSatur*seq(0.40,2.80,length.out=1000))

# Cut the regression solutions
PreyGLim=DataA[DataA$YPreyG >= 0.00 & DataA$YPreyG <= 0.18 & DataA$XPreyG >= 0.09 & DataA$XPreyG <= 0.36,]
PreyGYLimL=head(PreyGLim,1)[,4]
PreyGYLimU=tail(PreyGLim,1)[,4]
PreyGXLimL=head(PreyGLim,1)[,1]
PreyGXLimU=tail(PreyGLim,1)[,1]
AffinLim=DataA[DataA$YAffin >= 0.00 & DataA$YAffin <= 0.18 & DataA$XAffin >= 0.00 & DataA$XAffin <= 0.60,]
AffinYLimL=head(AffinLim,1)[,5]
AffinYLimU=tail(AffinLim,1)[,5]
AffinXLimL=head(AffinLim,1)[,2]
AffinXLimU=tail(AffinLim,1)[,2]
SaturLim=DataA[DataA$YSatur >= 0.00 & DataA$YSatur <= 0.18 & DataA$XSatur >= 0.40 & DataA$XSatur <= 2.80,]
SaturYLimL=head(SaturLim,1)[,6]
SaturYLimU=tail(SaturLim,1)[,6]
SaturXLimL=head(SaturLim,1)[,3]
SaturXLimU=tail(SaturLim,1)[,3]
DataA2=data.frame(Trait="Attack",PreyGYLimL,PreyGYLimU,PreyGXLimL,PreyGXLimU,AffinYLimL,AffinYLimU,AffinXLimL,AffinXLimU,SaturYLimL,SaturYLimU,SaturXLimL,SaturXLimU)

PreyGLim=DataH[DataH$YPreyG >= 0.00 & DataH$YPreyG <= 9.00 & DataH$XPreyG >= 0.09 & DataH$XPreyG <= 0.36,]
PreyGYLimL=head(PreyGLim,1)[,4]
PreyGYLimU=tail(PreyGLim,1)[,4]
PreyGXLimL=head(PreyGLim,1)[,1]
PreyGXLimU=tail(PreyGLim,1)[,1]
AffinLim=DataH[DataH$YAffin >= 0.00 & DataH$YAffin <= 9.00 & DataH$XAffin >= 0.00 & DataH$XAffin <= 0.60,]
AffinYLimL=head(AffinLim,1)[,5]
AffinYLimU=tail(AffinLim,1)[,5]
AffinXLimL=head(AffinLim,1)[,2]
AffinXLimU=tail(AffinLim,1)[,2]
SaturLim=DataH[DataH$YSatur >= 0.00 & DataH$YSatur <= 9.00 & DataH$XSatur >= 0.40 & DataH$XSatur <= 2.80,]
SaturYLimL=head(SaturLim,1)[,6]
SaturYLimU=tail(SaturLim,1)[,6]
SaturXLimL=head(SaturLim,1)[,3]
SaturXLimU=tail(SaturLim,1)[,3]
DataH2=data.frame(Trait="Handling",PreyGYLimL,PreyGYLimU,PreyGXLimL,PreyGXLimU,AffinYLimL,AffinYLimU,AffinXLimL,AffinXLimU,SaturYLimL,SaturYLimU,SaturXLimL,SaturXLimU)

PreyGLim=DataS[DataS$YPreyG >= 0.00 & DataS$YPreyG <= 6.00 & DataS$XPreyG >= 0.09 & DataS$XPreyG <= 0.36,]
PreyGYLimL=head(PreyGLim,1)[,4]
PreyGYLimU=tail(PreyGLim,1)[,4]
PreyGXLimL=head(PreyGLim,1)[,1]
PreyGXLimU=tail(PreyGLim,1)[,1]
AffinLim=DataS[DataS$YAffin >= 0.00 & DataS$YAffin <= 6.00 & DataS$XAffin >= 0.00 & DataS$XAffin <= 0.60,]
AffinYLimL=head(AffinLim,1)[,5]
AffinYLimU=tail(AffinLim,1)[,5]
AffinXLimL=head(AffinLim,1)[,2]
AffinXLimU=tail(AffinLim,1)[,2]
SaturLim=DataS[DataS$YSatur >= 0.00 & DataS$YSatur <= 6.00 & DataS$XSatur >= 0.40 & DataS$XSatur <= 2.80,]
SaturYLimL=head(SaturLim,1)[,6]
SaturYLimU=tail(SaturLim,1)[,6]
SaturXLimL=head(SaturLim,1)[,3]
SaturXLimU=tail(SaturLim,1)[,3]
DataS2=data.frame(Trait="Area",PreyGYLimL,PreyGYLimU,PreyGXLimL,PreyGXLimU,AffinYLimL,AffinYLimU,AffinXLimL,AffinXLimU,SaturYLimL,SaturYLimU,SaturXLimL,SaturXLimU)

PreyGLim=DataCN[DataCN$YPreyG >= 12 & DataCN$YPreyG <= 18 & DataCN$XPreyG >= 0.09 & DataCN$XPreyG <= 0.36,]
PreyGYLimL=head(PreyGLim,1)[,4]
PreyGYLimU=tail(PreyGLim,1)[,4]
PreyGXLimL=head(PreyGLim,1)[,1]
PreyGXLimU=tail(PreyGLim,1)[,1]
AffinLim=DataCN[DataCN$YAffin >= 12 & DataCN$YAffin <= 18 & DataCN$XAffin >= 0.00 & DataCN$XAffin <= 0.60,]
AffinYLimL=head(AffinLim,1)[,5]
AffinYLimU=tail(AffinLim,1)[,5]
AffinXLimL=head(AffinLim,1)[,2]
AffinXLimU=tail(AffinLim,1)[,2]
SaturLim=DataCN[DataCN$YSatur >= 12 & DataCN$YSatur <= 18 & DataCN$XSatur >= 0.40 & DataCN$XSatur <= 2.80,]
SaturYLimL=tail(SaturLim,1)[,6]
SaturYLimU=head(SaturLim,1)[,6]
SaturXLimL=tail(SaturLim,1)[,3]
SaturXLimU=head(SaturLim,1)[,3]
DataCN2=data.frame(Trait="Stoichio",PreyGYLimL,PreyGYLimU,PreyGXLimL,PreyGXLimU,AffinYLimL,AffinYLimU,AffinXLimL,AffinXLimU,SaturYLimL,SaturYLimU,SaturXLimL,SaturXLimU)

# Include limits in the dataset
DataLim=rbind(DataA2,DataH2,DataS2,DataCN2)
DataLim=DataLim[rep(seq_len(nrow(DataLim)),each=6),]
MeltData2=cbind(MeltData,DataLim[,-1])
SplitData2=split(MeltData2, list(MeltData2$Trait))

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=PreyGXLimL, xend=PreyGXLimU, y=PreyGYLimL, yend=PreyGYLimU, linetype=SigniPreyG, size=SigniPreyG), color="grey50") + 
    geom_errorbar(aes(PreyGM, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(PreyGM, Value, xmin=PreyGML, xmax=PreyGMU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(PreyGM, Value, color=Strain), fill="white", size=4, pch=16) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
    theme(axis.text.y.right=element_text(face="plain", colour="black", size=20)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
    theme(axis.title.y=element_blank()) + 
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=sprintf(seq(0.09,0.36,by=0.09), fmt="%.2f"), breaks=seq(0.09,0.36,by=0.09), limits=c(0.0684,0.3816)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_linetype_manual(values=c("Yes"="solid","No"=NA)) + scale_size_manual(values=c("Yes"=1,"No"=1.2)) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=1, nrow=4) +
    theme(legend.position="none")
}

Panel=lapply(SplitData2, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a)~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h)~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle size'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=0.18, x=0.36, color="black", size=6, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label=NA, y=9.0, x=0.36, color="black", size=6, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=0.36, color="black", size=6, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=0.36, color="black", size=6, fontface="italic")
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R1]), x=0, y=0.1800, color="mediumpurple3", size=5) + annotate("text", label=expression(C[R2]), x=0, y=0.1728, color="cornflowerblue", size=5)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R3]), x=0, y=0.1656, color="chartreuse3", size=5) + annotate("text", label=expression(C[R4]), x=0, y=0.1584, color="gold2", size=5)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R6]), x=0, y=0.036, color="darkorange1", size=5) + annotate("text", label=expression(C[R7]), x=0, y=0.1440, color="tomato2", size=5)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=20), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) 
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) 
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0,0,0.3),"cm"))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=20))
Panel[[4]]=Panel[[4]] + xlab(expression('Maximum growth rate'~italic(r[C~max])~'('*day^-1*')'))
Plot1=grid.arrange(grobs=Panel, ncol=1, nrow=4)

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=AffinXLimL, xend=AffinXLimU, y=AffinYLimL, yend=AffinYLimU, linetype=SigniAffin, size=SigniAffin), color="grey50") + 
    geom_errorbar(aes(Affin, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(Affin, Value, xmin=AffinL, xmax=AffinU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(Affin, Value, color=Strain), fill="white", size=4, pch=16) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=sprintf(seq(0,0.6,by=0.2), fmt="%.1f"), breaks=seq(0,0.6,by=0.2), limits=c(0,0.648)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_linetype_manual(values=c("Yes"="solid","No"=NA)) + scale_size_manual(values=c("Yes"=1,"No"=1.2)) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=1, nrow=4) +
    theme(legend.position="none")
}

Panel=lapply(SplitData2, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a)~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h)~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle size'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=0.18, x=0.6, color="black", size=6, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label=NA, y=9.0, x=0.6, color="black", size=6, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=0.6, color="black", size=6, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=0.6, color="black", size=6, fontface="italic")
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=20), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) 
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) 
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0,0,0.3),"cm")) 
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(vjust=1))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=20))
Panel[[4]]=Panel[[4]] + xlab(expression('Affinity constant'~italic(f[C])~'('*µM~NO[3]^{'-'}~L^-1~day^-1*')'))
Plot2=grid.arrange(grobs=Panel, ncol=1, nrow=4)

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=SaturXLimL, xend=SaturXLimU, y=SaturYLimL, yend=SaturYLimU, linetype=SigniSatur, size=SigniSatur), color="grey50") + 
    geom_errorbar(aes(Satur, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(Satur, Value, xmin=SaturL, xmax=SaturU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(Satur, Value, color=Strain), fill="white", size=4, pch=16) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=20)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=20)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=sprintf(seq(0.4,2.8,by=0.8), fmt="%.1f"), breaks=seq(0.4,2.8,by=0.8), limits=c(0.4,3.024)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_linetype_manual(values=c("Yes"="solid","No"=NA)) + scale_size_manual(values=c("Yes"=1,"No"=1.2)) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=1, nrow=4) +
    theme(legend.position="none")
}

Panel=lapply(SplitData2, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a)~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h)~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle size'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label="NS", y=0.18, x=2.8, color="black", size=6, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label="NS", y=9.0, x=2.8, color="black", size=6, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=2.8, color="black", size=6, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=2.8, color="black", size=6, fontface="italic")
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=20), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(2.0,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=20), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=20), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=20), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0,0,0.3),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0,0,0.3),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=20), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(vjust=1))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=20))
Panel[[4]]=Panel[[4]] + xlab(expression('Half-saturation constant'~italic(K[C])~'('*µM~NO[3]^{'-'}~L^-1*')'))
Plot3=grid.arrange(grobs=Panel, ncol=1, nrow=4)

tiff('[Labels] Trait Spaces.tiff', units="in", width=20, height=25, res=1000)
grid.arrange(Plot1, Plot2, Plot3, ncol=3, nrow=1)
dev.off()
