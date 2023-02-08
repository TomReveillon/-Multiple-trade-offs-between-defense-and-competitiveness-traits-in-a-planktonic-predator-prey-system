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

###################################################################################
###################################################################################
##### FUNCTIONAL RESPONSE FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
###################################################################################
###################################################################################

####################################################
### Estimation of functional response parameters ###
####################################################

# Import the dataset
DataFR=read.table("~/Activité Professionnelle/LIMNO 2019-2022/Experiments/Functional Response/Data_FRPODE.txt", h=T, dec=",")
summary(DataFR)
names(DataFR)

# Specify the variables as numeric or factor
DataFR[,c(2:5)] %<>% mutate_if(is.character,as.numeric)
DataFR$Strain=factor(DataFR$Strain, levels=unique(DataFR$Strain))

# Split the dataset
SplitDataFR=split(DataFR, list(DataFR$Strain))

# Extract combinations of names
Strain=unique(DataFR[,c("Strain")])

# Functional response model
FuncFR=function(x) {
  ModFR=nls(IngesP ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  ModFRL=nls(IngesPLSD ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  ModFRU=nls(IngesPUSD ~ (a * IDensP) / (1 + a * h * IDensP), start=c(a=1.0, h=0.1), data=x)
  Attack=data.frame(Attack=coef(ModFR)[1], AttackL=coef(ModFRL)[1], AttackU=coef(ModFRU)[1])
  Handling=data.frame(Handling=coef(ModFR)[2], HandlingL=coef(ModFRU)[2], HandlingU=coef(ModFRL)[2])
  Parameters=list(Attack=Attack, Handling=Handling)
}
OutFR=lapply(SplitDataFR, FuncFR)

# Calculate the attack rates
Attacks=bind_rows(lapply(OutFR, function (x) x[c("Attack")]))
Attacks=round(as.data.frame(do.call("rbind",Attacks)),4)
Attacks=cbind(Strain=unique(Strain),Attacks)
rownames(Attacks)=c()

# Calculate the handling times
Handlings=bind_rows(lapply(OutFR, function (x) x[c("Handling")]))
Handlings=round(as.data.frame(do.call("rbind",Handlings)),4)
Handlings=cbind(Strain=Strain,Handlings)
rownames(Handlings)=c()

# Calculate the minimum handling time
Handling=as.data.frame(setDT(subset(DataFR, IDensP == 15))[, .(Handling=round(1/IngesP,4), HandlingL=round(1/IngesPLSD,4), HandlingU=round(1/IngesPUSD,4)), by=list(Strain)])
HandlingMin=Handling[which.min(Handling$Handling),]
HandlingMin=as.data.frame(HandlingMin)

# Calculate the ingestion defenses
Ingestions=round(HandlingMin[,2]*data.frame(IngesP=setDT(subset(DataFR, IDensP == 15))$IngesP),4)
IngestionsL=round(HandlingMin[,3]*data.frame(IngesPLSD=setDT(subset(DataFR, IDensP == 15))$IngesPLSD),4)
IngestionsU=round(HandlingMin[,4]*data.frame(IngesPUSD=setDT(subset(DataFR, IDensP == 15))$IngesPUSD),4)

# Calculate the maximum attack rate
Attack=as.data.frame(setDT(subset(DataFR, IDensP <= 0.5))[, .(Attack=round(coef(lm(IngesP~IDensP))[2],4), AttackL=round(coef(lm(IngesPLSD~IDensP))[2],4), AttackU=round(coef(lm(IngesPUSD~IDensP))[2],4)), by=list(Strain)])
Attack[,2]=round(Attack[,2]/Ingestions[,1],4)
Attack[,3]=round(Attack[,3]/IngestionsL[,1],4)
Attack[,4]=round(Attack[,4]/IngestionsU[,1],4)
AttackMax=Attack[which.max(Attack$Attack),]
AttackMax=as.data.frame(AttackMax)

# Calculate the detection defenses
Detections=round(Attacks[,2]/(AttackMax[,2]*Ingestions),4)
DetectionsL=round(Attacks[,3]/(AttackMax[,3]*IngestionsU),4)
DetectionsU=round(Attacks[,4]/(AttackMax[,4]*IngestionsL),4)


#################################################
### Estimation of prey growth rate parameters ###
#################################################

# Import the dataset
DataAG=read.table("~/Activité Professionnelle/LIMNO 2019-2022/Experiments/Prey Growth/Data_AGP.txt", h=T, dec=",")
summary(DataAG)
names(DataAG)

# Specify the variables as numeric or factor
DataAG[,c(3:8)] %<>% mutate_if(is.character,as.numeric)
DataAG$Strain=factor(DataAG$Strain, levels=unique(DataAG$Strain))

# Select the time period
DataAG=as.data.frame(DataAG %>% group_by(Strain) %>% dplyr::slice(11:n()))

# Split the dataset
SplitDataAG=split(DataAG, list(DataAG$Strain))

# Extract combinations of names
Strain=unique(DataAG[,c("Strain")])

# Calculate the intrinsic growth rates
ModGR=function(x) {coef(summary(lm(log(DensP+1)~DayP, data=subset(x, DayP <= 3))))[2,1]}
ModGRSD=function(x) {coef(summary(lm(log(DensP+1)~DayP, data=subset(x, DayP <= 3))))[2,2]}

# Calculate the intrinsic growth rates
GR=round(c(do.call("rbind",lapply(SplitDataAG, ModGR))),4)
GRL=round(GR-c(do.call("rbind",lapply(SplitDataAG, ModGRSD))),4)
GRU=round(GR+c(do.call("rbind",lapply(SplitDataAG, ModGRSD))),4)
GrowAG=data.frame(Strain=Strain, GrowAG=GR, GrowAGL=GRL, GrowAGU=GRU)


#####################################################
### Estimation of predator growth rate parameters ###
#####################################################

# Import the dataset
DataRG=read.table("~/Activité Professionnelle/LIMNO 2019-2022/Experiments/Predator Growth/Data_RGP.txt", h=T, dec=",")
summary(DataRG)
names(DataRG)

# Specify the variables as numeric or factor
DataRG[,c(2:7)] %<>% mutate_if(is.character,as.numeric)
DataRG$Strain=factor(DataRG$Strain, levels=unique(DataRG$Strain))

# Select the time period
DataRG=as.data.frame(DataRG %>% group_by(Strain) %>% dplyr::slice(0:n()))

# Split the dataset
SplitDataRG=split(DataRG, list(DataRG$Strain))

# Extract combinations of names
Strain=unique(DataRG[,c("Strain")])

# Calculate the intrinsic growth rates
ModGR=function(x) {coef(summary(lm(DensP~DayP, data=subset(x, DayP <= 5))))[2,1]}
ModGRSD=function(x) {coef(summary(lm(DensP~DayP, data=subset(x, DayP <= 5))))[2,2]}

# Calculate the intrinsic growth rates
GR=round(c(do.call("rbind",lapply(SplitDataRG, ModGR))),4)
GRL=round(GR-c(do.call("rbind",lapply(SplitDataRG, ModGRSD))),4)
GRU=round(GR+c(do.call("rbind",lapply(SplitDataRG, ModGRSD))),4)
GrowRG=data.frame(Strain=Strain, GrowRG=GR, GrowRGL=GRL, GrowRGU=GRU)


########################################################
### Estimation of nitrogen affinities and morphology ###
########################################################

# Import the dataset
DataHS=read.table("~/Activité Professionnelle/LIMNO 2019-2022/Experiments/Prey Growth/Data_HSP.txt", h=T, dec=",")
summary(DataHS)
names(DataHS)

# Specify the variables as numeric or factor
DataHS[,c(2:10)] %<>% mutate_if(is.character,as.numeric)

# Morphological and functional parameters
Affin=DataHS$Affin; AffinL=DataHS$AffinL; AffinU=DataHS$AffinU
Satur=DataHS$Satur; SaturL=DataHS$SaturL; SaturU=DataHS$SaturU
GrowM=DataHS$Grow; GrowML=DataHS$GrowL; GrowMU=DataHS$GrowU
Area=c(3.23,1.78,1.71,1.08,1.21,1.36)
AreaL=c(0.93,0.60,0.60,0.73,0.67,0.91)
AreaU=c(5.53,2.96,2.82,1.42,1.74,1.81)

# Create the datasets
Data=data.frame(Strain=unique(DataFR$Strain),Attack=Attacks[,2],Handling=Handlings[,2],PreyG=GrowAG[,2],PredG=GrowRG[,2],Area=Area,Affin=Affin,Satur=Satur,PreyGM=GrowM) 
DataL=data.frame(Strain=unique(DataFR$Strain),AttackL=Attacks[,3],HandlingL=Handlings[,3],PreyGL=GrowAG[,3],PredGL=GrowRG[,3],AreaL=AreaL,AffinL=AffinL,SaturL=SaturL,PreyGML=GrowML)
DataU=data.frame(Strain=unique(DataFR$Strain),AttackU=Attacks[,4],HandlingU=Handlings[,4],PreyGU=GrowAG[,4],PredGU=GrowRG[,4],AreaU=AreaU,AffinU=AffinU,SaturU=SaturU,PreyGMU=GrowMU)
Data[,c(2:9)]=round(Data[,c(2:9)],4); DataL[,c(2:9)]=round(DataL[,c(2:9)],4); DataU[,c(2:9)]=round(DataU[,c(2:9)],4)

# Melt the dataset
MeltData=melt(Data, id.vars=c("Strain","PreyG","Affin","Satur","PreyGM")); colnames(MeltData)[6:7]=c("Trait","Value")
MeltDataL=melt(DataL, id.vars=c("Strain","PreyGL","AffinL","SaturL","PreyGML")); colnames(MeltDataL)[6:7]=c("TraitL","ValueL")
MeltDataU=melt(DataU, id.vars=c("Strain","PreyGU","AffinU","SaturU","PreyGMU")); colnames(MeltDataU)[6:7]=c("TraitU","ValueU")

# Combine the datasets
MeltData=cbind(MeltData[,c(1:7)],MeltDataL[,c(2:5,7)],MeltDataU[,c(2:5,7)])
MeltData=MeltData[,c("Strain","Trait","PreyG","PreyGL","PreyGU","PreyGM","PreyGML","PreyGMU","Affin","AffinL","AffinU","Satur","SaturL","SaturU","Value","ValueL","ValueU")]


##############################################################
### Correlation between defense and competitiveness traits ###
##############################################################

# Create a list of combined traits
Data2=transform(Data, Attack=rev(Attack), PredG=rev(PredG), Satur=Satur)
Data3=transform(Data, Attack=rev(Attack), PredG=rev(PredG), Satur=rev(Satur))
Data4=transform(Data, Handling=rev(Handling), Area=rev(Area), Satur=rev(Satur))
CorListI1=list(Data[,c(2,9,7,8)],Data[,c(3,9,7,8)],Data[,c(5,9,7,8)],Data[,c(6,9,7,8)])
CorListI2=list(Data3[,c(2,9,7,8)],Data4[,c(3,9,7,8)],Data3[,c(5,9,7,8)],Data4[,c(6,9,7,8)])
CorListS1=list(Data3[,c(2,9,7,8)],Data[,c(3,9,7,8)],Data3[,c(5,9,7,8)],Data[,c(6,9,7,8)])
CorListS2=list(Data2[,c(2,9,7,8)],Data[,c(3,9,7,8)],Data2[,c(5,9,7,8)],Data[,c(6,9,7,8)])

# Calculate regression intercepts
PreyGI=do.call("rbind",lapply(CorListI1, function(x) {coef(lm(x[,1]~x[,2]))}))[,1]
AffinI=do.call("rbind",lapply(CorListI1, function(x) {coef(lm(x[,1]~x[,3]))}))[,1]
SaturI=do.call("rbind",lapply(CorListI2, function(x) {coef(lm(x[,1]~x[,4]))}))[,1]

# Calculate regression slopes
PreyGS=do.call("rbind",lapply(CorListS1, function(x) {coef(lm(x[,1]~x[,2]))}))[,2]
AffinS=do.call("rbind",lapply(CorListS1, function(x) {coef(lm(x[,1]~x[,3]))}))[,2]
SaturS=do.call("rbind",lapply(CorListS2, function(x) {coef(lm(x[,1]~x[,4]))}))[,2]

# Order regression coefficients
PreyGI=c(-PreyGI[1],PreyGI[2],-PreyGI[3],PreyGI[4])
PreyGS=c(PreyGS[1],PreyGS[2],PreyGS[3],PreyGS[4])
AffinI=c(-AffinI[1],AffinI[2],-AffinI[3],AffinI[4])
AffinS=c(AffinS[1],AffinS[2],AffinS[3],AffinS[4])
SaturI=c(-SaturI[1],SaturI[2],-SaturI[3],SaturI[4])
SaturS=c(-SaturS[1],-SaturS[2],-SaturS[3],-SaturS[4])

# Include regression coefficients
MeltData$PreyGI=rep(round(PreyGI,4), each=6)
MeltData$PreyGS=rep(round(PreyGS,4), each=6)
MeltData$AffinI=rep(round(AffinI,4), each=6)
MeltData$AffinS=rep(round(AffinS,4), each=6)
MeltData$SaturI=rep(round(SaturI,4), each=6)
MeltData$SaturS=rep(round(SaturS,4), each=6)


######################################################
### Trade-off plots of defense and competitiveness ###
######################################################

# Split the dataset
SplitData=split(MeltData, list(MeltData$Trait))

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) +
    geom_abline(aes(intercept=PreyGI, slope=PreyGS), color="grey50", linetype="solid", size=0.6) +
    geom_errorbar(aes(PreyGM, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(PreyGM, Value, xmin=PreyGML, xmax=PreyGMU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(PreyGM, Value, color=Strain), fill="white", size=3, stroke=1.2, pch=21) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=sprintf(seq(0.09,0.36,by=0.09), fmt="%.2f"), breaks=seq(0.09,0.36,by=0.09), limits=c(0.08,0.373)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=2, nrow=2) +
    theme(legend.position="none")
}

tiff('Trait Spaces 1.tiff', units="in", width=12, height=12, res=1000)
Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_reverse(expression('Attack rate'~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(1.8,0,by=-0.6), fmt="%.1f"), breaks=c(seq(0.18,0,by=-0.06)), limits=c(0.18,0))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0,9.0))
Panel[[3]]=Panel[[3]] + scale_y_reverse(expression('Predator growth rate'~'('*day^-1*')'), labels=sprintf(seq(6.0,0,by=-2.0), fmt="%.1f"), breaks=c(seq(6.0,0,by=-2.0)), limits=c(6.0,-0.3))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Cell area'~'('*10^3~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0,6.0))
Xaxis=textGrob(expression('Prey maximum growth rate'~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=2)
dev.off()


PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) +
    geom_abline(aes(intercept=AffinI, slope=AffinS), color="grey50", linetype="solid", size=0.6) +
    geom_errorbar(aes(Affin, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(Affin, Value, xmin=AffinL, xmax=AffinU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(Affin, Value, color=Strain), fill="white", size=3, stroke=1.2, pch=21) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.x=element_blank()) +
    scale_x_continuous(labels=sprintf(seq(0,0.6,by=0.2), fmt="%.1f"), breaks=seq(0,0.6,by=0.2), limits=c(0,0.630)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=2, nrow=2) +
    theme(legend.position="none")
}

tiff('Trait Spaces 2.tiff', units="in", width=12, height=12, res=1000)
Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_reverse(expression('Attack rate'~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(1.8,0,by=-0.6), fmt="%.1f"), breaks=c(seq(0.18,0,by=-0.06)), limits=c(0.18,0))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0,9.0))
Panel[[3]]=Panel[[3]] + scale_y_reverse(expression('Predator growth rate'~'('*day^-1*')'), labels=sprintf(seq(6.0,0,by=-2.0), fmt="%.1f"), breaks=c(seq(6.0,0,by=-2.0)), limits=c(6.0,-0.3))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Cell area'~'('*10^3~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0,6.0))
Xaxis=textGrob(expression('Prey nitrate affinity'~'('*µM~NO[3]^{'-'}~L^-1~day^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=2)
dev.off()


PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) +
    geom_abline(aes(intercept=SaturI, slope=SaturS), color="grey50", linetype="solid", size=0.6) +
    geom_errorbar(aes(Satur, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_errorbar(aes(Satur, Value, xmin=SaturL, xmax=SaturU, color=Strain), linetype="solid", alpha=0.7, size=1.2, width=0) +
    geom_point(aes(Satur, Value, color=Strain), fill="white", size=3, stroke=1.2, pch=21) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
    theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
    theme(axis.title.x=element_blank()) +
    scale_x_reverse(labels=sprintf(seq(2.0,0.5,by=-0.5), fmt="%.1f"), breaks=seq(2.0,0.5,by=-0.5), limits=c(2.075,0.425)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=2, nrow=2) +
    theme(legend.position="none")
}

tiff('Trait Spaces 3.tiff', units="in", width=12, height=12, res=1000)
Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_reverse(expression('Attack rate'~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(1.8,0,by=-0.6), fmt="%.1f"), breaks=c(seq(0.18,0,by=-0.06)), limits=c(0.18,0))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0,9.0))
Panel[[3]]=Panel[[3]] + scale_y_reverse(expression('Predator growth rate'~'('*day^-1*')'), labels=sprintf(seq(6.0,0,by=-2.0), fmt="%.1f"), breaks=c(seq(6.0,0,by=-2.0)), limits=c(6.0,-0.3))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Cell area'~'('*10^3~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0,6.0))
Xaxis=textGrob(expression('Prey nitrate half-saturation'~'('*µM~NO[3]^{'-'}~L^-1*')'), gp=gpar(fontface="bold", fontsize=18), rot=0)
grid.arrange(grobs=Panel, bottom=Xaxis, ncol=2, nrow=2)
dev.off()
