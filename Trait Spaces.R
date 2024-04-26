setwd("~/LIMNO 2019-2023/Experiments")

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

####################################################
### Estimation of functional response parameters ###
####################################################

# Import the dataset
DataFR=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Predator Ingestion/Data_FRPODE.txt", h=T, dec=",")
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
Attacks=cbind(Strain=Strain,Attacks)
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
DataAG=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Prey Growth/Data_AGP.txt", h=T, dec=",")
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
ModGR=function(x) {coef(summary(lm(log(DensP+1)~DayP, data=subset(x, DayP <= 5))))[2,1]}
ModGRSD=function(x) {coef(summary(lm(log(DensP+1)~DayP, data=subset(x, DayP <= 5))))[2,2]}

# Calculate the intrinsic growth rates
GR=round(c(do.call("rbind",lapply(SplitDataAG, ModGR))),4)
GRL=round(GR-c(do.call("rbind",lapply(SplitDataAG, ModGRSD))),4)
GRU=round(GR+c(do.call("rbind",lapply(SplitDataAG, ModGRSD))),4)
GrowAG=data.frame(Strain=Strain, GrowAG=GR, GrowAGL=GRL, GrowAGU=GRU)


#####################################################
### Estimation of predator growth rate parameters ###
#####################################################

# Import the dataset
DataRG=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Predator Growth/Data_RGP.txt", h=T, dec=",")
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


####################################################
### Estimation of prey nitrogen absorption rates ###
####################################################

# Import the dataset
DataHS=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Prey Growth/Data_HSP.txt", h=T, dec=",")
summary(DataHS)
names(DataHS)

# Specify the variables as numeric or factor
DataHS[,c(2:10)] %<>% mutate_if(is.character,as.numeric)

# Calculate maximum growth rates
GrowM=DataHS$Grow
GrowML=DataHS$GrowL
GrowMU=DataHS$GrowU

# Calculate affinity constants
Affin=DataHS$Affin
AffinL=DataHS$AffinL
AffinU=DataHS$AffinU

# Calculate half-saturation constants
Satur=DataHS$Satur
SaturL=DataHS$SaturL
SaturU=DataHS$SaturU


#######################################################
### Estimation of prey morphology and stoichiometry ###
#######################################################

# Import the dataset
DataST=read.table("~/Activité Professionnelle/LIMNO 2019-2023/Experiments/Prey Stoichiometry/Data_SP.txt", h=T, dec=",")
summary(DataST)
names(DataST)

# Specify the variables as numeric or factor
DataST[,c(2,4)] %<>% mutate_if(is.character,as.numeric)

# Calculate mean element ratios
DataST=subset(DataST, Element=="CellCN")
DataST=setDT(DataST)[, .(Mass=mean(Mass),MassL=mean(Mass)-sd(Mass),MassU=mean(Mass)+sd(Mass)), by=list(Strain)]
DataST=as.data.frame(DataST)

# Calculate particle sizes
Area=c(3.23,1.78,1.71,1.08,1.21,1.36)
AreaL=c(0.93,0.60,0.60,0.73,0.67,0.91)
AreaU=c(5.53,2.96,2.82,1.42,1.74,1.81)

# Calculate stoichiometric ratios
Stoichio=c(round(DataST[,2],2))
StoichioL=c(round(DataST[,3],2))
StoichioU=c(round(DataST[,4],2))


############################
### Combine the datasets ###
############################

# Create the datasets
Data=data.frame(Strain=unique(DataFR$Strain),Attack=Attacks[,2],Handling=Handlings[,2],PreyG=GrowAG[,2],PredG=GrowRG[,2],Area=Area,Stoichio=Stoichio,PreyGM=GrowM,Affin=Affin,Satur=Satur) 
DataL=data.frame(Strain=unique(DataFR$Strain),AttackL=Attacks[,3],HandlingL=Handlings[,3],PreyGL=GrowAG[,3],PredGL=GrowRG[,3],AreaL=AreaL,StoichioL=StoichioL,PreyGML=GrowML,AffinL=AffinL,SaturL=SaturL)
DataU=data.frame(Strain=unique(DataFR$Strain),AttackU=Attacks[,4],HandlingU=Handlings[,4],PreyGU=GrowAG[,4],PredGU=GrowRG[,4],AreaU=AreaU,StoichioU=StoichioU,PreyGMU=GrowMU,AffinU=AffinU,SaturU=SaturU)
Data[,c(2:10)]=round(Data[,c(2:10)],4); DataL[,c(2:10)]=round(DataL[,c(2:10)],4); DataU[,c(2:10)]=round(DataU[,c(2:10)],4)

# Melt the dataset
MeltData=melt(Data, id.vars=c("Strain","PreyG","PreyGM","Affin","Satur")); colnames(MeltData)[6:7]=c("Trait","Value")
MeltDataL=melt(DataL, id.vars=c("Strain","PreyGL","PreyGML","AffinL","SaturL")); colnames(MeltDataL)[6:7]=c("TraitL","ValueL")
MeltDataU=melt(DataU, id.vars=c("Strain","PreyGU","PreyGMU","AffinU","SaturU")); colnames(MeltDataU)[6:7]=c("TraitU","ValueU")

# Combine the datasets
MeltData=cbind(MeltData[,c(1:7)],MeltDataL[,c(2:5,7)],MeltDataU[,c(2:5,7)])
MeltData=MeltData[,c("Strain","Trait","PreyG","PreyGL","PreyGU","PreyGM","PreyGML","PreyGMU","Affin","AffinL","AffinU","Satur","SaturL","SaturU","Value","ValueL","ValueU")]
MeltData=subset(MeltData, !Trait=="PredG")
MeltData=droplevels(MeltData)


##############################################################
### Correlation between defense and competitiveness traits ###
##############################################################

# Create a list of traits
ListData=list(Data[,c(2,8,9,10)],Data[,c(3,8,9,10)],Data[,c(6,8,9,10)],Data[,c(7,8,9,10)])

# Calculate regression intercepts
MeltData$InterPreyG=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,2]))}))[,1],4), each=6)
MeltData$InterAffin=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,3]))}))[,1],4), each=6)
MeltData$InterSatur=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,4]))}))[,1],4), each=6)

# Calculate regression slopes
MeltData$SlopePreyG=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,2]))}))[,2],4), each=6)
MeltData$SlopeAffin=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,3]))}))[,2],4), each=6)
MeltData$SlopeSatur=rep(round(do.call("rbind",lapply(ListData, function(x) {coef(lm(x[,1]~x[,4]))}))[,2],4), each=6)

# Calculate correlation coefficients
MeltData$CorPreyG=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,2])[[4]]})),4), each=6)
MeltData$CorAffin=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,3])[[4]]})),4), each=6)
MeltData$CorSatur=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,4])[[4]]})),4), each=6)
MeltData$SigPreyG=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,2])[[3]]})),4), each=6)
MeltData$SigAffin=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,3])[[3]]})),4), each=6)
MeltData$SigSatur=rep(round(do.call("rbind",lapply(ListData, function(x) {cor.test(x[,1],x[,4])[[3]]})),4), each=6)

# Identify non-Sigficant correlations
MeltData$SigPreyG=ifelse(MeltData$SigPreyG > 0.05, "No", "Yes")
MeltData$SigAffin=ifelse(MeltData$SigAffin > 0.05, "No", "Yes")
MeltData$SigSatur=ifelse(MeltData$SigSatur > 0.05, "No", "Yes")


################################################
### Trade-off of defense and competitiveness ###
################################################

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
MeltData=cbind(MeltData,DataLim[,-1])
SplitData=split(MeltData, list(MeltData$Trait))


######################################################
### Trade-off plots of defense and competitiveness ###
######################################################

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=PreyGXLimL, xend=PreyGXLimU, y=PreyGYLimL, yend=PreyGYLimU, linetype=SigPreyG, size=SigPreyG), color="grey50", size=1.5) + 
    geom_errorbar(aes(PreyGM, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_errorbar(aes(PreyGM, Value, xmin=PreyGML, xmax=PreyGMU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_point(aes(PreyGM, Value, color=Strain), fill="white", size=7, pch=20) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=25)) +  
    theme(axis.text.y.right=element_text(face="plain", colour="black", size=25)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=25)) + 
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

Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a[B])~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h[B])~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle area'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=0.18, x=0.36, color="black", size=8, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label=NA, y=9.0, x=0.36, color="black", size=8, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=0.36, color="black", size=8, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=0.36, color="black", size=8, fontface="italic")
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R1]), x=0.09, y=0.1800, color="mediumpurple3", size=7) + annotate("text", label=expression(C[R2]), x=0.09, y=0.1656, color="cornflowerblue", size=7)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R3]), x=0.09, y=0.1512, color="chartreuse3", size=7) + annotate("text", label=expression(C[R4]), x=0.09, y=0.1368, color="gold2", size=7)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R6]), x=0.09, y=0.1224, color="darkorange1", size=7) + annotate("text", label=expression(C[R7]), x=0.09, y=0.1080, color="tomato2", size=7)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=25), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) 
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) 
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0.2,0,0.35),"cm"))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=25))
Panel[[4]]=Panel[[4]] + xlab(expression('Maximum growth rate'~italic(r[C~max])~'('*day^-1*')'))
Plot1=grid.arrange(grobs=Panel, ncol=1, nrow=4)

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=AffinXLimL, xend=AffinXLimU, y=AffinYLimL, yend=AffinYLimU, linetype=SigAffin, size=SigAffin), color="grey50", size=1.5) + 
    geom_errorbar(aes(Affin, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_errorbar(aes(Affin, Value, xmin=AffinL, xmax=AffinU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_point(aes(Affin, Value, color=Strain), fill="white", size=7, pch=20) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=25)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=25)) + 
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

Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a)~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h)~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle area'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=0.18, x=0.6, color="black", size=8, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label=NA, y=9.0, x=0.6, color="black", size=8, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=0.6, color="black", size=8, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=0.6, color="black", size=8, fontface="italic")
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=25), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) 
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) 
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0.2,0,0.35),"cm")) 
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(vjust=1))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=25))
Panel[[4]]=Panel[[4]] + xlab(expression('Affinity constant'~italic(f[C])~'('*µM~NO[3]^{'-'}~L^-1~day^-1*')'))
Plot2=grid.arrange(grobs=Panel, ncol=1, nrow=4)

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) + coord_cartesian(clip="off") + 
    geom_segment(aes(x=SaturXLimL, xend=SaturXLimU, y=SaturYLimL, yend=SaturYLimU, linetype=SigSatur, size=SigSatur), color="grey50", size=1.5) + 
    geom_errorbar(aes(Satur, Value, ymin=ValueL, ymax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_errorbar(aes(Satur, Value, xmin=SaturL, xmax=SaturU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_point(aes(Satur, Value, color=Strain), fill="white", size=7, pch=20) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=25)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=25)) + 
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

Panel=lapply(SplitData, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_y_continuous(expression('Attack rate'~italic(a[B])~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0-(0.18-0.0)*0.08,0.18+(0.18-0.0)*0.08))
Panel[[2]]=Panel[[2]] + scale_y_continuous(expression('Handling time'~italic(h[B])~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0-(9.0-0.0)*0.08,9.0+(9.0-0.0)*0.08))
Panel[[3]]=Panel[[3]] + scale_y_continuous(expression('Particle area'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0-(6.0-0.0)*0.08,6.0+(6.0-0.0)*0.08))
Panel[[4]]=Panel[[4]] + scale_y_continuous(expression('Carbon to nitrogen ratio '~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12-(18-12)*0.08,18+(18-12)*0.08))
Panel[[1]]=Panel[[1]] + annotate("text", label="NS", y=0.18, x=2.8, color="black", size=8, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label="NS", y=9.0, x=2.8, color="black", size=8, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=6.0, x=2.8, color="black", size=8, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=18, x=2.8, color="black", size=8, fontface="italic")
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.19, ymax=0.18+(0.18-0.0)*0.19)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="first", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Competitiveness'), gp=gpar(fontface="bold", fontsize=25), rot=0), xmin=-Inf, xmax=Inf, ymin=0.18+(0.18-0.0)*0.26, ymax=0.18+(0.18-0.0)*0.26)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=25), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=25), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.2,0,0.1),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=25), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0.2,0,0.35),"cm")) + annotation_custom(grob=linesGrob(arrow=arrow(type="closed", ends="last", length=unit(0.4,"cm")), gp=gpar(col="black", fill="black", lty="solid", lwd=2)), xmin=2.8+(2.8-0.09)*0.18, xmax=2.8+(2.8-0.09)*0.18, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0.2,0,0.35),"cm")) + annotation_custom(grob=textGrob(expression('Defense'), gp=gpar(fontface="bold", fontsize=25), rot=270), xmin=2.8+(2.8-0.09)*0.23, xmax=2.8+(2.8-0.09)*0.23, ymin=-Inf, ymax=Inf)
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(vjust=1))
Panel[[4]]=Panel[[4]] + theme(axis.title.x=element_text(face="plain", colour="black", size=25))
Panel[[4]]=Panel[[4]] + xlab(expression('Half-saturation constant'~italic(K[C])~'('*µM~NO[3]^{'-'}~L^-1*')'))
Plot3=grid.arrange(grobs=Panel, ncol=1, nrow=4)

tiff('[Labels] Trait Spaces.tiff', units="in", width=21, height=27.5, res=1000)
grid.arrange(Plot1, Plot2, Plot3, ncol=3, nrow=1)
dev.off()


######################################
### Correlation traits and fitness ###
######################################

MeltData2$Strain=gsub("CR6", "CR7", MeltData2$Strain)
MeltData2$Strain=gsub("CR5", "CR6", MeltData2$Strain)
Labels=c(expression(C[R1]),expression(C[R2]),expression(C[R3]),expression(C[R4]),expression(C[R6]),expression(C[R7]))
SplitData2=split(MeltData2, list(MeltData2$Trait))[c(1:4)]
SplitData3=split(MeltData2, list(MeltData2$Trait))[c(5:7)]

# Simulate regression solutions
DataA=data.frame(
  XGrow=seq(0.00,0.18,length.out=1000),
  YGrow=SplitData2[[1]]$Inter + SplitData2[[1]]$Slope*seq(0.00,0.18,length.out=1000))

DataH=data.frame(
  XGrow=seq(0.00,9.00,length.out=1000),
  YGrow=SplitData2[[2]]$Inter + SplitData2[[2]]$Slope*seq(0.00,9.00,length.out=1000))

DataS=data.frame(
  XGrow=seq(0.00,6.00,length.out=1000),
  YGrow=SplitData2[[3]]$Inter + SplitData2[[3]]$Slope*seq(0.00,6.00,length.out=1000))

DataCN=data.frame(
  XGrow=seq(12,18,length.out=1000),
  YGrow=SplitData2[[4]]$Inter + SplitData2[[4]]$Slope*seq(12,18,length.out=1000))

DataG=data.frame(
  XGrow=seq(0.09,0.36,length.out=1000),
  YGrow=SplitData3[[1]]$Inter + SplitData3[[1]]$Slope*seq(0.09,0.36,length.out=1000))

DataF=data.frame(
  XGrow=seq(0.0,0.6,length.out=1000),
  YGrow=SplitData3[[2]]$Inter + SplitData3[[2]]$Slope*seq(0.0,0.6,length.out=1000))

DataK=data.frame(
  XGrow=seq(0.4,2.8,length.out=1000),
  YGrow=SplitData3[[3]]$Inter + SplitData3[[3]]$Slope*seq(0.4,2.8,length.out=1000))

# Cut the regression solutions
GrowLim=DataA[DataA$YGrow >= -0.05 & DataA$YGrow <= 5.5 & DataA$XGrow >= 0.00 & DataA$XGrow <= 0.18,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataA2=data.frame(Trait="Attack",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataH[DataH$YGrow >= -0.05 & DataH$YGrow <= 5.5 & DataH$XGrow >= 0.00 & DataH$XGrow <= 9.00,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataH2=data.frame(Trait="Handling",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataS[DataS$YGrow >= -0.05 & DataS$YGrow <= 5.5 & DataS$XGrow >= 0.00 & DataS$XGrow <= 6.00,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataS2=data.frame(Trait="Area",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataCN[DataCN$YGrow >= -0.05 & DataCN$YGrow <= 5.5 & DataCN$XGrow >= 12 & DataCN$XGrow <= 18,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataCN2=data.frame(Trait="Stoichio",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataG[DataG$YGrow >= 0.2 & DataG$YGrow <= 0.8 & DataG$XGrow >= 0.09 & DataG$XGrow <= 0.36,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataG2=data.frame(Trait="Grow",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataF[DataF$YGrow >= 0.2 & DataF$YGrow <= 0.8 & DataF$XGrow >= 0.0 & DataF$XGrow <= 0.6,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataF2=data.frame(Trait="Affin",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

GrowLim=DataK[DataK$YGrow >= 0.2 & DataK$YGrow <= 0.8 & DataK$XGrow >= 0.4 & DataK$XGrow <= 2.8,]
GrowYLimL=head(GrowLim,1)[,2]
GrowYLimU=tail(GrowLim,1)[,2]
GrowXLimL=head(GrowLim,1)[,1]
GrowXLimU=tail(GrowLim,1)[,1]
DataH2=data.frame(Trait="Satur",GrowYLimL,GrowYLimU,GrowXLimL,GrowXLimU)

# Include limits in the dataset
DataLim=rbind(DataA2,DataH2,DataS2,DataCN2,DataG2,DataF2,DataH2)
DataLim=DataLim[rep(seq_len(nrow(DataLim)),each=6),]
MeltData2=cbind(MeltData2,DataLim[,-1])
SplitData2=split(MeltData2, list(MeltData2$Trait))[c(1:4)]
SplitData3=split(MeltData2, list(MeltData2$Trait))[c(5:7)]


###############################################
### Correlation plots of traits and fitness ###
###############################################

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) +
    geom_segment(aes(x=GrowXLimL, xend=GrowXLimU, y=GrowYLimL, yend=GrowYLimU, linetype=Sig, size=Sig), color="grey50", size=1.5) + 
    geom_errorbar(aes(Value, PredG, ymin=PredGL, ymax=PredGU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_errorbar(aes(Value, PredG, xmin=ValueL, xmax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_point(aes(Value, PredG, color=Strain), fill="white", size=7, pch=20) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_text(face="plain", colour="black", size=24)) + 
    scale_y_continuous(labels=sprintf(seq(-0.5,5.5,by=2.0), fmt="%.1f"), breaks=seq(-0.5,5.5,by=2.0), limits=c(-0.5,5.912)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_linetype_manual(values=c("Yes"="solid","No"=NA)) + scale_size_manual(values=c("Yes"=1,"No"=1.2)) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=2, nrow=2) +
    theme(legend.position="none")
}

Panel=lapply(SplitData2, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_x_continuous(expression('Attack rate'~italic(a[B])~'('*10^-6~mL~sec^-1*')'), labels=sprintf(seq(0,1.8,by=0.6), fmt="%.1f"), breaks=c(seq(0,0.18,by=0.06)), limits=c(0,0.18))
Panel[[2]]=Panel[[2]] + scale_x_continuous(expression('Handling time'~italic(h[B])~'('*sec*')'), labels=sprintf(seq(0,9.0,by=3.0), fmt="%.1f"), breaks=c(seq(0,9.0,by=3.0)), limits=c(0,9.0))
Panel[[3]]=Panel[[3]] + scale_x_continuous(expression('Particle area'~italic(s[C])~'('*10^2~µm^2*')'), labels=sprintf(seq(0,6.0,by=2.0), fmt="%.1f"), breaks=c(seq(0,6.0,by=2.0)), limits=c(0,6.0))
Panel[[4]]=Panel[[4]] + scale_x_continuous(expression('Carbon to nitrogen ratio'~italic('C:N'[C])), labels=sprintf(seq(12,18,by=2.0), fmt="%.0f"), breaks=c(seq(12,18,by=2.0)), limits=c(12,18))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=5.5, x=0.18, color="black", size=8, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label="NS", y=5.5, x=9.0, color="black", size=8, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=5.5, x=6.0, color="black", size=8, fontface="italic")
Panel[[4]]=Panel[[4]] + annotate("text", label="NS", y=5.5, x=18, color="black", size=8, fontface="italic")
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R1]), x=-0.05, y=0.1800, color="mediumpurple3", size=7) + annotate("text", label=expression(C[R2]), x=-0.05, y=0.1656, color="cornflowerblue", size=7)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R3]), x=-0.05, y=0.1512, color="chartreuse3", size=7) + annotate("text", label=expression(C[R4]), x=-0.05, y=0.1368, color="gold2", size=7)
Panel[[1]]=Panel[[1]] + annotate("text", label=expression(C[R6]), x=-0.05, y=0.1224, color="darkorange1", size=7) + annotate("text", label=expression(C[R7]), x=-0.05, y=0.1080, color="tomato2", size=7)
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.5,0,0.1),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.5,0,0.1),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.5,0,0.1),"cm"))
Panel[[4]]=Panel[[4]] + theme(plot.margin=unit(c(0.5,0.5,0,0.1),"cm"))
Yaxis=textGrob(expression('Fitness growth rate'~italic(r[B])~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=24), rot=90)
Plot1=grid.arrange(grobs=Panel, left=Yaxis, ncol=2, nrow=2)

PlotFunc=function(x) {
  ggplot(x, aes(group=Strain)) +
    geom_segment(aes(x=GrowXLimL, xend=GrowXLimU, y=GrowYLimL, yend=GrowYLimU, linetype=Sig, size=Sig), color="grey50", size=1.5) + 
    geom_errorbar(aes(Value, PreyG, ymin=PreyGL, ymax=PreyGU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_errorbar(aes(Value, PreyG, xmin=ValueL, xmax=ValueU, color=Strain), linetype="solid", alpha=0.7, size=1.5, width=0) +
    geom_point(aes(Value, PreyG, color=Strain), fill="white", size=7, pch=20) + 
    theme(axis.text.y=element_text(face="plain", colour="black", size=24)) +  
    theme(axis.text.x=element_text(face="plain", colour="black", size=24)) + 
    theme(axis.title.y=element_blank()) +
    theme(axis.title.x=element_text(face="plain", colour="black", size=24)) + 
    scale_y_continuous(labels=sprintf(seq(0.2,0.8,by=0.2), fmt="%.1f"), breaks=seq(0.2,0.8,by=0.2), limits=c(0.2,0.86)) +
    theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
    scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR6"="darkorange1","CR7"="tomato2")) +
    theme(strip.background=element_blank(), strip.text.x=element_blank()) +
    facet_wrap(~Trait, ncol=2, nrow=2) +
    theme(legend.position="none")
}

Panel=lapply(SplitData3, PlotFunc)
Panel[[1]]=Panel[[1]] + scale_x_continuous(expression('Maximum growth rate'~italic(r[C~max])~'('*day^-1*')'), labels=sprintf(seq(0.09,0.36,by=0.09), fmt="%.2f"), breaks=seq(0.09,0.36,by=0.09), limits=c(0.09,0.36))
Panel[[2]]=Panel[[2]] + scale_x_continuous(expression('Affinity constant'~italic(f[C])~'('*µM~NO[3]^{'-'}~L^-1~day^-1*')'), labels=sprintf(seq(0,0.6,by=0.2), fmt="%.1f"), breaks=seq(0,0.6,by=0.2), limits=c(0,0.6))
Panel[[3]]=Panel[[3]] + scale_x_continuous(expression('Half-saturation constant'~italic(K[C])~'('*µM~NO[3]^{'-'}~L^-1*')'), labels=sprintf(seq(0.4,2.8,by=0.8), fmt="%.1f"), breaks=seq(0.4,2.8,by=0.8), limits=c(0.4,2.8))
Panel[[1]]=Panel[[1]] + annotate("text", label=NA, y=1.6, x=0.18, color="black", size=8, fontface="italic")
Panel[[2]]=Panel[[2]] + annotate("text", label=NA, y=1.6, x=9.0, color="black", size=8, fontface="italic")
Panel[[3]]=Panel[[3]] + annotate("text", label="NS", y=1.6, x=6.0, color="black", size=8, fontface="italic")
Panel[[1]]=Panel[[1]] + theme(plot.margin=unit(c(0.5,0.5,0,0.40),"cm"))
Panel[[2]]=Panel[[2]] + theme(plot.margin=unit(c(0.5,0.5,0,0.40),"cm"))
Panel[[3]]=Panel[[3]] + theme(plot.margin=unit(c(0.5,0.5,0,0.40),"cm"))
Yaxis=textGrob(expression('Fitness growth rate'~italic(r[C])~'('*day^-1*')'), gp=gpar(fontface="bold", fontsize=24), rot=90)
Plot2=grid.arrange(grobs=Panel, left=Yaxis, ncol=2, nrow=2)

tiff('[Labels] Correlation Traits.tiff', units="in", width=13.7, height=27.5, res=1000)
grid.arrange(Plot1, Plot2, ncol=1, nrow=2)
dev.off()
