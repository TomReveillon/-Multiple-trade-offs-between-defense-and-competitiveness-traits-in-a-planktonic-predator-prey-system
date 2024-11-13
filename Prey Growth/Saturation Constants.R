setwd("~/LIMNO 2019-2023/Experiments/Prey Growth")

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

###############################################################################
###############################################################################
##### HALF-SATURATION FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
###############################################################################
###############################################################################

# Import the dataset
Data=read.table("Data_HS.txt", h=T, dec=",")
names(Data)
summary(Data)

# Specify the variables as numeric or factor
Data[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Sort the dataset by strain and nitrogen concentrations
Data=Data[order(Data$Strain,Data$Nitro),]
Data$Nitro=factor(Data$Nitro, levels=unique(Data$Nitro))

# Calculate the densities
Data$Dens=Data$Cells*Data$Volu*Data$Site*Data$Dilu*Data$Cove

# Calculate mean densities
Data2=setDT(na.omit(Data))[, .(MeanDens=mean(Dens)), by=list(Strain,Nitro,Day)]
Data2=as.data.frame(Data2)


#########################################
### Plots of replicated growth curves ###
#########################################

ggplot(subset(Data, Strain=="CR1"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="mediumpurple3", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR2"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="cornflowerblue", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR3"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="chartreuse3", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR4"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="gold2", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR5"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="darkorange1", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR6"), aes(Day, Dens/10^5, group=Nitro)) + 
  geom_smooth(method="loess", color="tomato2", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,8,by=1), limits=c(1,8)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")


##########################
### Growth rate models ###
##########################

# Logistic model
SetLo=function(Coef){
  K=Coef[1]; r=1/(Coef[3]); b=K/(1 + exp(Coef[2]/Coef[3]))
  Out=c(b,K,r); names(Out)=c("b","K","r")
  return(Out)}

FuncLo=function(Model,Day) {
  Coef=coef(Model)
  Params=SetLo(Coef)
  b=Params[1] 
  K=Params[2] 
  r=Params[3]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,K=K,r=r)
  Rates=data.frame(
    Time = Day,
    DensP = (b * K) / (b + (K - b) * exp(-r * Day)),
    AGR = (r * b * K * (K - b) * exp(-r * Day)) / (b + (K - b)*exp(-r * Day))^2)
  
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}


##################################
### Fitting growth rate models ###
##################################

# Remove outlier densities
Data2=subset(Data2, !Strain=="CR1"|!Nitro=="20"|!Day=="5")
Data2=subset(Data2, !Strain=="CR1"|!Nitro=="30"|!Day=="5")

# Split the dataset
SplitData2=split(Data2, list(Data2$Nitro,Data2$Strain))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Nitro","Day")])
Names=Names[order(Names$Strain,Names$Nitro),]
Strain=Names$Strain; Nitro=Names$Nitro

# Logistic model
ModLo=function(x) {
  FitLo=nls(log(MeanDens) ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)}
OutLo=lapply(SplitData2, ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain,Nitro=Nitro,Model="Logistic",RateLo)
rownames(RateLo)=c()
ParamLo=bind_rows(lapply(OutLo, function (x) x[c("Parameters")]))
ParamLo=round(as.data.frame(do.call("rbind",ParamLo)),4)
rownames(ParamLo)=c()
SummaLo=bind_rows(lapply(OutLo, function (x) x[c("Summary")]))
SummaLo=round(as.data.frame(do.call("rbind",SummaLo)),4)
rownames(SummaLo)=c()
ModelLo=unlist(lapply(OutLo, function (x) x[c("Model")]),recursive=F)

# Set the predicted dataset
Strain=rep(unique(Data2$Strain), each=8*51)
Nitro=rep(rep(unique(Data2$Nitro), each=51),6)
DayP=as.numeric(rep(seq(0,10,by=0.2),6*8))
Data3=data.frame(Strain,Nitro,DayP)

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,1]}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(1:2448),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")


################################
### Calculating growth rates ###
################################

# Create a dataset
Data4=rbind(DataLo[,c(1:9)])
Data4[,c(4:8)]=exp(Data4[,c(4:8)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],0)
Data4[,c(4:8)][Data4[,c(4:8)]<0]=0
Data4=Data4[order(Data4$Strain,Data4$Nitro),]

# Select the time period
Data4=as.data.frame(Data4 %>% group_by(Strain,Nitro) %>% dplyr::slice(11:n()))

# Extract combinations of names
Names=unique(Data4[,c("Strain","Nitro")])
Names=Names[order(Names$Nitro,Names$Strain),]
Strain=Names$Strain; Nitro=Names$Nitro

# Split the dataset
SplitData4=split(Data4, list(Data4$Strain,Data4$Nitro))

# Mean per capita growth rate models
ModelGR=function(x) {summary(lm(log(DensP+1)~DayP, data=x))}
OutGR=lapply(SplitData4, ModelGR)

# Calculate mean per capita growth rates
ModelGR2=function(x) {c(MeanGR=coef(x)[2],MeanGRLSD=coef(x)[2]-coef(x)[4],MeanGRUSD=coef(x)[2]+coef(x)[4])}
OutGR2=round(as.data.frame(do.call("rbind",lapply(OutGR, ModelGR2))),4)
OutGR2=cbind(Strain=Strain,Nitro=Nitro,OutGR2)
rownames(OutGR2)=c()

# Subset important columns
SplitData5=lapply(SplitData4, "[", c("DensP","DayP"))

# Consecutive per capita growth rates function
FuncGR=function(x) {
  FitGR=dynlm(formula=log(DensP)~L(DayP), data=as.data.frame(x))
  OutGR=summary(FitGR)$coef}

# Consecutive per capita growth rates model
ModelGR=function(y) {rollapply(y, width=3, FUN=FuncGR, by.column=F)}

# Calculate consecutive per capita growth rates
OutGR=lapply(SplitData5, ModelGR)
RateGR=round(as.data.frame(do.call("rbind",OutGR))[,c(2,4)],4)
colnames(RateGR)=c("GrowP","GrowPSD")
RateGR$GrowPLSD=RateGR[,1]-RateGR[,2]
RateGR$GrowPUSD=RateGR[,1]+RateGR[,2]

# Select important columns
Strain=c(do.call("rbind",lapply(SplitData4, "[", c("Strain"))))[[1]]
Nitro=c(do.call("rbind",lapply(SplitData4, "[", c("Nitro"))))[[1]]
DayP=c(do.call("rbind",lapply(SplitData4, "[", c("DayP"))))[[1]]
Model=c(do.call("rbind",lapply(SplitData4, "[", c("Model"))))[[1]]

# Create a dataset
Data5=subset(data.frame(Strain,Nitro,DayP,Model), !DayP %in% c("9.8","10"))
Data5=data.frame(Data5,RateGR[,c(1,3,4)]); Data5=Data5[,c(1,2,3,5,6,7,4)]

# Extract combinations of names
Names=unique(Data4[,c("Strain","Nitro")])
Names=Names[order(Names$Nitro,Names$Strain),]
Strain=Names$Strain; Nitro=Names$Nitro

# Find maximum consecutive per capita growth rates
OutGR3=data.frame(Strain=Strain, Nitro=Nitro)
OutGR3$GrowCM=as.data.frame(setDT(Data5)[, .SD[which.max(GrowP)], by=list(Strain,Nitro)])[,4]
OutGR3$GrowCMLSD=as.data.frame(setDT(Data5)[, .SD[which.max(GrowPLSD)], by=list(Strain,Nitro)])[,5]
OutGR3$GrowCMUSD=as.data.frame(setDT(Data5)[, .SD[which.max(GrowPUSD)], by=list(Strain,Nitro)])[,6]


#####################################################
### Calculating affinity and saturation constants ###
#####################################################

# Create datasets
Data6=cbind(unique(OutGR2[,c(1:2)]), Grow=OutGR2[,3], GrowL=OutGR2[,4], GrowU=OutGR2[,5], GrowC=OutGR3[,3], GrowCL=OutGR3[,4], GrowCU=OutGR3[,5])
Data0=cbind(Strain=unique(OutGR2[,1]), Nitro=rep(0,6), Grow=rep(0,6), GrowL=rep(0,6), GrowU=rep(0,6), GrowC=rep(0,6), GrowCL=rep(0,6), GrowCU=rep(0,6))

# Combine datasets
Data6=rbind(Data0,Data6)
Data6[,c(2:8)] %<>% mutate_if(is.character,as.numeric)
Data6[,c(3:8)]=round(Data6[,c(3:8)],4)
Data6=Data6[order(Data6$Strain,Data6$Nitro),]
Data6=subset(Data6, Nitro < 50)

# Extract combinations of names
Names=unique(Data6[,c("Strain")])
Names=Names[order(Names$Strain),]
Strain=Names$Strain

# Split the dataset
SplitData6=split(Data6, list(Data6$Strain))

# Monod model
FuncH=function(x) {OutH=summary(nls(Grow~(GrowM*Nitro)/(Nitro + K), start=c(GrowM=0.1, K=5), data=x))}
OutH=lapply(SplitData6, FuncH)

# Extract half-saturation constants
FuncH2=function(x) {c(Satur=coef(x)[2],SaturLSD=coef(x)[2]-coef(x)[4],SaturUSD=coef(x)[2]+coef(x)[4])}
OutH2=round(as.data.frame(do.call("rbind",lapply(OutH, FuncH2))),4)
OutH2=cbind(Strain=Strain,OutH2)
rownames(OutH2)=c()

# Extract maximum intrinsic growth rates
FuncH3=function(x) {c(Grow=coef(x)[1],GrowLSD=coef(x)[1]-coef(x)[3],GrowUSD=coef(x)[1]+coef(x)[3])}
OutH3=round(as.data.frame(do.call("rbind",lapply(OutH, FuncH3))),4)
OutH3=cbind(Strain=Strain,OutH3)
rownames(OutH3)=c()

# Extract affinity constants
FuncH4=function(x) {c(Affin=coef(x)[1]/coef(x)[2],AffinLSD=(coef(x)[1]-coef(x)[3])/coef(x)[2],AffinUSD=(coef(x)[1]+coef(x)[3])/coef(x)[2])}
OutH4=round(as.data.frame(do.call("rbind",lapply(OutH, FuncH4))),4)
OutH4=cbind(Strain=Strain,OutH4)
rownames(OutH4)=c()

# Export the dataset
Data7=cbind(OutH3[,c(1:4)],OutH4[,c(2:4)],OutH2[,c(2:4)])
write.table(Data7, file="Data_HSP.txt", sep="\t", row.names=F)

# Monod model
FuncH=function(x) {OutH=nls(Grow~(GrowM*Nitro)/(Nitro + K), start=c(GrowM=0.1, K=5), data=x)}
ModelH=lapply(SplitData6, FuncH)

# Set the predicted dataset
Strain=rep(unique(Data6$Strain), each=151)
NitroP=as.numeric(rep(seq(0,30,by=0.2),6))
Data8=data.frame(Strain,NitroP)

# Saturation model
GrowP=list(); GrowPSD=list(); GrowPLCI=list(); GrowPUCI=list()
for (i in 1:length(ModelH)) {GrowP[[i]]=predict(ModelH[[i]], newdata=data.frame(Nitro=unique(NitroP)))}
for (i in 1:length(ModelH)) {GrowPSD[[i]]=predictNLS(ModelH[[i]], newdata=data.frame(Nitro=unique(NitroP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelH)) {GrowPLCI[[i]]=predictNLS(ModelH[[i]], newdata=data.frame(Nitro=unique(NitroP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelH)) {GrowPUCI[[i]]=predictNLS(ModelH[[i]], newdata=data.frame(Nitro=unique(NitroP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataH=data.frame(Data8[c(1:906),], GrowP=unlist(GrowP), GrowPLSD=unlist(GrowP)-unlist(GrowPSD), GrowPUSD=unlist(GrowP)+unlist(GrowPSD), GrowPLCI=unlist(GrowPLCI), GrowPUCI=unlist(GrowPUCI), Model="Monod")


##################################################
### Plotting saturation and affinity constants ###
##################################################

# Create a dataset
Data9=rbind(DataH[,c(1:8)])
Data9[,c(3:7)]=round(Data9[,c(3:7)],4)
Data9[,c(3:7)][Data9[,c(3:7)]<0]=0
Data9=Data9[order(Data9$Strain),]

# Correct standard errors and confidence intervals
Data9$GrowPLSD=ifelse(Data9$GrowPLSD==0, Data9$GrowP-Data9$GrowP*((lead(Data9$GrowP)-lead(Data9$GrowPLSD))/lead(Data9$GrowP)), Data9$GrowPLSD)
Data9$GrowPUSD=ifelse(Data9$GrowPUSD==0, Data9$GrowP+Data9$GrowP*((lead(Data9$GrowPUSD)-lead(Data9$GrowP))/lead(Data9$GrowP)), Data9$GrowPUSD)
Data9$GrowPLCI=ifelse(Data9$GrowPLCI==0, Data9$GrowP-Data9$GrowP*((lead(Data9$GrowP)-lead(Data9$GrowPLCI))/lead(Data9$GrowP)), Data9$GrowPLCI)
Data9$GrowPUCI=ifelse(Data9$GrowPUCI==0, Data9$GrowP+Data9$GrowP*((lead(Data9$GrowPUCI)-lead(Data9$GrowP))/lead(Data9$GrowP)), Data9$GrowPUCI)

# Correct standard errors and confidence intervals
Data9$GrowPLSD=ifelse(Data9$GrowPLSD < Data9$GrowP*0.80, Data9$GrowP*0.80, Data9$GrowPLSD)
Data9$GrowPUSD=ifelse(Data9$GrowPUSD > Data9$GrowP*1.20, Data9$GrowP*1.20, Data9$GrowPUSD)
Data9$GrowPLCI=ifelse(Data9$GrowPLCI < Data9$GrowP*0.80, Data9$GrowP*0.80, Data9$GrowPLCI)
Data9$GrowPUCI=ifelse(Data9$GrowPUCI > Data9$GrowP*1.20, Data9$GrowP*1.20, Data9$GrowPUCI)

tiff('Saturation Constants.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data9, aes(NitroP, GrowP, group=Strain)) +
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_point(data=Data6, aes(Nitro, Grow, color=Strain), size=2, pch=16) +
  ylab(expression(italic('C. reinhardtii')~'growth rate'~'('*day^-1*')')) + 
  xlab(expression('Nitrate concentration'~'('*µM~NO[3]^{'-'}~L^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,0.4,by=0.1), limits=c(0,0.4)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()
