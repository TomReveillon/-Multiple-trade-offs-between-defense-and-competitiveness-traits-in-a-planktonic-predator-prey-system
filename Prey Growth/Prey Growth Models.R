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

#######################################################################################
#######################################################################################
##### PREY GROWTH RATE MODELS FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
#######################################################################################
#######################################################################################

# Import the dataset
Data=read.table("Data_AG.txt", h=T, dec=",")
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
Data2=setDT(Data)[, .(MeanDens=mean(Dens)), by=list(Strain,Nitro,Day)]
Data2=as.data.frame(Data2)


##########################
### Growth rate models ###
##########################

# Linear model
SetLi=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=min(xy$y); a=coef(lm(y~x,xy))[2] 
  Out=c(b,a); names(Out)=mCall[c("b","a")]
  return(Out)}
SSline=selfStart(as.formula("~b + a*x"), initial=SetLi, parameters=c("b","a"))

FuncLi=function(Model, Day) {
  Params=coef(Model)
  names(Params)=NULL
  b=Params[1]
  a=Params[2]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,a=a)
  Rates=data.frame(
    Time = Day,
    DensP = b + a * Day,
    AGR = a)
  
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

# Exponential model
SetEx=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  b=min(xy$y); r=coef(lm(log(y)~x,xy))[2]
  Out=c(b,r); names(Out)=mCall[c("b","r")]
  return(Out)}
SSexpo=selfStart(as.formula("~b * exp(r * x)"), initial=SetEx, parameters=c("b","r"))

FuncEx=function(Model,Day) {
  Params=coef(Model)
  names(Params)=NULL
  b=Params[1]
  r=Params[2]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(b=b,r=r)
  Rates=data.frame(
    Time = Day,
    DensP = b * exp(r * Day),
    AGR = r * b * exp(r * Day))
    
  Rates$RGR = r
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

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

# Gompertz model
SetGo=function(Coef){
  K=Coef[1]; b=K/exp(Coef[2]); r=-log(Coef[3])
  Out=c(b,K,r); names(Out)=c("b","K","r") 
  return(Out)}

FuncGo=function(Model,Day) {
  Coef=coef(Model)
  Params=SetGo(Coef)
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
    DensP = K * ((b / K)^exp(-r * Day)),
    AGR = r * K * exp(-r * Day) * log(K / b) * (b / K )^exp(-r * Day))
    
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}

# Mortality model
SetMo=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  K=max(xy$y); m=-0.1; t=5.0; r=0.3
  Out=c(K,m,t,r); names(Out)=mCall[c("K","m","t","r")]
  return(Out)}
SSmort=selfStart(as.formula("~K * exp(m * (x - t) - (m / r) * (1 - exp(-r * (x - t))))"), initial=SetMo, parameters=c("K","m","t","r"))

FuncMo=function(Model,Day) {
  Params=coef(Model)
  names(Params)=NULL
  K=Params[1]
  m=Params[2]
  t=Params[3]
  r=Params[4]
  
  Fit=Model$m$fitted()
  Res=Model$m$resid()
  MSS=sum((Fit-mean(Fit))^2)
  RSS=sum(Res^2)
  R2=MSS/(MSS+RSS)
  AIC=AIC(Model)
  
  Summary=data.frame(R2,AIC)
  Parameters=data.frame(K=K,m=m,t=t,r=r)
  Rates=data.frame(
    Time = Day,
    DensP = K * exp(m * (Day - t) - (m / r) * (1 - exp(-r * (Day - t)))),
    AGR = r)
  
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}


############################
### Fitting growth rates ###
############################

# Split the dataset per strain
SplitData2=split(Data2, list(Data2$Strain))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Day")])
Names=Names[order(Names$Strain),]
Strain=Names$Strain

# Linear model
ModLi=function(x) {
  FitLi=nls(log(MeanDens) ~ SSline(Day, a, b), data=x)
  OutLi=FuncLi(FitLi, x$Day)} 
OutLi=lapply(SplitData2, ModLi)

RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(as.data.frame(do.call("rbind",RateLi)),4)
RateLi=cbind(Strain=Strain[c(1:60)],Nitro="100",Model="Linear",RateLi)
rownames(RateLi)=c()
ParamLi=bind_rows(lapply(OutLi, function (x) x[c("Parameters")]))
ParamLi=round(as.data.frame(do.call("rbind",ParamLi)),4)
rownames(ParamLi)=c()
SummaLi=bind_rows(lapply(OutLi, function (x) x[c("Summary")]))
SummaLi=round(as.data.frame(do.call("rbind",SummaLi)),4)
rownames(SummaLi)=c()
ModelLi=unlist(lapply(OutLi, function (x) x[c("Model")]),recursive=F)

# Exponential model
ModEx=function(x) {
  FitEx=nls(log(MeanDens) ~ SSexpo(Day, r, b), data=x)
  OutEx=FuncEx(FitEx, x$Day)}
OutEx=lapply(SplitData2, ModEx)

RateEx=bind_rows(lapply(OutEx, function (x) x[c("Rates")]))
RateEx=round(as.data.frame(do.call("rbind",RateEx)),4)
RateEx=cbind(Strain=Strain[c(1:60)],Nitro="100",Model="Exponential",RateEx)
rownames(RateEx)=c()
ParamEx=bind_rows(lapply(OutEx, function (x) x[c("Parameters")]))
ParamEx=round(as.data.frame(do.call("rbind",ParamEx)),4)
rownames(ParamEx)=c()
SummaEx=bind_rows(lapply(OutEx, function (x) x[c("Summary")]))
SummaEx=round(as.data.frame(do.call("rbind",SummaEx)),4)
rownames(SummaEx)=c()
ModelEx=unlist(lapply(OutEx, function (x) x[c("Model")]),recursive=F)

# Logistic model
ModLo=function(x) {
  FitLo=nls(log(MeanDens) ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)}
OutLo=lapply(SplitData2, ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain[c(1:60)],Nitro="100",Model="Logistic",RateLo)
rownames(RateLo)=c()
ParamLo=bind_rows(lapply(OutLo, function (x) x[c("Parameters")]))
ParamLo=round(as.data.frame(do.call("rbind",ParamLo)),4)
rownames(ParamLo)=c()
SummaLo=bind_rows(lapply(OutLo, function (x) x[c("Summary")]))
SummaLo=round(as.data.frame(do.call("rbind",SummaLo)),4)
rownames(SummaLo)=c()
ModelLo=unlist(lapply(OutLo, function (x) x[c("Model")]),recursive=F)

# Gompertz model
ModGo=function(x) {
  FitGo=nls(log(MeanDens) ~ SSgompertz(Day, b, K, r), data=x)
  OutGo=FuncGo(FitGo, x$Day)}
OutGo=lapply(SplitData2[-c(2,3)], ModGo)

RateGo=bind_rows(lapply(OutGo, function (x) x[c("Rates")]))
RateGo=round(as.data.frame(do.call("rbind",RateGo)),4)
RateGo=cbind(Strain=Strain[-c(11:30)],Nitro="100",Model="Gompertz",RateGo)
rownames(RateGo)=c()
ParamGo=bind_rows(lapply(OutGo, function (x) x[c("Parameters")]))
ParamGo=round(as.data.frame(do.call("rbind",ParamGo)),4)
rownames(ParamGo)=c()
SummaGo=bind_rows(lapply(OutGo, function (x) x[c("Summary")]))
SummaGo=round(as.data.frame(do.call("rbind",SummaGo)),4)
rownames(SummaGo)=c()
ModelGo=unlist(lapply(OutGo, function (x) x[c("Model")]),recursive=F)

# Mortality model
ModMo=function(x) {
  FitMo=nls(log(MeanDens) ~ SSmort(Day, K, t, m, r), data=x)
  OutMo=FuncMo(FitMo, x$Day)}
OutMo=lapply(SplitData2[-c(1:3,5:6)], ModMo)

RateMo=bind_rows(lapply(OutMo, function (x) x[c("Rates")]))
RateMo=round(as.data.frame(do.call("rbind",RateMo)),4)
RateMo=cbind(Strain=Strain[-c(1:30,41:60)],Nitro="100",Model="Mortality",RateMo)
rownames(RateMo)=c()
ParamMo=bind_rows(lapply(OutMo, function (x) x[c("Parameters")]))
ParamMo=round(as.data.frame(do.call("rbind",ParamMo)),4)
rownames(ParamMo)=c()
SummaMo=bind_rows(lapply(OutMo, function (x) x[c("Summary")]))
SummaMo=round(as.data.frame(do.call("rbind",SummaMo)),4)
rownames(SummaMo)=c()
ModelMo=unlist(lapply(OutMo, function (x) x[c("Model")]),recursive=F)

# Set the predicted dataset
Strain=rep(unique(Data2$Strain), each=101)
Nitro=rep(unique(Data2$Nitro), each=6*101)
DayP=as.numeric(rep(seq(0,10,by=0.1),6))
Data3=data.frame(Strain,Nitro,DayP)

# Linear model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLi)) {DensP[[i]]=predict(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLi)) {DensPSD[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLi)) {DensPLCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLi)) {DensPUCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLi=data.frame(Data3[c(1:606),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Linear")

# Exponential model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEx)) {DensP[[i]]=predict(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEx)) {DensPSD[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEx)) {DensPLCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEx)) {DensPUCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEx=data.frame(Data3[c(1:606),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Exponential")

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predict(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(1:606),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")

# Gompertz model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelGo)) {DensP[[i]]=predict(ModelGo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelGo)) {DensPSD[[i]]=predictNLS(ModelGo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelGo)) {DensPLCI[[i]]=predictNLS(ModelGo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelGo)) {DensPUCI[[i]]=predictNLS(ModelGo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataGo=data.frame(Data3[-c(102:303),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Gompertz")

# Mortality model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelMo)) {DensP[[i]]=predict(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelMo)) {DensPSD[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelMo)) {DensPLCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelMo)) {DensPUCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataMo=data.frame(Data3[-c(1:303,405:606),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Mortality")


####################################
### Comparing growth rate models ###
####################################

# Extract combinations of names
Strain=unique(Data2[,c("Strain")])
Strain1=Strain[-c(2,3)]
Strain2=Strain[-c(1,2,3,5,6)]

# Summaries of models
SummaLi$Strain=Strain; SummaLi$Model="Linear"
SummaEx$Strain=Strain; SummaEx$Model="Exponential"
SummaLo$Strain=Strain; SummaLo$Model="Logistic"
SummaGo$Strain=Strain1; SummaGo$Model="Gompertz"
SummaMo$Strain=Strain2; SummaMo$Model="Mortality"

# Combine AIC coefficients
DataAIC=rbind(SummaLi,SummaEx,SummaLo,SummaGo,SummaMo)
DataAIC=DataAIC[order(DataAIC$Strain),]

# Likelihood ratios
LR1=list(); LR2=list(); LR3=list(); LR4=list()
for (i in 1:length(Strain)) {
  for (j in 1:length(Strain1)) {
    for (k in 1:length(Strain2)) {
      LR1[[i]]=lrtest(ModelEx[[i]],ModelLi[[i]])
      LR2[[i]]=lrtest(ModelLo[[i]],ModelLi[[i]])
      LR3[[j]]=lrtest(ModelGo[[j]],ModelLi[[j]])
      LR4[[k]]=lrtest(ModelMo[[k]],ModelLi[[k]])
    }
  }
}


#######################################
### Plotting predicted growth rates ###
#######################################

# Create a dataset
Data4=rbind(DataLi[,c(1:9)],DataEx[,c(1:9)],DataLo[,c(1:9)],DataGo[,c(1:9)],DataMo[,c(1:9)])
Data4[,c(4:8)]=exp(Data4[,c(4:8)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][Data4[,c(4:8)]<0]=0
Data4=Data4[order(Data4$Strain),]

# Include model selection
Data4$AIC=rep(DataAIC[,2], each=101)
Data5=as.data.frame(Data4 %>% group_by(Strain) %>% slice(which.min(AIC)))
Data4$Selection=ifelse(Data4[,1] %in% Data5[,1] & Data4[,10] %in% Data5[,10], "Accepted", "Rejected")

tiff('Growth Models.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP/10^5)) +
  geom_line(data=subset(Data4, Selection=="Rejected"), aes(linetype=Model, size=Model), color="grey60") +
  geom_line(data=subset(Data4, Selection=="Accepted"), aes(color=Strain, linetype=Model, size=Model)) +
  geom_point(data=Data, aes(Day, Dens/10^5, color=Strain), size=2, pch=16) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("Linear"="dotted","Exponential"="dotted","Logistic"="solid","Gompertz"="dashed","Mortality"="longdash")) +
  scale_size_manual(values=c("Linear"=1,"Exponential"=1.4,"Logistic"=1,"Gompertz"=1,"Mortality"=1)) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
