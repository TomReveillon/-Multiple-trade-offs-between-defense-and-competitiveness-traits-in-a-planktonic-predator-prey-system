setwd("~/LIMNO 2019-2022/Experiments/Predator Growth")

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

###########################################################################################
###########################################################################################
##### PREDATOR GROWTH RATE MODELS FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
###########################################################################################
###########################################################################################

# Import the dataset
Data=read.table("Data_RG.txt", h=T, dec=",")
names(Data)
summary(Data)

# Specify the variables as numeric or factor
Data[,c(2,4)] %<>% mutate_if(is.character,as.numeric)

# Sort the dataset by strain
Data=Data[order(Data$Strain),]

# Calculate mean densities of trials
Data2=setDT(Data)[, .(MeanDens=mean(Count)), by=list(Strain,Day)]
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
    DensP = K * (b * K) / (b + (K - b) * exp(-r * Day)),
    AGR = (r * b * K * (K - b) * exp(-r * Day)) / (b + (K - b)*exp(-r * Day))^2)
    
  Rates$RGR = Rates$AGR / Rates$DensP
  
  Out=list(Model=Model, Summary=Summary, Parameters=Parameters, Rates=Rates)
  return(Out)
}


# Mortality model
SetMo=function(mCall,LHS,data) {
  xy=sortedXyData(mCall[["x"]],LHS,data)
  K=max(xy$y); m=-0.5; t=2.0; r=0.5
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


##################################
### Fitting growth rate models ###
##################################

# Split the dataset
SplitData2=split(Data2, list(Data2$Strain))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Day")])
Names=Names[order(Names$Strain),]
Strain=Names$Strain; Time=Names$Day

# Linear model
ModLi=function(x) {
  FitLi=nls(MeanDens ~ SSline(Day, a, b), data=x)
  OutLi=FuncLi(FitLi, x$Day)} 
OutLi=lapply(SplitData2, ModLi)

RateLi=bind_rows(lapply(OutLi, function (x) x[c("Rates")]))
RateLi=round(as.data.frame(do.call("rbind",RateLi)),4)
RateLi=cbind(Strain=Strain[c(1:36)],Nitro="100",Model="Linear",RateLi)
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
  FitEx=nls(MeanDens ~ SSexpo(Day, r, b), data=x)
  OutEx=FuncEx(FitEx, x$Day)}
OutEx=lapply(SplitData2, ModEx)

RateEx=bind_rows(lapply(OutEx, function (x) x[c("Rates")]))
RateEx=round(as.data.frame(do.call("rbind",RateEx)),4)
RateEx=cbind(Strain=Strain[c(1:36)],Nitro="100",Model="Exponential",RateEx)
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
  FitLo=nls(MeanDens ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)}
OutLo=lapply(SplitData2[c(2:6)], ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain[c(7:36)],Model="Logistic",RateLo)
rownames(RateLo)=c()
ParamLo=bind_rows(lapply(OutLo, function (x) x[c("Parameters")]))
ParamLo=round(cbind(do.call("rbind",ParamLo)),4)
rownames(ParamLo)=c()
SummaLo=bind_rows(lapply(OutLo, function (x) x[c("Summary")]))
SummaLo=round(cbind(do.call("rbind",SummaLo)),4)
rownames(SummaLo)=c()
ModelLo=unlist(lapply(OutLo, function (x) x[c("Model")]),recursive=F)


# Mortality model
ModMo=function(x) {
  FitMo=nls(MeanDens ~ SSmort(Day, K, t, m, r), data=x)
  OutMo=FuncMo(FitMo, x$Day)}
OutMo=lapply(SplitData2[c(1:4)], ModMo)

RateMo=bind_rows(lapply(OutMo, function (x) x[c("Rates")]))
RateMo=round(as.data.frame(do.call("rbind",RateMo)),4)
RateMo=cbind(Strain=Strain[c(1:24)],Model="Mortality",RateMo)
rownames(RateMo)=c()
ParamMo=bind_rows(lapply(OutMo, function (x) x[c("Parameters")]))
ParamMo=round(cbind(do.call("rbind",ParamMo)),4)
rownames(ParamMo)=c()
SummaMo=bind_rows(lapply(OutMo, function (x) x[c("Summary")]))
SummaMo=round(cbind(do.call("rbind",SummaMo)),4)
rownames(SummaMo)=c()
ModelMo=unlist(lapply(OutMo, function (x) x[c("Model")]),recursive=F)

# Set the predicted dataset
Strain=rep(unique(Data2$Strain), each=51)
DayP=as.numeric(rep(seq(0,5,by=0.1),6))
Data3=data.frame(Strain,DayP)

# Linear model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLi)) {DensP[[i]]=predict(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLi)) {DensPSD[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLi)) {DensPLCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLi)) {DensPUCI[[i]]=predictNLS(ModelLi[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLi=data.frame(Data3[c(1:306),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Linear")

# Exponential model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelEx)) {DensP[[i]]=predict(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelEx)) {DensPSD[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelEx)) {DensPLCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelEx)) {DensPUCI[[i]]=predictNLS(ModelEx[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataEx=data.frame(Data3[c(1:306),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Exponential")

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predict(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(52:306),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")

# Mortality model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelMo)) {DensP[[i]]=predict(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelMo)) {DensPSD[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelMo)) {DensPLCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelMo)) {DensPUCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataMo=data.frame(Data3[c(1:204),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Mortality")


####################################
### Comparing growth rate models ###
####################################

# Extract combinations of names
Strain=unique(Data2[,c("Strain")])
Strain1=Strain[-c(1)]
Strain2=Strain[-c(5,6)]

# Summaries of models
SummaLi$Strain=Strain; SummaLi$Model="Linear"
SummaEx$Strain=Strain; SummaEx$Model="Exponential"
SummaLo$Strain=Strain1; SummaLo$Model="Logistic"
SummaMo$Strain=Strain2; SummaMo$Model="Mortality"

# Combine AIC coefficients
DataAIC=rbind(SummaLi,SummaEx,SummaLo,SummaMo)
DataAIC=DataAIC[order(DataAIC$Strain),]

# Likelihood ratios
LR1=list(); LR2=list(); LR3=list()
for (i in 1:length(Strain)) {
  for (j in 1:length(Strain1)) {
    for (k in 1:length(Strain2)) {
      LR1[[i]]=lrtest(ModelEx[[i]],ModelLi[[i]])
      LR2[[j]]=lrtest(ModelLo[[j]],ModelLi[[j]])
      LR3[[k]]=lrtest(ModelMo[[k]],ModelLi[[k]])
    }
  }
}


#######################################
### Plotting predicted growth rates ###
#######################################

# Create a dataset
Data4=rbind(DataLi[,c(1:8)],DataEx[,c(1:8)],DataLo[,c(1:8)],DataMo[,c(1:8)])
Data4[,c(3:7)]=round(Data4[,c(3:7)],4)
Data4[,c(3:7)][Data4[,c(3:7)]<0]=0
Data4=Data4[order(Data4$Strain),]

# Include model selection
Data4$AIC=rep(DataAIC[,2], each=51)
Data5=as.data.frame(Data4 %>% group_by(Strain) %>% slice(which.min(AIC)))
Data4$Selection=ifelse(Data4[,1] %in% Data5[,1] & Data4[,9] %in% Data5[,9], "Accepted", "Rejected")

tiff('Growth Models.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP)) +
  geom_line(data=subset(Data4, Selection=="Rejected"), aes(linetype=Model, size=Model), color="grey70") +
  geom_line(data=subset(Data4, Selection=="Accepted"), aes(color=Strain, linetype=Model, size=Model)) +
  geom_point(data=Data, aes(Day, Count, color=Strain), size=2, pch=16) +
  ylab(expression(italic('B. calyciflorus')~'density'~'('*rotifers~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,40,by=10), limits=c(0,40)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,5,by=1), limits=c(0,5)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_linetype_manual(values=c("Linear"="dotted","Exponential"="dotted","Logistic"="solid","Gompertz"="dashed","Mortality"="longdash")) +
  scale_size_manual(values=c("Linear"=1,"Exponential"=1.4,"Logistic"=1,"Gompertz"=1,"Mortality"=1)) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
