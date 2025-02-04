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

################################################################################
################################################################################
##### PREY GROWTH RATE FOR PREDATOR-PREY SYSTEM WITH DIFFERENT CLONE TYPES #####
################################################################################
################################################################################

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


#########################################
### Plots of replicated growth curves ###
#########################################

ggplot(subset(Data, Strain=="CR1"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="mediumpurple3", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR2"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="cornflowerblue", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR3"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="chartreuse3", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR4"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="gold2", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR5"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="darkorange1", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
  theme(axis.line=element_line(colour="black", size=0.7, linetype="solid")) +
  theme(panel.background=element_rect(fill="white", colour="white")) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(legend.position="none")

ggplot(subset(Data, Strain=="CR6"), aes(Day, Dens/10^5, group=Strain)) + 
  geom_smooth(method="loess", color="tomato2", size=1, se=F) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,8,by=1), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(1,10,by=1), limits=c(1,10)) +
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


##################################
### Fitting growth rate models ###
##################################

# Split the dataset
SplitData2=split(Data2, list(Data2$Strain))

# Extract combinations of names
Names=unique(Data2[,c("Strain","Day")])
Names=Names[order(Names$Strain),]
Strain=Names$Strain

# Logistic model
ModLo=function(x) {
  FitLo=nls(log(MeanDens) ~ SSlogis(Day, b, K, r), data=x)
  OutLo=FuncLo(FitLo, x$Day)}
OutLo=lapply(SplitData2[c(1:3,5:6)], ModLo)

RateLo=bind_rows(lapply(OutLo, function (x) x[c("Rates")]))
RateLo=round(as.data.frame(do.call("rbind",RateLo)),4)
RateLo=cbind(Strain=Strain[c(1:30,41:60)],Nitro="100",Model="Logistic",RateLo)
rownames(RateLo)=c()
ParamLo=bind_rows(lapply(OutLo, function (x) x[c("Parameters")]))
ParamLo=round(as.data.frame(do.call("rbind",ParamLo)),4)
rownames(ParamLo)=c()
SummaLo=bind_rows(lapply(OutLo, function (x) x[c("Summary")]))
SummaLo=round(as.data.frame(do.call("rbind",SummaLo)),4)
rownames(SummaLo)=c()
ModelLo=unlist(lapply(OutLo, function (x) x[c("Model")]),recursive=F)

# Mortality model
ModMo=function(x) {
  FitMo=nls(log(MeanDens) ~ SSmort(Day, K, t, m, r), data=x)
  OutMo=FuncMo(FitMo, x$Day)}
OutMo=lapply(SplitData2[c(4)], ModMo)

RateMo=bind_rows(lapply(OutMo, function (x) x[c("Rates")]))
RateMo=round(as.data.frame(do.call("rbind",RateMo)),4)
RateMo=cbind(Strain=Strain[c(31:40)],Nitro="100",Model="Mortality",RateMo)
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

# Logistic model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelLo)) {DensP[[i]]=predict(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelLo)) {DensPSD[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelLo)) {DensPLCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelLo)) {DensPUCI[[i]]=predictNLS(ModelLo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataLo=data.frame(Data3[c(1:303,405:606),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Logistic")

# Mortality model
DensP=list(); DensPSD=list(); DensPLCI=list(); DensPUCI=list()
for (i in 1:length(ModelMo)) {DensP[[i]]=predict(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)))}
for (i in 1:length(ModelMo)) {DensPSD[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,3]}
for (i in 1:length(ModelMo)) {DensPLCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,5]}
for (i in 1:length(ModelMo)) {DensPUCI[[i]]=predictNLS(ModelMo[[i]], newdata=data.frame(Day=unique(DayP)), interval="confidence", alpha=0.05, nsim=10000)[[1]][,6]}
DataMo=data.frame(Data3[c(304:404),], DensP=unlist(DensP), DensPLSD=unlist(DensP)-unlist(DensPSD), DensPUSD=unlist(DensP)+unlist(DensPSD), DensPLCI=unlist(DensPLCI), DensPUCI=unlist(DensPUCI), Model="Mortality")


########################################
### Plotting predicted growth curves ###
########################################

# Create a dataset
Data4=rbind(DataLo[,c(1:9)],DataMo[,c(1:9)])
Data4[,c(4:8)]=exp(Data4[,c(4:8)])
Data4[,c(4:8)]=round(Data4[,c(4:8)],4)
Data4[,c(4:8)][Data4[,c(4:8)]<0]=0
Data4=Data4[order(Data4$Strain),]

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLSD))/lead(Data4$DensP)), Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUSD)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI==0, Data4$DensP-Data4$DensP*((lead(Data4$DensP)-lead(Data4$DensPLCI))/lead(Data4$DensP)), Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI==0, Data4$DensP+Data4$DensP*((lead(Data4$DensPUCI)-lead(Data4$DensP))/lead(Data4$DensP)), Data4$DensPUCI)

# Correct standard errors and confidence intervals
Data4$DensPLSD=ifelse(Data4$DensPLSD < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLSD)
Data4$DensPUSD=ifelse(Data4$DensPUSD > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUSD)
Data4$DensPLCI=ifelse(Data4$DensPLCI < Data4$DensP*0.80, Data4$DensP*0.80, Data4$DensPLCI)
Data4$DensPUCI=ifelse(Data4$DensPUCI > Data4$DensP*1.20, Data4$DensP*1.20, Data4$DensPUCI)

# Export the dataset
Data4[,c(5:8)]=replace(Data4[,c(5:8)],Data4[,c(5:8)]<0,0)
write.table(Data4, file="Data_AGP.txt", sep="\t", row.names=F)

tiff('Growth Curves.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP/10^5, group=Strain)) + 
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_point(data=Data, aes(Day, Dens/10^5, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=1), limits=c(0,10)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()

tiff('Growth Curves Intervals.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data4, aes(DayP, DensP/10^5, group=Strain)) + 
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_ribbon(aes(ymin=DensPLSD/10^5, ymax=DensPUSD/10^5, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
  geom_point(data=Data, aes(Day, Dens/10^5, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  ylab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=1), limits=c(0,10)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()


################################
### Calculating growth rates ###
################################

# Split the dataset
SplitData4=split(Data4, list(Data4$Strain))

# Subset important columns
SplitData5=lapply(SplitData4, "[", c("DensP","DayP"))

# Consecutive per capita growth rates function
FuncGR=function(x) {
  FitGR=dynlm(formula=log(DensP)~L(DayP), data=as.data.frame(x))
  OutGR=c(GrowP=summary(FitGR)$coef[2,1], GrowPSD=summary(FitGR)$coef[2,2])}

# Consecutive per capita growth rates model
ModelGR=function(y) {rollapply(y, width=3, FUN=FuncGR, by.column=F)}

# Calculate consecutive per capita growth rates
OutGR=lapply(SplitData5, ModelGR)
RateGR=round(as.data.frame(do.call("rbind",OutGR)),4)
RateGR$GrowPLSD=RateGR[,1]-RateGR[,2]
RateGR$GrowPUSD=RateGR[,1]+RateGR[,2]

# Select important columns
Strain=c(do.call("rbind",lapply(SplitData4, "[", c("Strain"))))[[1]]
Nitro=c(do.call("rbind",lapply(SplitData4, "[", c("Nitro"))))[[1]]
DayP=c(do.call("rbind",lapply(SplitData4, "[", c("DayP"))))[[1]]
Model=c(do.call("rbind",lapply(SplitData4, "[", c("Model"))))[[1]]

# Create a dataset
DataGR=subset(data.frame(Strain,Nitro,DayP,Model), !DayP %in% c("9.9","10"))
DataGR=data.frame(DataGR,RateGR[,c(1,3,4)]); DataGR=DataGR[,c(1,2,3,5,6,7,4)]


#######################################
### Plotting predicted growth rates ###
#######################################

# Create a dataset
Data5=DataGR[,c(1:7)]
Data5[,c(4:6)]=round(Data5[,c(4:6)],4)
Data5=Data5[order(Data5$Strain),]

tiff('Growth Rates.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(DayP, GrowP, group=Strain)) + 
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  ylab(expression(italic('C. reinhardtii')~'growth rate'~'('*day^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,8,by=2), limits=c(0,8)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=1), limits=c(0,10)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()
