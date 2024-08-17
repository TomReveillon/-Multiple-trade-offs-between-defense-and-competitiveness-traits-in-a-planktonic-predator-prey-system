setwd("~/LIMNO 2019-2023/Experiments/Predator Ingestion")

rm(list=ls())

library(bbmle)
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

# Import the dataset
DataI=read.table("Data_FRI.txt", h=T, dec=",")
names(DataI)
summary(DataI)

# Specify the variables as numeric or factor
DataI[,c(2,7)] %<>% mutate_if(is.numeric,as.character())
DataI[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Calculate the initial densities
DataI$IDens=round(DataI$Cells*DataI$Volu*DataI$Site*DataI$Dilu*DataI$Cove,0)
DataI=ddply(DataI,.(Strain,Conc), summarize, IDens=round(mean(IDens),0))

# Import the dataset
DataF=read.table("Data_FRF.txt", h=T, dec=",")
names(DataF)
summary(DataF)

# Specify the variables as numeric numbers
DataF[,c(2,7)] %<>% mutate_if(is.numeric,as.character())
DataF[,c(3:6)] %<>% mutate_if(is.character,as.numeric)

# Calculate the final densities
DataF$FDens=round(DataF$Cells*DataF$Volu*DataF$Site*DataF$Dilu*DataF$Cove,0)
DataF=ddply(DataF,.(Strain,Conc,Trial), summarize, FDens=round(mean(FDens),0))

# Create a complete dataset
Data=data.frame(Strain=DataF[,1], Conc=DataF[,2], Trial=DataF[,3], Time=rep(8/24,198), Pred=rep(4,198), IDens=rep(DataI[,3], each=3), FDens=DataF[,4])

# Calculate densities depletion
Data$DDens=(Data$IDens-Data$FDens)
Data$DDens[Data$DDens<0]=0

# Convert density into ingestion rate
Data$Inges=Data$DDens/Data$Time

# Calculate ingestion rate per rotifer
Data$Inges=Data$Inges/4

# Correct ingestion rate per volume
Data$Inges=Data$Inges/5

# Replace negative values by 0 values
Data$Inges=round(Data$Inges,4)
Data$Inges[Data$Inges<0]=0
Data[is.na(Data)]=0


#######################################
### Ordinary differential equations ###
#######################################

# Rescale densities
Data[,c(6:9)]=Data[,c(6:9)]/10^5

# Split the dataset
SplitData=split(Data, list(Data$Strain))

# Extract combinations of names
Strain=lapply(SplitData, function(x) {unique(x$Strain)})
Strain=as.character(unlist(Strain))
Param=unique(c("a","h","sigma"))

# Functional response function
Inges=function(t, x, parms){
  with(as.list(parms),{
    dA = (-a*x[1]/(1 + a*h*x[1]))*P
    return(list(c(dA)))
  })
}

# Densities depletion function
DDens=c()
DensEaten=function(IDens, a, h, P, Time, steps=100) {
  for (i in 1:length(IDens)){
    DDens[i] = IDens[i] - lsoda(y=IDens[i], times=seq(0,Time[i],length=steps), func=Inges, parms=c(a=a, h=h, P=P[i]))[100,2]
  }
  return(DDens)
}

# Maximum likelihood function
Likelihood=function(DDens, IDens, a, h, sigma, P, Time, steps=100){
  if(a <= 0 || h <= 0 || sigma <= 0) return(Inf)
  Func=DensEaten(IDens=IDens, a=a, h=h, P=P, Time=Time, steps=steps)
  LR=-1*sum(dnorm(x=log(DDens+1), mean=log(Func+1), sd=sigma, log=T))
  return(LR)
}

# Fitting the model
FuncHII=function(x) {
  ModHII=mle2(Likelihood, start=list(a=0.1, h=5, sigma=1), control=list(maxit=1000), data=list(IDens=x$IDens, DDens=x$DDens, P=x$Pred, Time=x$Time))}
OutHII=lapply(SplitData, FuncHII)

CoefHII=lapply(OutHII, summary)
CoefHII=lapply(CoefHII, coef)
CoefHII=round(as.data.frame(do.call("rbind",CoefHII)),4)
CoefHII=cbind(Strain=rep(Strain[c(1:6)], each=3),Param=rep(Param,6),Value=CoefHII[,c(1)],Error=CoefHII[,c(2)])
CoefHII=as.data.frame(CoefHII)
rownames(CoefHII)=c()


##########################################
### Fitting functional response models ###
##########################################

# Set the predicted datasets
Strain=rep(Strain, each=150)
IDensP=as.numeric(rep(seq(0.1,15,by=0.1),6))
Time=round(rep(8/24,150*6),4); Pred=rep(4,150*6)
Data2=data.frame(Strain,IDensP,Time,Pred)
Data3=data.frame(Strain,IDensP,Time,Pred)
Data4=data.frame(Strain,IDensP,Time,Pred)

# Specify the variables as numeric numbers
CoefHII[,c(3:4)] %<>% mutate_if(is.character,as.numeric)
CoefHII=CoefHII %>% arrange(factor(Strain, levels=c("CR1","CR2","CR3","CR4","CR5","CR6")))
CoefHII$Error[CoefHII$Error=="NaN"]=0

# Include the coefficients
Data2$a=rep(subset(CoefHII, Param=="a")$Value, each=150)
Data3$a=rep(subset(CoefHII, Param=="a")$Value-subset(CoefHII, Param=="a")$Error, each=150)
Data4$a=rep(subset(CoefHII, Param=="a")$Value+subset(CoefHII, Param=="a")$Error, each=150)
Data2$h=rep(subset(CoefHII, Param=="h")$Value, each=150)
Data3$h=rep(subset(CoefHII, Param=="h")$Value+subset(CoefHII, Param=="h")$Error, each=150)
Data4$h=rep(subset(CoefHII, Param=="h")$Value-subset(CoefHII, Param=="h")$Error, each=150)

# Correct the coefficients
Data2[,c(5:6)][Data2[,c(5:6)] < 0]=0
Data3[,c(5:6)][Data3[,c(5:6)] < 0]=0
Data4[,c(5:6)][Data4[,c(5:6)] < 0]=0

# Split the dataset
SplitData2=split(Data2, list(Data2$Strain))
SplitData3=split(Data3, list(Data3$Strain))
SplitData4=split(Data4, list(Data4$Strain))

# Holling II model
Func2HII=function(x) {
  Mod2HII=DensEaten(IDens=x$IDensP, a=x$a[1], h=x$h[1], Time=rep(x$Time[1],length(IDensP)), P=rep(x$Pred[1],length(IDensP)))}
Out2HII=lapply(SplitData2, Func2HII)
Rate2HII=round(as.data.frame(bind_rows(Out2HII)),4)
Rate2HII=stack(Rate2HII[1:6])[,c(2:1)]
Out3HII=lapply(SplitData3, Func2HII)
Rate3HII=round(as.data.frame(bind_rows(Out3HII)),4)
Rate3HII=stack(Rate3HII[1:6])[,c(2:1)]
Out4HII=lapply(SplitData4, Func2HII)
Rate4HII=round(as.data.frame(bind_rows(Out4HII)),4)
Rate4HII=stack(Rate4HII[1:6])[,c(2:1)]


###############################################
### Plotting predicted functional responses ###
###############################################

# Create a dataset
Data5=data.frame(Strain=Data2[,1], IDensP=Data2[,2], IngesP=Rate2HII[,2], IngesPLSD=Rate3HII[,2], IngesPUSD=Rate4HII[,2])
Data5[,c(3:5)]=round(Data5[,c(3:5)],4)
Data5[,c(3:5)][Data5[,c(3:5)]<0]=0

# Correct ingestion rate per rotifer
Data5[,c(3:5)]=Data5[,c(3:5)]/4
Data5[,c(3:5)]=round(Data5[,c(3:5)],4)

# Include experimental points
Data6=read.table("Data_FRP.txt", h=T, dec=",")
Data6[,c(2:7)] %<>% mutate_if(is.character,as.numeric)

# Correct standard errors and confidence intervals
Data5$IngesPLSD=ifelse(Data5$IngesPLSD < Data5$IngesP*0.80, Data5$IngesP*0.80, Data5$IngesPLSD)
Data5$IngesPUSD=ifelse(Data5$IngesPUSD > Data5$IngesP*1.20, Data5$IngesP*1.20, Data5$IngesPUSD)

# Export the dataset
Data5[,c(3:5)]=replace(Data5[,c(3:5)],Data5[,c(3:5)]<0,0)
write.table(Data5, file="Data_FRPODE.txt", sep="\t", row.names=F)

tiff('Functional Responses ODE.tiff', units="in", width=8, height=8, res=1000)
ggplot(Data5, aes(IDensP, IngesP, group=Strain)) +
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_line(data=Data6, aes(color=Strain), linetype="dotted", size=1) +
  geom_point(data=Data, aes(IDens, Inges, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(legend.position="none")
dev.off()

tiff('Functional Responses ODE Intervals.tiff', units="in", width=15, height=8, res=1000)
ggplot(Data5, aes(IDensP, IngesP, group=Strain)) +
  geom_line(aes(color=Strain), linetype="solid", size=1) +
  geom_line(data=Data6, aes(color=Strain), linetype="dotted", size=1) +
  geom_ribbon(aes(ymin=IngesPLSD, ymax=IngesPUSD, fill=Strain, color=Strain), linetype="solid", size=0.5, alpha=0.3) +
  geom_point(data=Data, aes(IDens, Inges, color=Strain), size=2, pch=16, position=position_jitter(h=0, w=0.2)) +
  ylab(expression(italic('B. calyciflorus')~'ingestion rate'~'('*cells~sec^-1~rotifer^-1*')')) +
  xlab(expression(italic('C. reinhardtii')~'density'~'('*10^5~cells~mL^-1*')')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=18)) + 
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,1.2,by=0.3), limits=c(0,1.2)) +
  scale_x_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(0,15,by=3), limits=c(0,15)) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  scale_fill_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  scale_color_manual(values=c("CR1"="mediumpurple3","CR2"="cornflowerblue","CR3"="chartreuse3","CR4"="gold2","CR5"="darkorange1","CR6"="tomato2")) +
  theme(strip.text.x=element_blank()) +
  facet_wrap(~Strain, scales="free", ncol=3, nrow=2) +
  theme(legend.position="none")
dev.off()
