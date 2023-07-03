#different sample size, different splitting time
library(latex2exp)
library(ivreg)
library(MASS)
library(ggplot2)
library(hdi)
library(glmnet)
library(dplyr)
library(qqman)
library(screening)
#Y outcome, D exposure
#Z instruments, X covariates
J=300 #number of SNPs
J1=5 #useful SNPs
N=c(500,1000,3000) #sample size
mul=c(100,100)
#mul=c(80,60)
MAF=0.3 #minor allele frequency
Times=1000
rho=0.1
squareh=c(0.15,0.2,0.3)#extremely weak, weak, moderate
ratio=matrix(0,nrow=5,ncol=2)
ratio[,1]=c(0.4,0.4,0.1,0.05,0.05)
ratio[,2]=c(0.2,0.2,0.2,0.2,0.2)
beta=0.25 #causal effect
k=2
set.seed(1)
rr=1
for(ii in 1:1){
  for(hh in 2:2){
    file=paste0("Ty_multiple_","h_",squareh[hh],"N_",N[ii],".txt")
    data=read.table(file)
    random=1:mul[ii]
    p=matrix(0,nrow = Times,ncol=length(random))
    meanbeta=matrix(0,nrow = Times,ncol=length(random))
    set.seed(1)
    for(ran in 1:length(random)){
      Data=matrix(0,nrow=random[ran]*Times,ncol=2)
      for(i in 1:Times){
        sam=sample(1:mul[ii],size=random[ran])
        Data[(1:random[ran])+(i-1)*random[ran],]=as.matrix(data[(1:mul[ii])+(i-1)*mul[ii],][sam,])
      }
      Data=as.data.frame(Data)
      for(i in 1:Times){
        p1=Data[(1:random[ran])+(i-1)*random[ran],]
        T1=mean(tan((0.5-p1[,2])*pi))
        p[i,ran]=0.5-atan(T1)/pi
        meanbeta[i,ran]=mean(p1[,1])
      }
    }
    write.table(meanbeta,file=paste0("bresmul","h_",squareh[hh],"N_",N[ii],".txt"),row.names = FALSE,col.names = FALSE)
    write.table(p,file=paste0("presmul","h_",squareh[hh],"N_",N[ii],".txt"),row.names = FALSE,col.names = FALSE)
     }
}


hh=3
file=paste0("multiple_","h_",squareh[hh],"N_",N[ii],".txt")
data=read.table(file)
random=1:mul[ii]
p=matrix(0,nrow = Times,ncol=length(random))
meanbeta=matrix(0,nrow = Times,ncol=length(random))
set.seed(1)
for(ran in 1:length(random)){
  Data=matrix(0,nrow=random[ran]*Times,ncol=2)
  for(i in 1:Times){
    sam=sample(1:mul[ii],size=random[ran])
    Data[(1:random[ran])+(i-1)*random[ran],]=as.matrix(data[(1:mul[ii])+(i-1)*mul[ii],][sam,])
  }
  Data=as.data.frame(Data)
  for(i in 1:Times){
    p1=Data[(1:random[ran])+(i-1)*random[ran],]
    T1=mean(tan((0.5-p1[,2])*pi))
    p[i,ran]=0.5-atan(T1)/pi
    meanbeta[i,ran]=mean(p1[,1])
  }
}
power=NULL
for(i in 1:length(p[1,])){
  power[i]=sum(p[,i]<0.05)/length(p[,1])
}
split.time=seq(1,100,1)
df=data.frame(x=split.time,y=power[split.time])
plot(df,type = "l",xlab = "split times",ylab="Power")

boxplot(meanbeta)

abline(
  h=0.08,col="red"
)
par(mar=c(5,5,2,2))

plot(colMeans(meanbeta),type = "l",xlab = "split times",ylab=TeX("\\hat{\\beta}"))
plot(colMedians(meanbeta),type = "l",xlab = "split times",ylab=TeX("\\hat{\\beta}"))

plot(apply(meanbeta, 2, min),type = "l",xlab = "split times",ylab=TeX("min(\\hat{\\beta})"))
plot(apply(meanbeta, 2, max),type = "l",xlab = "split times",ylab=TeX("max(\\hat{\\beta})"))
