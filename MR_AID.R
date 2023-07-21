#This file is to do MR-AID using one sample individual data

#load necessary packages
library(trio)
library(ivreg)
library(MASS)
library(ggplot2)
library(hdi)
library(glmnet)
library(dplyr)
library(qqman)
library(screening)
library(ivmodel)
library(genetics)

#functions
#fill NAs of gene data
fill_missing_values <- function(column) {
  G <- as.factor(column)
  na_indices <- which(is.na(G))
  if (length(na_indices) > 0) {
    k <- summary(G)
    q <- as.factor(names(k)[1:(length(k) - 1)])
    probabilities <- k[1:(length(k) - 1)] / sum(k[1:(length(k) - 1)])
    G[na_indices] <- sample(q, length(na_indices), replace = TRUE, prob = probabilities)
    return(as.numeric(as.character(G)))
  }
  return(column)
}
#calculate the weights of selected weak IVs
weightcal=function(weight){
  len=length(weight)
  a=(weight>=0)
  b=abs(weight)/sum(abs(weight))
  c=NULL
  for(i in 1:len){
    if(a[i])c[i]=b[i]
    if(!a[i])c[i]=-1*b[i]
  }
  return(c)
}

#load data
#X: exposure variable. Should be a vector
#Y: outcome variable. Should be a vector
#gene: gene matrix. Should be a data frame. 
#      can has NA's. 
#cova: covariates
#N: number of sample size
#I: number of gene variables

#fill NAs of gene data
gene <- as.data.frame(apply(gene, 2, fill_missing_values))

#sample splitting times
split.time=60 #can change
#results
point.est=matrix(0,nrow=split.time,ncol=4)
colnames(point.est)=c("Estimate","Std.error","t value","P.value")
genemajorselected=list()
geneweakselected=list()
for(j in 1:split.time){
  #split sample
  selected=list()
  unselected=1:N
  n=sample(unselected,size=N/2,replace = FALSE)
  selected[[1]]=n
  selected[[2]]=unselected[-n]
  
  selectIV=list()
  #use SIS to first screening IVs
  for(p in 1:2){
    n=selected[[p]]
    fit=screening(gene[-n,],X[-n],method = "sis")
    selectIV[[p]]=gene[,fit$screen]
  }
  
  #Use LASSO to select IVs and get the estimated weights for weak IV
  weights=list()
  for(p in 1:2){
    n=selected[[p]]
    #select IVs 
    lambda=cv.glmnet(as.matrix(selectIV[[p]][-n,]),X[-n],alpha=1)$lambda.min
    lasso=glmnet(selectIV[[p]][-n,],X[-n],lambda=lambda)
    #coef of all IVs
    a=coef(lasso)[-1,]
    weights[[p]]=a[a!=0]
    #selected IVs
    selectIV[[p]]=selectIV[[p]][,a!=0]
  }
  
  #calculate partial R2 to select major IVs
  partialR2=list()
  for(p in 1:2){
    n=selected[[p]]
    full=lm(X[-n]~.,data=selectIV[[p]][-n,])
    sse_full=sum((full$residuals)^2)
    R2=NULL
    for(q in 1:length(selectIV[[p]][1,])){
      reduced=lm(X[-n]~.,data=selectIV[[p]][-n,-q])
      sse_reduced=sum((reduced$residuals)^2)
      R2[q]=(sse_reduced-sse_full)/sse_reduced
    }
    partialR2[[p]]=R2
  }
  
  selectmajor=list()#SNPs selected as major IVs
  selectweak=list()#SNPs selected as weak IVs
  for(p in 1:2){
    selectmajor[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]>=0.1])
    selectweak[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]<0.1])
    genemajorselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]>=0.1])
    geneweakselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]<0.1])
  }
  
  #combine weak IVs 
  combineMajor=list()
  for(p in 1:2){
    n=selected[[p]]
    weakwei=weights[[p]][partialR2[[p]]<0.1]
    weakwei=weightcal(weakwei)
    combineMajor[[p]]=cbind(selectmajor[[p]][n,],as.matrix(selectweak[[p]][n,])%*%weakwei)
  }
  
  #MR-AID
  residual=NULL
  for(p in 1:2){
    n=selected[[p]]
    fit=lm(X[n]~.,data=as.data.frame(combineMajor[[p]]))
    Fstat[j,p+2]=summary(fit)$fstat[1]
    residual=c(residual,fit$residuals)
  }
  ord=c(selected[[1]],selected[[2]])
  fit2=lm(Y[ord]~residual+X[ord]+cova[ord])
  point.est[j,]=summary(fit2)$coef[3,1:4]
}
#output of estimates at each splitting time
write.table(point.est,file="MR-AID.txt")

#Names and times of SNPs' that selected as major and weak IVs
majortime=NULL
for(l in 1:(2*split.time)){
  num=length(genemajorselected[[l]])
  if(num>0){
    for(ll in 1:num){
      posi=match(genemajorselected[[l]][[ll]],colnames(majortime))
      if(is.na(posi)){
        majortime=cbind(majortime,1)
        colnames(majortime)[length(majortime)]=genemajorselected[[l]][[ll]]
      }else{
        majortime[posi]=majortime[posi]+1
      }
    }
  }
}
weaktime=NULL
for(l in 1:(2*split.time)){
  num=length(geneweakselected[[l]])
  if(num>0){
    for(ll in 1:num){
      posi=match(geneweakselected[[l]][[ll]],colnames(weaktime))
      if(is.na(posi)){
        weaktime=cbind(weaktime,1)
        colnames(weaktime)[length(weaktime)]=geneweakselected[[l]][[ll]]
      }else{
        weaktime[posi]=weaktime[posi]+1
      }
    }
  }
}
alltime=matrix(0,nrow=1,ncol = 3)
colnames(alltime)=c("IVs","major","weak")
sel=NULL
for(i in 1:length(majortime)){
  posi=match(colnames(majortime)[i],colnames(weaktime))
  if(!is.na(posi)){
    sel=c(sel,posi)
    alltime=rbind(alltime,c(majortime[i]+weaktime[posi],majortime[i],weaktime[posi]))
    rownames(alltime)[i+1]=colnames(majortime)[i]
  }
  if(is.na(posi)){
    alltime=rbind(alltime,c(majortime[i],majortime[i],0))
    rownames(alltime)[i+1]=colnames(majortime)[i]
  }
}
q=cbind(weaktime[,-sel],0,weaktime[-sel])
rownames(q)=colnames(weaktime)[-sel]
alltime=rbind(alltime,q)
alltime=alltime[-1,]
#output
write.table(alltime,file = "count time.txt")

#cauchy combination test to summarize the results
T1=mean(tan((0.5-point.est[,4])*pi))
p=0.5-atan(T1)/pi
print(p)#aggregated p-values
print(mean(point.est[,1]))#aggregated estimates

#draw figures of pvalues and estimates of multiple splitting times
par(mfrow=c(1,2))
hist(point.est[,4],breaks = 20,main = NULL,xlab = "p-values of multiple split times")
abline(v = 0.05, col = "red", lty = 2)
legend("topright", legend = "p-value=0.05", lty = 2, col = "red", cex = 0.8, bty = "n")
hist(point.est[,1], breaks=10,main = NULL, xlab = "Estimates of multiple split times")
abline(v = mean(point.est1[,1]), col = "red", lty = 2)
legend("topright", legend = paste0("Estimate=",mean(point.est1[,1])), lty = 2, col = "red", cex = 0.8, bty = "n")
