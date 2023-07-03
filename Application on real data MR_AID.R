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
setwd("C:/Users/shirx/OneDrive/Desktop/Major/Application on real data")
Y=read.table("finaloutcome.txt",header = TRUE)
cova=Y$AGE_.YRS.
Y=Y$CYP2E1_activity
#CYP2E1
X=CYP2E1$CYP2E1_exposure
gene<- read.table("CY_cleaned_gene_data.txt",header = TRUE)
N=nrow(gene)
I=ncol(gene)
#169,2987

a=NULL
for(i in 1:N){
  a[i]=sum(is.na(gene[i,]))
}
summary(a/I)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.005357 0.014731 0.020422 0.023900 0.029461 0.090392 
set.seed(1)
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
gene <- as.data.frame(apply(gene, 2, fill_missing_values))
filter_ld <- function(data1,data2) {
  ld_result <- LD(as.genotype.allele.count(data1), as.genotype.allele.count(data2))
  if (ld_result$`R^2` > r2_threshold) {
    return(1)
  }
  if (ld_result$`R^2` <= r2_threshold) {
    return(0)
  }
}

kk=2
split.time=60
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
r2_threshold <- 0.8
point.est1=matrix(0,nrow=split.time,ncol=4)
point.est2=matrix(0,nrow=split.time,ncol=4)
colnames(point.est1)=c("Estimate","Std.error","t value","P.value")
colnames(point.est2)=c("Estimate","Std.error","t value","P.value")
genemajorselected=list()
geneweakselected=list()
Fstat=matrix(0,nrow=split.time,ncol=kk*2)
set.seed(100)
for(j in 1:split.time){
  #split sample
  selected=list()
  unselected=1:N
  n=sample(unselected,size=N/kk,replace = FALSE)
  selected[[1]]=n
  selected[[2]]=unselected[-n]
  
  selectIV=list()
  for(p in 1:kk){
    n=selected[[p]]
    fit=screening(gene[-n,],X[-n],method = "sis")
    selectIV[[p]]=gene[,fit$screen]
  }
  
  for(p in 1:kk){
    gene_LD=selectIV[[p]]
    i=1
    while(i<ncol(gene_LD)){
      snps_to_remove=apply(as.matrix(gene_LD[,(i+1):ncol(gene_LD)]),2,function(x)filter_ld(x,gene_LD[,i]))
      if((i+1)==ncol(gene_LD)){
        if(snps_to_remove==1)gene_LD=gene_LD[,1:i]
      }else{
        gene_LD[,(i+1):ncol(gene_LD)]=gene_LD[,(i+1):ncol(gene_LD)][,!snps_to_remove]
      }
      i=i+1
    }
    selectIV[[p]]=gene_LD
  }
  #LASSO
  weights=list()
  for(p in 1:kk){
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
  length(selectIV[[1]][1,])
  #28
  length(selectIV[[2]][1,])
  #32
  
  fit1=lm(X[selected[[1]]]~.,data=selectIV[[1]][selected[[1]],])
  fit2=lm(X[selected[[2]]]~.,data=selectIV[[2]][selected[[2]],])
  Fstat[j,1]=summary(fit1)$fstat[1]
  Fstat[j,2]=summary(fit2)$fstat[1]
  
  #partial R2
  partialR2=list()
  for(p in 1:kk){
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
  hist(partialR2[[1]])
  hist(partialR2[[2]])
  
  selectmajor=list()
  selectweak=list()
  for(p in 1:kk){
    selectmajor[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]>=0.1])
    selectweak[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]<0.1])
    genemajorselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]>=0.1])
    geneweakselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]<0.1])
  }
  
  #Get combineMajor
  
  combineMajor=list()
  for(p in 1:kk){
    n=selected[[p]]
    weakwei=weights[[p]][partialR2[[p]]<0.1]
    weakwei=weightcal(weakwei)
    combineMajor[[p]]=cbind(selectmajor[[p]][n,],as.matrix(selectweak[[p]][n,])%*%weakwei)
  }
  #MAJOR
  residual=NULL
  for(p in 1:kk){
    n=selected[[p]]
    fit=lm(X[n]~.,data=as.data.frame(combineMajor[[p]]))
    Fstat[j,p+2]=summary(fit)$fstat[1]
    residual=c(residual,fit$residuals)
  }
  ord=c(selected[[1]],selected[[2]])
  fit2=lm(Y[ord]~residual+X[ord]+cova[ord])
  point.est1[j,]=summary(fit2)$coef[3,1:4]
  
  if(FALSE){ #CFMR
    combineCFMR=list()
    for(p in 1:kk){
      n=selected[[p]]
      combineCFMR[[p]]=as.matrix(selectIV[[p]][n,])%*%weightcal(weights[[p]])
    }
    combineCFMR=c(combineCFMR[[1]],combineCFMR[[2]])
    
    first1=lm(X[ord]~combineCFMR)
    fit2=lm(Y[ord]~first1$residuals+X[ord]+cova[ord])
    point.est2[j,]=summary(fit2)$coef[3,1:4]
  }
}

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
#write.table(alltime,file = "count time.txt")
#cauchy combination test
T1=mean(tan((0.5-point.est1[,4])*pi))
T2=mean(tan((0.5-point.est2[,4])*pi))
p1=0.5-atan(T1)/pi
p2=0.5-atan(T2)/pi
print(p1)
#0.001183231
print(p2)
#0.000655879
write.table(point.est1,file="CY_Major.txt")
write.table(point.est2,file="CY_CFMR.txt")

par(mfrow=c(2,2))
for(j in 1:4){
  hist(Fstat[,j],main=NULL,xlab = paste0("Fstatistics",j))
}


pvalue=data.frame(Major=point.est1[,4],
                  CFMR=point.est2[,4])
boxplot(pvalue,ylab="p-values")
abline(h=0.05,col="red")
beta=data.frame(Major=point.est1[,1],
                CFMR=point.est2[,1])
boxplot(beta,ylab=TeX("\\hat{\\beta}"))
par(mfrow=c(2,2),mar=c(4,4,1,2))
hist(point.est1[,4],breaks = 20,ylim=c(0,60),main = NULL,xlab = "p-values of 60 split times")
abline(v = 0.05, col = "red", lty = 2)
legend("topright", legend = "p-value=0.05", lty = 2, col = "red", cex = 0.8, bty = "n")
hist(pvalue[,2],breaks = 20,ylim=c(0,60),xlab="pvalues of CFMR",main=NULL)
abline(v = 0.05, col = "red", lty = 2)
legend("topright", legend = "p-value=0.05", lty = 2, col = "red", cex = 0.8, bty = "n")
hist(point.est1[,1], breaks=10,main = NULL, xlab = "Estimates of 60 split times")
abline(v = mean(point.est1[,1]), col = "red", lty = 2)
legend("topright", legend = expression(hat(beta) == 0.2778), lty = 2, col = "red", cex = 0.8, bty = "n")
hist(beta[,2], breaks = 20, xlab = "Estimates of CFMR", main = NULL)
abline(v = mean(beta[,2]), col = "red", lty = 2)
legend("topright", legend = expression(hat(beta) == 0.7932), lty = 2, col = "red", cex = 0.8, bty = "n")

#AOX1
X=AOX1$AOX1_exposure
gene<- read.table("AO_cleaned_gene_data.txt",header = TRUE)
N=nrow(gene)
I=ncol(gene)
#169,2886
cova=Y$AGE_.YRS.
Y=Y$CYP2E1_activity
set.seed(1)
gene <- as.data.frame(apply(gene, 2, fill_missing_values))

kk=2
split.time=60
r2_threshold <- 0.8
point.est1=matrix(0,nrow=split.time,ncol=4)
#point.est2=matrix(0,nrow=split.time,ncol=4)
colnames(point.est1)=c("Estimate","Std.error","t value","P.value")
#colnames(point.est2)=c("Estimate","Std.error","t value","P.value")
genemajorselected=list()
geneweakselected=list()
Fstat=matrix(0,nrow=split.time,ncol=kk*2)
set.seed(100)
for(j in 1:split.time){
  
  #split sample
  selected=list()
  unselected=1:N
  n=sample(unselected,size=N/kk,replace = FALSE)
  selected[[1]]=n
  selected[[2]]=unselected[-n]
  
  selectIV=list()
  for(p in 1:kk){
    n=selected[[p]]
    fit=screening(gene[-n,],X[-n],method = "sis")
    selectIV[[p]]=gene[,fit$screen]
  }
  
  
  for(p in 1:kk){
    gene_LD=selectIV[[p]]
    i=1
    while(i<ncol(gene_LD)){
      snps_to_remove=apply(as.matrix(gene_LD[,(i+1):ncol(gene_LD)]),2,function(x)filter_ld(x,gene_LD[,i]))
      if((i+1)==ncol(gene_LD)){
        if(snps_to_remove==1)gene_LD=gene_LD[,1:i]
      }else{
        gene_LD[,(i+1):ncol(gene_LD)]=gene_LD[,(i+1):ncol(gene_LD)][,!snps_to_remove]
      }
      i=i+1
    }
    selectIV[[p]]=gene_LD
  }
  #LASSO
  weights=list()
  for(p in 1:kk){
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
  length(selectIV[[1]][1,])
  #28
  length(selectIV[[2]][1,])
  #32
  
  fit1=lm(X[selected[[1]]]~.,data=selectIV[[1]][selected[[1]],])
  fit2=lm(X[selected[[2]]]~.,data=selectIV[[2]][selected[[2]],])
  Fstat[j,1]=summary(fit1)$fstat[1]
  Fstat[j,2]=summary(fit2)$fstat[1]
  
  #partial R2
  partialR2=list()
  for(p in 1:kk){
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
  hist(partialR2[[1]])
  hist(partialR2[[2]])
  
  selectmajor=list()
  selectweak=list()
  for(p in 1:kk){
    selectmajor[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]>=0.1])
    selectweak[[p]]=as.matrix(selectIV[[p]][,partialR2[[p]]<0.1])
    genemajorselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]>=0.1])
    geneweakselected[[(j-1)*2+p]]=c(colnames(selectIV[[p]])[partialR2[[p]]<0.1])
  }
  
  #Get combineMajor
  
  combineMajor=list()
  for(p in 1:kk){
    n=selected[[p]]
    weakwei=weights[[p]][partialR2[[p]]<0.1]
    weakwei=weightcal(weakwei)
    combineMajor[[p]]=cbind(selectmajor[[p]][n,],as.matrix(selectweak[[p]][n,])%*%weakwei)
  }
  #MAJOR
  residual=NULL
  for(p in 1:kk){
    n=selected[[p]]
    fit=lm(X[n]~.,data=as.data.frame(combineMajor[[p]]))
    Fstat[j,p+2]=summary(fit)$fstat[1]
    residual=c(residual,fit$residuals)
  }
  ord=c(selected[[1]],selected[[2]])
  fit2=lm(Y[ord]~residual+X[ord]+cova[ord])
  point.est1[j,]=summary(fit2)$coef[3,1:4]
  
  #CFMR
  if(FALSE){combineCFMR=list()
  for(p in 1:kk){
    n=selected[[p]]
    combineCFMR[[p]]=as.matrix(selectIV[[p]][n,])%*%weightcal(weights[[p]])
  }
  combineCFMR=c(combineCFMR[[1]],combineCFMR[[2]])
  
  first1=lm(X[ord]~combineCFMR)
  fit2=lm(Y[ord]~first1$residuals+X[ord]+cova[ord])
  point.est2[j,]=summary(fit2)$coef[3,1:4]
  }
}

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
write.table(alltime,file = "AO_count time.txt")
#cauchy combination test
T1=mean(tan((0.5-point.est1[,4])*pi))
T2=mean(tan((0.5-point.est2[,4])*pi))
p1=0.5-atan(T1)/pi
p2=0.5-atan(T2)/pi
print(p1)
#0.066
print(p2)
#0.060
write.table(point.est1,file="AO_Major.txt")
write.table(point.est2,file="AO_CFMR.txt")

pvalue=data.frame(Major=point.est1[,4],
                  CFMR=point.est2[,4])
boxplot(pvalue,ylab="p-values")
abline(h=0.05,col="red")
beta=data.frame(Major=point.est1[,1],
                CFMR=point.est2[,1])
boxplot(beta,ylab=TeX("\\hat{\\beta}"))
par(mfrow=c(2,2),mar=c(4,4,1,2))
hist(point.est1[,4],breaks = 20,ylim=c(0,10),main = NULL,xlab = "p-values")
abline(v = 0.05, col = "red", lty = 2)
legend("topright", legend = "p-value=0.05", lty = 2, col = "red", cex = 0.8, bty = "n")
hist(pvalue[,2],breaks = 20,ylim=c(0,60),xlab="pvalues of CFMR",main=NULL)
abline(v = 0.05, col = "red", lty = 2)
legend("topright", legend = "p-value=0.05", lty = 2, col = "red", cex = 0.8, bty = "n")
hist(point.est1[,1], breaks = 20, main = NULL, xlab = "Estimates")
abline(v = mean(point.est1[,1]), col = "red", lty = 2)
legend("topright", legend = expression(hat(beta) == 0.2778), lty = 2, col = "red", cex = 0.8, bty = "n")
hist(beta[,2], breaks = 20, xlab = "Estimates of CFMR", main = NULL)
abline(v = mean(beta[,2]), col = "red", lty = 2)
legend("topright", legend = expression(hat(beta) == 0.2643), lty = 2, col = "red", cex = 0.8, bty = "n")
