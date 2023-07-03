#different sample size, different splitting time
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
N=c(1000,3000) #sample size
mul=c(150,100)
MAF=0.3 #minor allele frequency
Times=1000
rho=0.1
squareh=c(0.15,0.2,0.3)#extremely weak, weak, moderate
ratio=matrix(0,nrow=5,ncol=2)
ratio[,1]=c(0.4,0.4,0.1,0.05,0.05)
ratio[,2]=c(0.2,0.2,0.2,0.2,0.2)
beta=0.08 #causal effect
k=2
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
set.seed(1)
sumresult=list()
rr=1
for(ii in 1:length(N)){
  for(hh in 1:1){
    #MAJOR
    point.est=matrix(NA,ncol=2,nrow=mul[ii]*Times)
    for(t in 1:Times){
      G=matrix(0,nrow=N[ii],ncol=J)
      for(i in 1:J){
        maf=MAF[i%%length(MAF)]
        if(i%%length(MAF)==0)maf=MAF[length(MAF)]
        GG1=runif(N[ii],min=0,max=1)
        GG1[GG1<=(1-maf)^2]=0
        GG1[GG1>(1-maf^2)]=2
        GG1[which(GG1!=0&GG1!=2)]=1
        G[,i]=GG1
      }
      Sigma=matrix(c(1,rho,rho,1),nrow=2,ncol=2)
      sumsquarepi=squareh[hh]*Sigma[2,2]/(1-squareh[hh])/(2*MAF*(1-MAF))
      pi=sqrt(sumsquarepi*ratio[,rr])
      error=mvrnorm(N[ii],mu=c(0,0),Sigma=Sigma)
      real=sample(c(1:J),size=J1,replace = FALSE)
      Greal=G[,real]
      X=as.matrix(Greal)%*%as.matrix(pi)+error[,1]
      Y=X*beta+error[,2]
      for(tt in 1:mul[ii]){
        #split sample
        selected=NULL
        unselected=1:N[ii]
        n=sample(unselected,size=N[ii]/k,replace = FALSE)
        selected=cbind(selected,n)
        for(p in 1:k){
          if(p<k){
            n=sample(unselected[-selected],size=N[ii]/k,replace = FALSE)
            selected=cbind(selected,n)
          }
        }
        #select IVs using SIS
        selectIV=list()
        for(p in 1:k){
          n=selected[,p]
          if(J<N[ii]/2){
            fit=screening(G[-n,],X[-n],method = "sis",num.select = 100)
          }
          if(J>N[ii]/2){
            fit=screening(G[-n,],X[-n],method = "sis")
          }
          selectIV[[p]]=fit$screen
        }
        #LASSO
        weights=list()
        for(p in 1:k){
          n=selected[,p]
          #select IVs 
          lambda=cv.glmnet(G[-n,selectIV[[p]]],X[-n],alpha=1)$lambda.min
          lasso=glmnet(G[-n,selectIV[[p]]],X[-n],lambda=lambda)
          #coef of all IVs
          a=coef(lasso)[-1,]
          while(length(a[a!=0])==0){
            lambda=cv.glmnet(G[-n,selectIV[[p]]],X[-n],alpha=1)$lambda.min
            lasso=glmnet(G[-n,selectIV[[p]]],X[-n],lambda=lambda)
            #coef of all IVs
            a=coef(lasso)[-1,]
          }
          weights[[p]]=a[a!=0]
          #selected IVs
          selectIV[[p]]=selectIV[[p]][a!=0]
        }
        #select major and weak based on partial R2 statistics
        partialR2=list()
        for(p in 1:k){
          n=selected[,p]
          full=lm(X[-n]~G[-n,selectIV[[p]]])
          sse_full=sum((full$residuals)^2)
          R2=NULL
          if(length(selectIV[[p]])==1){
            R2=summary(lm(X[-n]~G[-n,selectIV[[p]]]))$r.squared
          }
          if(length(selectIV[[p]])>1){
            for(q in 1:length(selectIV[[p]])){
              reduced=lm(X[-n]~G[-n,selectIV[[p]][-q]])
              sse_reduced=sum((reduced$residuals)^2)
              R2[q]=(sse_reduced-sse_full)/sse_reduced
            }
          }
          partialR2[[p]]=R2
        }
        selectmajor=list()
        selectweak=list()
        for(p in 1:k){
          selectmajor[[p]]=selectIV[[p]][partialR2[[p]]>=0.1]
          selectweak[[p]]=selectIV[[p]][partialR2[[p]]<0.1]
        }
        
        #Get coef
        combineMajor=list()
        for(p in 1:k){
          n=selected[,p]
          weakwei=weights[[p]][partialR2[[p]]<0.1]
          if(length(weakwei)>1){
            weakwei=weightcal(weakwei)
            combineMajor[[p]]=cbind(G[n,selectmajor[[p]]],G[n,selectweak[[p]]]%*%weakwei)
          } 
          if(length(weakwei)<=1) combineMajor[[p]]=G[n,selectIV[[p]]]
        }
        #MAJOR
        residual=NULL
        for(p in 1:k){
          n=selected[,p]
          fit=lm(X[n]~combineMajor[[p]])
          residual=c(residual,fit$residuals)
        }
        fit2=lm(Y[c(selected)]~residual+X[c(selected)])
        point.est[tt+(t-1)*mul[ii],]=summary(fit2)$coef[3,c(1,4)]
      }
    }
    write.table(point.est,file=paste0("multiple_","h_",squareh[hh],"N_",N[ii],".txt"),row.names = FALSE,col.names = FALSE)
   
    }
}


