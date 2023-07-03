library(ivreg)
library(MASS)
library(ggplot2)
library(hdi)
library(glmnet)
library(dplyr)
library(qqman)
library(screening)
library(ivmodel)
#Y outcome, D exposure
#Z instruments, X covariates
J=300#number of SNPs
J1=5 #useful SNPs
N=1000 #sample size
MAF=0.3 #minor allele frequency
Times=1000
sig=c(0.1,0.2)
squareh=c(0.15,0.3,0.5)
#pi=pi*c(1,-1,1,-1,1)
beta=c(-0.08,0,0.08) #causal effect
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
xun=1
for(ss in 1:2){
  for(hh in 1:3){
    allresult=matrix(0,nrow=6,ncol=10)
    rownames(allresult)=rep(c("Major", "TSLS", "LIML"),2)
    colnames(allresult)=c("squareh","sigma","Real.beta","mean","bias","sd","RMSE","cp","power","TyI")
    
    qq=1
    for(jj in 1:3){
      #TSLS
      point.est1=matrix(NA,ncol=2,nrow=Times)
      #LIML
      point.est3=matrix(NA,ncol=2,nrow=Times)
      #MAJOR
      point.est2=matrix(NA,ncol=2,nrow=Times)
      for(t in 1:Times){
        #generate SNPs
        G=matrix(0,nrow=N,ncol=J)
        for(i in 1:J){
          maf=MAF[i%%length(MAF)]
          if(i%%length(MAF)==0)maf=MAF[length(MAF)]
          GG1=runif(N,min=0,max=1)
          GG1[GG1<=(1-maf)^2]=0
          GG1[GG1>(1-maf^2)]=2
          GG1[which(GG1!=0&GG1!=2)]=1
          G[,i]=GG1
        }
        Sigma=matrix(c(1,sig[ss],sig[ss],1),nrow=2,ncol=2)
        sumsquarepi=squareh[hh]*Sigma[2,2]/(1-squareh[hh])/(2*MAF*(1-MAF))
        pi=sqrt(sumsquarepi*c(0.4,0.4,0.1,0.05,0.05))
        error=mvrnorm(N,mu=c(0,0),Sigma=Sigma)
        #randomly choose SNPs that generates the esposure
        real=sample(c(1:J),size=J1,replace = FALSE)
        Greal=G[,real]
        #generate exposure
        X=as.matrix(Greal)%*%as.matrix(pi)+error[,1]
        #generate outcome
        Y=X*beta[jj]+error[,2]
        
        #split sample
        selected=NULL
        unselected=1:N
        n=sample(unselected,size=N/k,replace = FALSE)
        selected=cbind(selected,n)
        for(p in 1:k){
          if(p<k){
            n=sample(unselected[-selected],size=N/k,replace = FALSE)
            selected=cbind(selected,n)
          }
        }
        #select IVs using SIS
        selectIV=list()
        for(p in 1:k){
          n=selected[,p]
          if(J<N/2){
            fit=screening(G[-n,],X[-n],method = "sis",num.select = 100)
          }
          if(J>N/2){
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
          orderR=order(-partialR2[[p]])
          selectmajor[[p]]=selectIV[[p]][partialR2[[p]]>=0.1]
          selectweak[[p]]=selectIV[[p]][partialR2[[p]]<0.1]
        }
        #Get coef
        combineMajor=list()
        for(p in 1:k){
          orderF=order(-partialR2[[p]])
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
        point.est2[t,]=summary(fit2)$coef[3,1:2]
        
        #TSLS and LIML
        #select IVs using SIS
        selectIV=NULL
        if(J<N/2){
          fit=screening(G,X,method = "sis",num.select = 100)
        }
        if(J>N/2){
          fit=screening(G,X,method = "sis")
        }
        selectIV=fit$screen
        #LASSO
        #select IVs 
        lambda=cv.glmnet(G[,selectIV],X,alpha=1)$lambda.min
        lasso=glmnet(G[,selectIV],X,lambda=lambda)
        #coef of all IVs
        a=coef(lasso)[-1,]
        #selected IVs
        selectIV=selectIV[a!=0]
        SNP=G[,selectIV]
        fit1=ivmodel(Y=Y,D=X,Z=SNP)
        point.est1[t,]=c(fit1$kClass$point.est[2],fit1$kClass$std.err[2])
        point.est3[t,]=c(fit1$LIML$point.est,fit1$LIML$std.err)
      }
      
      result=data.frame("squareh"=rep(squareh[hh],3),
                        "sigma"=rep(Sigma[1,2],3),
                        "Real.beta"=rep(beta[[jj]],3),
                        "mean"=c(mean(point.est2[,1]),mean(point.est1[,1]),mean(point.est3[,1])),
                        "bias"=c(mean(abs(point.est2[,1]-beta[jj])),mean(abs(point.est1[,1]-beta[jj])),mean(abs(point.est3[,1]-beta[jj]))),
                        "sd"=c(mean(point.est2[,2]),mean(point.est1[,2]),mean(point.est3[,2])),
                        "RMSE"=c(sqrt(sum((point.est2[,1]-beta[[jj]])^2)/Times),sqrt(sum((point.est1[,1]-beta[[jj]])^2)/Times),sqrt(sum((point.est3[,1]-beta[[jj]])^2)/Times)))
      rownames(result)=c("MAJOR","TSLS","LIML")
      #95% confidence interval
      CI.l1=point.est1[,1]-qnorm(0.975)*point.est1[,2]
      CI.l2=point.est2[,1]-qnorm(0.975)*point.est2[,2]
      CI.l3=point.est3[,1]-qnorm(0.975)*point.est3[,2]
      CI.U1=point.est1[,1]+qnorm(0.975)*point.est1[,2]
      CI.U2=point.est2[,1]+qnorm(0.975)*point.est2[,2]
      CI.U3=point.est3[,1]+qnorm(0.975)*point.est3[,2]
      reject=rep(0,3)
      for(i in 1:Times){
        if(0<CI.l1[i]||0>CI.U1[i])reject[2]=reject[2]+1
        if(0<CI.l2[i]||0>CI.U2[i])reject[1]=reject[1]+1
        if(0<CI.l3[i]||0>CI.U3[i])reject[3]=reject[3]+1
      }
      reject=reject/Times
      result$reject=reject
      
      cp=rep(0,3)
      for(i in 1:Times){
        if(beta[[jj]]>CI.l1[i]&&beta[[jj]]<CI.U1[i])cp[2]=cp[2]+1
        if(beta[[jj]]>CI.l2[i]&&beta[[jj]]<CI.U2[i])cp[1]=cp[1]+1
        if(beta[[jj]]>CI.l3[i]&&beta[[jj]]<CI.U3[i])cp[3]=cp[3]+1
      }
      result$cp=cp/Times
      writefile=paste0("squareh_",hh,"_sigma",ss,"beta_",beta[jj],".txt")
      writedata=data.frame("MAJOR"=point.est2[,1],"LIML"=point.est3[,1],"TSLS"=point.est1[,1])
      write.table(writedata,file=writefile)#save estimates results
      
      if(beta[[jj]]!=0){
        allresult[qq:(qq+2),1:9]=as.matrix(cbind(result[,1:7],result$cp,result$reject))[,1:9]
        qq=qq+3
      }
      if(beta[[jj]]==0){
        allresult[,10]=rep(result$reject,2)
      }
    }
    sumresult[[xun]]=allresult
    print(sumresult[[xun]])
    xun=xun+1
  }
}
write.table(sumresult,file="sumresult.txt")

