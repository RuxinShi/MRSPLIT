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
J=c(300,3000,30000) #number of SNPs
J1=5 #useful SNPs
N=c(1000,3000,5000) #sample size
MAF=0.3 #minor allele frequency
Times=1000
sig=0.8
squareh=c(0.15,0.2,0.3)#extremely weak, weak, strong
ratio=matrix(0,nrow=5,ncol=2)
ratio[,1]=c(0.4,0.4,0.1,0.05,0.05)#weights of SNPs
ratio[,2]=c(0.2,0.2,0.2,0.2,0.2)
beta=c(-0.08,-0.05,0,0.05,0.08) #causal effect
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

sumresult=list()
snp=1
xun=1
for(hh in 1:3){
  for(rr in 1:2){
    for(ii in 1:3){
      allresult=matrix(0,nrow=8,ncol=9)
      rownames(allresult)=rep(c("Major", "CFMR"),4)
      colnames(allresult)=c("ratiotype","Real.beta","mean","bias","sd","RMSE","cp","power","TyI")
      qq=1
      for(jj in 1:5){
        #CFMR
        point.est1=matrix(NA,ncol=2,nrow=Times)
        #MAJOR
        point.est2=matrix(NA,ncol=2,nrow=Times)
        for(t in 1:Times){
          G=matrix(0,nrow=N[ii],ncol=J[snp])
          for(i in 1:J[snp]){
            maf=MAF[i%%length(MAF)]
            if(i%%length(MAF)==0)maf=MAF[length(MAF)]
            GG1=runif(N[ii],min=0,max=1)
            GG1[GG1<=(1-maf)^2]=0
            GG1[GG1>(1-maf^2)]=2
            GG1[which(GG1!=0&GG1!=2)]=1
            G[,i]=GG1
          }
          Sigma=matrix(c(5,sig,sig,5),nrow=2,ncol=2)
          sumsquarepi=squareh[hh]*Sigma[2,2]/(1-squareh[hh])/(2*MAF*(1-MAF))
          pi=sqrt(sumsquarepi*ratio[,rr])
          error=mvrnorm(N[ii],mu=c(0,0),Sigma=Sigma)
          real=sample(c(1:J[snp]),size=J1,replace = FALSE)
          Greal=G[,real]
          X=as.matrix(Greal)%*%as.matrix(pi)+error[,1]
          Y=X*beta[jj]+error[,2]
          
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
            if(J[snp]<N[ii]/2){
              fit=screening(G[-n,],X[-n],method = "sis",num.select = 100)
            }
            if(J[snp]>N[ii]/2){
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
          #get coef of CFMR
          combineCFMR=matrix(NA,ncol=k,nrow=(N[ii]-N[ii]/k))
          for(p in 1:k){
            n=selected[,p]
            if(length(selectIV[[p]])>=2){
              combineCFMR[,p]=G[n,selectIV[[p]]]%*%weightcal(weights[[p]])
            }
            if(length(selectIV[[p]])<2){
              combineCFMR[,p]=G[n,selectIV[[p]]]
            }
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
          #CFMR
          combineCFMR=c(combineCFMR)
          first1=lm(X[c(selected)]~combineCFMR)
          fit1=lm(Y[c(selected)]~first1$residuals+X[c(selected)])
          point.est1[t,]=summary(fit1)$coef[3,1:2]
          #MAJOR
          residual=NULL
          for(p in 1:k){
            n=selected[,p]
            fit=lm(X[n]~combineMajor[[p]])
            residual=c(residual,fit$residuals)
          }
          fit2=lm(Y[c(selected)]~residual+X[c(selected)])
          point.est2[t,]=summary(fit2)$coef[3,1:2]
        }
        
        result=data.frame("ratiotype"=rep(rr,2),
                          "Real.beta"=rep(beta[[jj]],2),
                          "mean"=c(mean(point.est2[,1]),mean(point.est1[,1])),
                          "bias"=c(mean(abs(point.est2[,1]-beta[jj])),mean(abs(point.est1[,1]-beta[jj]))),
                          "sd"=c(mean(point.est2[,2]),mean(point.est1[,2])),
                          "RMSE"=c(sqrt(sum((point.est2[,1]-beta[[jj]])^2)/Times),sqrt(sum((point.est1[,1]-beta[[jj]])^2)/Times)))
        rownames(result)=c("MAJOR","CFMR")
        #95% confidence interval
        CI.l1=point.est1[,1]-qnorm(0.975)*point.est1[,2]
        CI.l2=point.est2[,1]-qnorm(0.975)*point.est2[,2]
        CI.U1=point.est1[,1]+qnorm(0.975)*point.est1[,2]
        CI.U2=point.est2[,1]+qnorm(0.975)*point.est2[,2]
        reject=rep(0,2)
        for(i in 1:Times){
          if(0<CI.l1[i]||0>CI.U1[i])reject[2]=reject[2]+1
          if(0<CI.l2[i]||0>CI.U2[i])reject[1]=reject[1]+1
        }
        reject=reject/Times
        result$reject=reject
        
        cp=rep(0,2)
        for(i in 1:Times){
          if(beta[[jj]]>CI.l1[i]&&beta[[jj]]<CI.U1[i])cp[2]=cp[2]+1
          if(beta[[jj]]>CI.l2[i]&&beta[[jj]]<CI.U2[i])cp[1]=cp[1]+1
        }
        result$cp=cp/Times
        writefile=paste0("squareh",squareh[hh],"size_",ii,"ratio_",rr,"beta_",beta[jj],".txt")
        writedata=data.frame("MAJOR"=point.est2[,1],"CFMR"=point.est1[,1])
        write.table(writedata,file=writefile)#save estimates results
        
        if(beta[[jj]]!=0){
          allresult[qq:(qq+1),1:8]=as.matrix(cbind(result[,1:6],result$cp,result$reject))[,1:8]
          qq=qq+2
        }
        if(beta[[jj]]==0){
          allresult[,9]=rep(result$reject,4)
        }
      }
      sumresult[[xun]]=allresult
      print(sumresult[[xun]])
      xun=xun+1
    }   
  }
}
write.table(sumresult,file="sumresultCFMR.txt")
