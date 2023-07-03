library(HardyWeinberg)
setwd("C:/Users/shirx/OneDrive/Desktop/Major/Application on real data")
data=read.table("liverdata/curatedGenotype/genotype.txt",row.names = 1,header = TRUE)
transposed_data <- t(data)
transposed_data <- as.data.frame(transposed_data)
rownames(transposed_data) <- colnames(data)
colnames(transposed_data) <- rownames(data)
gene=transposed_data
geneID=gsub("X", "", rownames(gene))
makeequal <- function(A, B) {
  a=NULL
  b=NULL
  i <- 1
  q=1
  while (i <= length(A)) {
    B_id=grep(A[i],B)
    if (length(B_id)>0) {
      a[q]=i
      b[q]=B_id
      q=q+1
    }
    i=i+1
  }
  return(list(a,b))
}
Y=read.table("finaloutcome.txt",header = TRUE)
equallist=makeequal(Y$individual_id,geneID)
gene=gene[equallist[[2]],]

#delete X
features=read.table("liverdata/curatedGenotype/features.txt",header = TRUE,fill = TRUE)
genename=c(features$feature_id[features$chrom=="X"],features$feature_id[features$chrom=="A"])
position=match(genename,colnames(gene))
position=position[!is.na(position)]
gene=gene[,-position]
#441386

if(FALSE){
  Data=read.table("exposure and outcome.txt",header = TRUE)
  ID=Data[,1]#ID
  geneID=gsub("X", "", rownames(gene))
  ID=as.character(ID)
  makeequal <- function(A, B) {
    a=NULL
    b=NULL
    i <- 1
    q=1
    while (i <= length(A)) {
      B_id=grep(A[i],B)
      if (length(B_id)>0) {
        a[q]=i
        b[q]=B_id
        q=q+1
      }
      i=i+1
    }
    return(list(a,b))
  }
  equallist=makeequal(ID,geneID)
  data=Data[equallist[[1]],]
  X=data[,3]#exposure
  Y=data[,2]#outcome
  gene=gene[equallist[[2]],]
  #169, 449699
  data$self_reported_ethnicity=as.factor(data$self_reported_ethnicity)
  summary(data$self_reported_ethnicity)
  #population
  data$inferred_population=as.factor(data$inferred_population)
  summary(data$inferred_population)
  #all equal
  data=data[,-8]
  data=data[,-4]
  #gender
  data$GENDER=as.factor(data$GENDER)
  cova=data[,c(1,4:10)]
  pval=NULL
  for(i in 2:ncol(cova)){
    fit=lm(Y~cova[,i])
    pval[i-1]=summary(fit)$coef[2,4]
  }
  #only age significant
  cova=cova[,1:2]
  data=data.frame(cova,X,Y)
  write.table(data,"exposure and outcome_after gene.txt",row.names = FALSE)
  
}

#clean gene
is_same <- function(str) {
  substr(str, 1, 1) == substr(str, 2, 2)
}
get_larger_pos <- function(x) {
  if (x[1] > x[2]) {
    return(1)
  } else {
    return(2)
  }
}
convert_data <- function(x) {
  counts <- table(x)
  site <- is_same(names(counts))
  new <- rep(NA,length(x))
  if (length(counts) == 3) {
    new[x == names(counts)[!site]] <- 1
    larger <- get_larger_pos(counts[site])
    new[x == names(counts)[site][larger]] <- 0
    new[x == names(counts)[site][-larger]] <- 2
  }
  if (length(counts) == 2) {
    if (sum(site) == 1) {
      new[x == names(counts)[!site]] <- 1
      new[x == names(counts)[site]] <- 0
    }
    if (sum(site) == 0) {
      larger <- get_larger_pos(counts[site])
      new[x == names(counts)[site][larger]] <- 0
      new[x == names(counts)[site][-larger]] <- 2
    }
  }
  if (length(counts) == 1) {
    if (sum(site) == 1) {
      new[x == names(counts)[site]] <- 0
    }
    if (sum(site) == 0) {
      new[x == names(counts)[!site]] <- 1
    }
  }
  return(new)
}
clean_data <- apply(gene, 2, convert_data)
convert_data <- data.frame(clean_data)

#delete na>0.1, and unique value
delete=function(x){
  if (sum(is.na(x)) > length(x) * 0.1) {
    return(1)
  } else if (length(unique(na.omit(x))) == 1) {
    return(1) 
  }
  else{
    return(0)
  }
}
index.delete <- apply(convert_data, 2, delete)
data=convert_data[,!index.delete]
#169 obs. of 361322

calc_maf <- function(column) {
  freq <- table(column) / length(column)
  maf=ifelse(sum(names(freq)=="2")==1,freq[names(freq)=="2"],0)+ifelse(sum(names(freq)=="0")==1,0.5*freq[names(freq)=="0"],0)
  return(maf)
}
calc_hwp <- function(column) {
  obs <- table(column)
  if (length(obs) == 3) {
    n=sum(obs)
    freq=obs/n
    freq_exp=c((freq[1]+freq[2]/2)^2,2*(freq[1]+freq[2]/2)*(freq[3]+freq[2]/2),(freq[3]+freq[2]/2)^2)
    result=1-pchisq(sum((freq*n-freq_exp*n)^2/freq_exp/n), df=2)
    return(result)
  } else {
    n <- length(column)
    k <- sum(na.omit(column) == 1)
    p_value <- binom.test(k, n, p = 0.5)$p.value
    return(p_value)
  }
}
maf_values <- apply(data, 2, calc_maf)
keep_columns <- which(maf_values >= 0.05)
genotype_data <- data[, keep_columns]
#361321

#select gene
X=CYP2E1$CYP2E1_exposure
getpval=function(column){
  fit=lm(X~column)
  pval=summary(fit)$coef[2,4]
  return(pval)
}
pval <- apply(genotype_data, 2, getpval)
keep_columns <- which(pval<0.01)
length(keep_columns)
gene<- genotype_data[, keep_columns]

hwp_values <- apply(gene, 2, calc_hwp)
keep_columns <- which(hwp_values >= 1e-5)
gene <- gene[, keep_columns]
write.table(gene,file="CY_cleaned_gene_data.txt",row.names = FALSE)

X=AOX1$AOX1_exposure
getpval=function(column){
  fit=lm(X~column)
  pval=summary(fit)$coef[2,4]
  return(pval)
}
pval <- apply(genotype_data, 2, getpval)
keep_columns <- which(pval<0.01)
length(keep_columns)
gene<- genotype_data[, keep_columns]

hwp_values <- apply(gene, 2, calc_hwp)
keep_columns <- which(hwp_values >= 1e-5)
gene <- gene[, keep_columns]
write.table(gene,file="AO_cleaned_gene_data.txt",row.names = FALSE)
