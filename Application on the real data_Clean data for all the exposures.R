setwd("C:/Users/shirx/OneDrive/Desktop/Major/Application on real data")
#outcome
data=read.table("liverdata/curatedPhenotype/phenotype.txt",header = FALSE,row.names = 1)
transposed_data <- t(data)
transposed_data <- as.data.frame(transposed_data)
rownames(transposed_data) <- colnames(data)
colnames(transposed_data) <- rownames(data)
Y=transposed_data
Y$individual_id=as.character(Y$individual_id)
Y$CYP2E1_activity=as.numeric(Y$CYP2E1_activity)
na_rows=which(is.na(Y$CYP2E1_activity))
Y=Y[-na_rows,]
Y=Y[!(Y$CYP2E1_activity==0),]
Y$CYP2E1_activity=log(Y$CYP2E1_activity)

equallist=makeequal(Y$individual_id,finalID)
Y=Y[equallist[[1]],]


#gene
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
equallist=makeequal(Y$individual_id,geneID)
gene=gene[equallist[[2]],]
Y=Y[equallist[[1]],]

#exposure
data=read.table("liverdata/curatedExpressionLiver/expression.txt",header = FALSE,row.names = 1)
transposed_data <- t(data)
transposed_data <- as.data.frame(transposed_data)
rownames(transposed_data) <- colnames(data)
colnames(transposed_data) <- rownames(data)
data=transposed_data
feature=read.table("liverdata/curatedExpressionLiver/features.txt",header = TRUE,fill = TRUE)

#target: CYP2E1
ex=grep("CYP2E1", feature$genesymbol)
featureid=feature[ex,1]
X=matrix(0,nrow = length(data[,1]),ncol = length(ex))
for(i in 1:length(ex)){
  j=grep(featureid[i],colnames(data))
  X[,i]=data[,j]
  X[,i]=10^(X[,i])
}
X=rowMeans(X)
CYP2E1=data.frame("ID"=data$feature_id,
             "CYP2E1_exposure"=X)
na_rows=which(is.na(CYP2E1$CYP2E1_exposure))
CYP2E1=CYP2E1[-na_rows,]
equallist=makeequal(Y$individual_id,CYP2E1$ID)
gene=gene[equallist[[1]],]
Y=Y[equallist[[1]],]
CYP2E1=CYP2E1[equallist[[2]],]
finalID=CYP2E1$ID
write.table(finalID,"finalID.txt",row.names = FALSE,col.names = FALSE)
Y=Y[,c("individual_id","AGE_(YRS)","CYP2E1_activity")]
Y$`AGE_(YRS)`=as.numeric(Y$`AGE_(YRS)`)
write.table(Y,"finaloutcome.txt",row.names = FALSE)

fit=lm(Y$CYP2E1_activity~Y$`AGE_(YRS)`+CYP2E1$CYP2E1_exposure)
summary(fit)

#AOX1
ex=grep("AOX1", feature$genesymbol)
featureid=feature[ex,1]
X=matrix(0,nrow = length(data[,1]),ncol = length(ex))
for(i in 1:length(ex)){
  j=grep(featureid[i],colnames(data))
  X[,i]=data[,j]
  X[,i]=10^(X[,i])
}
X=rowMeans(X)
AOX1=data.frame("ID"=data$feature_id,
                "AOX1_exposure"=X)
na_rows=which(is.na(AOX1$AOX1_exposure))
if(length(na_rows)>0)AOX1=AOX1[-na_rows,]
Y=read.table("finaloutcome.txt",header = TRUE)
equallist=makeequal(Y$individual_id,AOX1$ID)
AOX1=AOX1[equallist[[2]],]

fit=lm(Y$CYP2E1_activity~Y$`AGE_(YRS)`+AOX1$AOX1_exposure)
summary(fit)
