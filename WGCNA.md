##### Step1-Load data
```r
### R packages
rm(list=ls())
options(stringsAsFactors = F)
library(WGCNA)

### Load data
RNAexpr <- read.csv("Gene_expression.csv",sep=",",header = T)
fpkm <- RNAexpr
rownames(fpkm)=fpkm[,1]
fpkm=fpkm[,-1]

### Sample info
subname=sapply(colnames(fpkm),function(x) strsplit(x,"_")[[1]][1])
datTraits = data.frame(gsm=names(fpkm),
                       subtype=subname)
rownames(datTraits)=datTraits[,1]

### Proprocess
RNAseq_voom <- fpkm 
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr <- WGCNA_matrix  ## top 5000 mad genes
```

##### Step2-Determine soft power
```r
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
best_beta=sft$powerEstimate
best_beta
```









