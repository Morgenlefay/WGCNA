##### Load data
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
