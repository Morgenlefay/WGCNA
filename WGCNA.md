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
##### Step3-Run WGCNA
```r
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = F, 
  verbose = 3
)
table(net$colors) 
```

##### Step4-Plot modules
```r
mergedColors = labels2colors(net$colors)
table(mergedColors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

##### Step5-Module-Trait correlation
```r
table(datTraits$subtype)
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  design=model.matrix(~0+ datTraits$subtype)
  colnames(design)=levels(datTraits$subtype)
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  tiff("step5-Module-trait-relationships.tiff",width = 800,height = 1200,res = 150)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  }
```


















