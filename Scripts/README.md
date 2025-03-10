# R Analysis
Philip Peffer
2025-01-25

- [<span class="toc-section-number">1</span> Exploratory
  analysis](#exploratory-analysis)
- [<span class="toc-section-number">2</span> Differential gene
  expression analysis](#differential-gene-expression-analysis)
- [<span class="toc-section-number">3</span> GSEA overalp
  analysis](#gsea-overalp-analysis)
- [<span class="toc-section-number">4</span> Gene Co-expression (WGCNA)
  Analysis](#gene-co-expression-wgcna-analysis)
- [<span class="toc-section-number">5</span> ssGSEA and survival
  analysis](#ssgsea-and-survival-analysis)
- [<span class="toc-section-number">6</span> Random forest prediction
  model](#random-forest-prediction-model)
- [<span class="toc-section-number">7</span> Model validation on
  independent dataset](#model-validation-on-independent-dataset)

In this project, I analyzed a gene expression dataset from tumor tissues
of muscle-invasive bladder cancer (MIBC) patients that received
platinum-based neoadjuvant chemotherapy (NAC) prior to radical
cystectomy (RC). The goal was to identify biomarkers that differentiate
patients that achieve a pathologic response to NAC from those that are
NAC resistant.

The analysis was performed using R within a Quarto markdown file
(Final_markdown_Peffer.qmd). Quarto is Posit’s enhanced version of R
Markdown that offers multi-language and multi-engine support. The
original rendered HTML output from the Quarto markdown can be viewed
[here](https://html-preview.github.io/?url=https://github.com/PhilipPeffer/MIBC_transcriptomic_biomarkers/blob/main/Scripts/Final_markdown_Peffer.html).
This README was created from a version of the Quarto file modified to
output to GitHub Flavored Markdown. The original code from the Quarto
markdown file has also been copied to a standalone R script
(Analysis_script.R) for easier viewing on GitHub.

# Exploratory analysis

``` r
# read in normalized expression data
# data was RMA-normalized, adjusted for labeling batch using ComBat24 and merged by median probe set values per gene symbol.
norm_data <- read.table("../Data/Discovery/GSE169455_normalized_by_gene.txt", header=T, sep='\t', row.names=1)
norm_data[1:10,1:10]
```

                 BLCA_Cx_73 BLCA_Cx_101 BLCA_Cx_118 BLCA_Cx_123 BLCA_Cx_130
    LOC100287497   5.637191    4.510258    4.726225    4.947787    4.745383
    SAMD11         4.791535    5.025635    5.748496    4.895198    4.924544
    KLHL17         4.573935    5.027887    5.248277    4.733811    4.924045
    PLEKHN1        4.881261    5.453513    5.570384    4.791831    5.310653
    ISG15          4.947811    4.886619    5.463361    5.698818    5.611286
    AGRN           5.377204    5.703746    6.023374    5.508075    5.915948
    MIR200B        4.484680    4.700068    4.842152    4.401712    5.262108
    MIR429         4.181502    4.361442    4.452980    4.275360    4.345431
    TTLL10         5.102119    5.340365    5.804617    5.301467    5.441642
    B3GALT6        4.796563    5.235698    5.555382    5.262584    5.269063
                 BLCA_Cx_131 BLCA_Cx_134 BLCA_Cx_143 BLCA_Cx_156 BLCA_Cx_178
    LOC100287497    5.864897    4.934421    4.689926    4.253382    4.725751
    SAMD11          4.376405    4.775536    5.258095    5.094424    5.136319
    KLHL17          4.293647    4.870542    4.910858    4.957031    4.820698
    PLEKHN1         4.637073    5.320277    5.107343    5.521946    5.395375
    ISG15           5.326671    4.546132    5.717496    7.392147    5.747426
    AGRN            5.448073    5.697307    5.518359    5.413662    5.941170
    MIR200B         4.725618    4.616541    4.731856    4.581271    5.109073
    MIR429          4.466166    4.806912    4.206000    4.457725    4.747474
    TTLL10          5.105598    5.322659    5.907872    5.664893    5.860805
    B3GALT6         5.131162    5.308575    5.342370    5.328338    5.350261

``` r
# read in patient phenotype data
p_data <- read.table("../Data/Discovery/GSE169455_pdata.txt", header=T, sep='\t', row.names=NULL)
p_data[1:10,1:10]
```

       GEO.Sample.ID         CEL.file Age.group..years.    Sex      Cohort
    1     BLCA_Cx_73     p1404_73.CEL             65-70   Male neoadjuvant
    2    BLCA_Cx_101    p1404_101.CEL             60-65   Male neoadjuvant
    3    BLCA_Cx_118    p1404_118.CEL             45-50 Female neoadjuvant
    4    BLCA_Cx_123 p1404_123_re.CEL             70-75   Male neoadjuvant
    5    BLCA_Cx_130    p1404_130.CEL             65-70   Male neoadjuvant
    6    BLCA_Cx_131 p1404_131_re.CEL             60-65   Male neoadjuvant
    7    BLCA_Cx_134 p1404_134_re.CEL             55-60   Male   induktion
    8    BLCA_Cx_143 p1404_143_re.CEL             65-70   Male neoadjuvant
    9    BLCA_Cx_156    p1404_156.CEL             65-70   Male   induktion
    10   BLCA_Cx_178    p1404_178.CEL             65-70   Male neoadjuvant
       cTNM.chemostart TURB.path.grade..WHO99. TURB.path.LVI
    1          cT3N0M0                      G3             1
    2         cT3bN0M0                      G2             0
    3          cT3N0M0                      G3           N/A
    4          cT2N0M0                      G3             0
    5          cT2N0M0                     N/A             1
    6          cT2N0M0                      G3             1
    7         cT4bN0M0                      G3             1
    8         cT3bN0M0                      G3             0
    9          cT3N2M0                      G3             0
    10         cT2N0M0                      G2             0
       TURB.path.keratinization TURB.path.cis
    1                         0             1
    2                         1             0
    3                         0             0
    4                         0             0
    5                       N/A             1
    6                         0             1
    7                         1             0
    8                         1             0
    9                         1             0
    10                        0             0

``` r
# remove induction chemo cohort (only interested in neoadjuvant chemo)
induc <- p_data$Cohort=="induktion"
neo_pdata <- p_data[!induc,]
neo_eset <- norm_data[,!induc]
dim(neo_pdata)[1]
```

    [1] 125

``` r
par(mfrow=c(1,1))
boxplot(neo_eset)
```

![](README_files/figure-commonmark/1.4-1.png)

``` r
library(ggplot2)

PCA <- prcomp(t(neo_eset), scale = FALSE)

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    Subtype = as.factor(neo_pdata$Consensus.subtype),
                    Response = as.factor(neo_pdata$Pathologic.response))

ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(shape = Subtype, colour = Response)) +
  ggtitle("PCA plot of the expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))
```

![](README_files/figure-commonmark/1.5-1.png)

``` r
# Heatmap of top 1000 genes with highest variance
library(ComplexHeatmap)

v <- apply(neo_eset, 1, var)
select=order(v,decreasing=T)[1:1000]

scl.expr <- t(scale(t(neo_eset[select,]), center=T, scale=T))

set.seed(50)
ha <- HeatmapAnnotation(sex=neo_pdata$Sex,
                        age=neo_pdata$Age.group..years.,
                        grade=neo_pdata$cTNM.chemostart,
                        survival=as.numeric(neo_pdata$CSS.time..months.),
                        chemo_courses=neo_pdata$Number.of.courses,
                        regimen=neo_pdata$Regimen.type,
                        subtype=neo_pdata$Consensus.subtype,
                        response=neo_pdata$Pathologic.response)

# row (gene) and column (sample) hierarchical clustering
library(dendsort)
corr.dist=function(x) { as.dist(1-cor(t(x))) }
row_dend = dendsort(hclust(corr.dist(scl.expr), method="ward.D2"))
col_dend = dendsort(hclust(corr.dist(t(scl.expr)), method="ward.D2"))

# generate heatmap
par(mfrow=c(1,1))
Heatmap(scl.expr,
        cluster_rows = row_dend,
        cluster_columns = col_dend,
        row_names_gp = gpar(fontsize = 1),
        show_column_names = FALSE,
        top_annotation = ha,
        heatmap_legend_param = list(title="expr"),
        show_row_names = FALSE)
```

![](README_files/figure-commonmark/1.6-1.png)

# Differential gene expression analysis

``` r
# DEGs associated with Complete Response (CR) using SAM

#install.packages("samr")
library(samr)

# group partial responders with complete responders
completeResponse <- function(x) {
  if(x=='no pR'){
    return(2)
  }
  else{
    return(1)
  }
}

neo_pdata$CR2 <- as.numeric(sapply(neo_pdata$Pathologic.response, completeResponse))


# Perform differential gene expression analysis using SAM
samfit <- SAM(x=as.matrix(neo_eset),
              y=neo_pdata$CR2, resp.type="Two class unpaired",
              genenames=rownames(neo_eset),
              nperms=100,
              testStatistic="standard",
              regression.method="standard",
              fdr.output = 0.05,
              random.seed=123)
```

    perm= 1
    perm= 2
    perm= 3
    perm= 4
    perm= 5
    perm= 6
    perm= 7
    perm= 8
    perm= 9
    perm= 10
    perm= 11
    perm= 12
    perm= 13
    perm= 14
    perm= 15
    perm= 16
    perm= 17
    perm= 18
    perm= 19
    perm= 20
    perm= 21
    perm= 22
    perm= 23
    perm= 24
    perm= 25
    perm= 26
    perm= 27
    perm= 28
    perm= 29
    perm= 30
    perm= 31
    perm= 32
    perm= 33
    perm= 34
    perm= 35
    perm= 36
    perm= 37
    perm= 38
    perm= 39
    perm= 40
    perm= 41
    perm= 42
    perm= 43
    perm= 44
    perm= 45
    perm= 46
    perm= 47
    perm= 48
    perm= 49
    perm= 50
    perm= 51
    perm= 52
    perm= 53
    perm= 54
    perm= 55
    perm= 56
    perm= 57
    perm= 58
    perm= 59
    perm= 60
    perm= 61
    perm= 62
    perm= 63
    perm= 64
    perm= 65
    perm= 66
    perm= 67
    perm= 68
    perm= 69
    perm= 70
    perm= 71
    perm= 72
    perm= 73
    perm= 74
    perm= 75
    perm= 76
    perm= 77
    perm= 78
    perm= 79
    perm= 80
    perm= 81
    perm= 82
    perm= 83
    perm= 84
    perm= 85
    perm= 86
    perm= 87
    perm= 88
    perm= 89
    perm= 90
    perm= 91
    perm= 92
    perm= 93
    perm= 94
    perm= 95
    perm= 96
    perm= 97
    perm= 98
    perm= 99
    perm= 100

    Computing delta table
    1
    2
    3
    4
    5
    6
    7
    8
    9
    10
    11
    12
    13
    14
    15
    16
    17
    18
    19
    20
    21
    22
    23
    24
    25
    26
    27
    28
    29
    30
    31
    32
    33
    34
    35
    36
    37
    38
    39
    40
    41
    42
    43
    44
    45
    46
    47
    48
    49
    50

``` r
print(samfit)
```

    Call:
    SAM(x = as.matrix(neo_eset), y = neo_pdata$CR2, resp.type = "Two class unpaired", 
        genenames = rownames(neo_eset), nperms = 100, testStatistic = "standard", 
        regression.method = "standard", random.seed = 123, fdr.output = 0.05)

    Genes up
          Gene ID  Gene Name Score(d) Numerator(r) Denominator(s+s0) Fold Change
     [1,] HOXD10   7663      3.611    0.75         0.208             1.167      
     [2,] APP      8732      2.944    0.378        0.128             1.059      
     [3,] PARP6    4456      2.885    0.292        0.101             1.054      
     [4,] SPP1     10010     2.872    0.739        0.257             1.131      
     [5,] MGC34796 318       2.596    0.183        0.071             1.057      
     [6,] LCE3A    1164      2.545    0.185        0.073             1.041      
     [7,] GRP      6173      2.514    0.364        0.145             1.091      
     [8,] ZNF532   6171      2.449    0.242        0.099             1.051      
     [9,] MAP1B    10373     2.443    0.397        0.163             1.095      
    [10,] RTL8B    13497     2.429    0.173        0.071             1.037      
          q-value(%)
     [1,] 0         
     [2,] 0         
     [3,] 0         
     [4,] 0         
     [5,] 0         
     [6,] 0         
     [7,] 0         
     [8,] 20.07     
     [9,] 20.07     
    [10,] 20.07     

    Genes down
    NULL

``` r
plot(samfit)
```

![](README_files/figure-commonmark/2.1-1.png)

``` r
# Extract list of DEGs associated with CR (FDR < 0.05)
deg.tab <- samfit$siggenes.table

up.genes <- data.frame(deg.tab$genes.up)
up.genes$q.value... <- as.numeric(up.genes$q.value...)
up.genes$Fold.Change <- as.numeric(up.genes$Fold.Change)
down.genes <- data.frame(deg.tab$genes.lo)
down.genes$q.value... <- as.numeric(down.genes$q.value...)
down.genes$Fold.Change <- -as.numeric(down.genes$Fold.Change)

deg.tab <- rbind(up.genes, down.genes)
deg.tab$q.value... <- deg.tab$q.value.../100

sig.up <- up.genes[up.genes$q.value...<5,]
sig.down <- down.genes[down.genes$q.value...<5,]
sig.genes <- deg.tab[deg.tab$q.value...<0.05,]

gene.list <- sig.genes$Gene.ID
print(gene.list)
```

    [1] "HOXD10"   "APP"      "PARP6"    "SPP1"     "MGC34796" "LCE3A"    "GRP"     

``` r
# write significant gene list to file for gene set overlap analysis
write.table(gene.list, file="../Results/CR_DEG_list.txt", quote=F, sep="", row.names=F, col.names=F)
```

``` r
# Heatmap of DEGs w/ annotations
library(ComplexHeatmap)

scl.expr <- t(scale(t(neo_eset[gene.list,]), center=T, scale=T))

# Heatmap column annotations
set.seed(100)
CRsort_pdata <- neo_pdata
CRsort_pdata$CR <- ifelse(CRsort_pdata$CR == 1, "responder", "non-responder")
CRsort_pdata <- CRsort_pdata[order(CRsort_pdata$CR2),]
ha <- HeatmapAnnotation(sex=CRsort_pdata$Sex,
                        age=CRsort_pdata$Age.group..years.,
                        grade=CRsort_pdata$cTNM.chemostart,
                        survival=as.numeric(CRsort_pdata$CSS.time..months.),
                        chemo_courses=CRsort_pdata$Number.of.courses,
                        regimen=CRsort_pdata$Regimen.type,
                        subtype=CRsort_pdata$Consensus.subtype,
                        response=CRsort_pdata$CR)

# Hierarchical clustering
library(dendsort)
corr.dist=function(x) { as.dist(1-cor(t(x))) }
row_dend = dendsort(hclust(corr.dist(scl.expr), method="ward.D2"))

# Heatmap
par(mfrow=c(1,1))
Heatmap(scl.expr[,CRsort_pdata$GEO.Sample.ID],
        cluster_rows = row_dend,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE,
        top_annotation = ha,
        heatmap_legend_param = list(title="expr"))
```

![](README_files/figure-commonmark/2.3-1.png)

# GSEA overalp analysis

GSEA overlap results (FDR \< 0.05):  
![GSEA overlap results - significant gene sets
table](https://github.com/PhilipPeffer/MIBC_transcriptomic_biomarkers/blob/main/Results/CR_GS_overlap1.png?raw=true)  
![GSEA overlap results - gene/geneset overlap
matrix](https://github.com/PhilipPeffer/MIBC_transcriptomic_biomarkers/blob/main/Results/CR_GS_overlap2.png?raw=true)  

# Gene Co-expression (WGCNA) Analysis

``` r
# install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel","matrixStats") ) 
# BiocManager::install("impute")
# install.packages("WGCNA")
library(matrixStats)
library(WGCNA)
allowWGCNAThreads() 
```

    Allowing multi-threading with up to 24 threads.

``` r
# Select most variable genes
v <- apply(neo_eset, 1, var)
select=order(v,decreasing=T)[1:1000]

dataExpr2<-t(neo_eset[select,])
#we need to transpose our matrix for the format WGCNA expects 


#Create Network these defaults will work for many situations, but you can try 
#changing the minimum module size if you get many, small modules. You can also 
#change the split level to anything between 1 and 4. Merged cut height determines 
#where in the dendrogram you cut a branch into a module

net2 <- blockwiseModules(dataExpr2, power = 7,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold =10, mergeCutHeight = 0.10, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
                       verbose=3, ds=3)
```

     Calculating module eigengenes block-wise from all genes
       Flagging genes and samples with too many missing values...
        ..step 1
     ..Working on block 1 .
        TOM calculation: adjacency..
        ..will not use multithreading.
         Fraction of slow calculations: 0.000000
        ..connectivity..
        ..matrix multiplication (system BLAS)..
        ..normalization..
        ..done.
     ....clustering..
     ....detecting modules..
     ....calculating module eigengenes..
     ....checking kME in modules..
         ..removing 5 genes from module 3 because their KME is too low.
         ..removing 3 genes from module 4 because their KME is too low.
      ..reassigning 20 genes from module 1 to modules with higher KME.
      ..reassigning 4 genes from module 2 to modules with higher KME.
      ..reassigning 6 genes from module 3 to modules with higher KME.
      ..reassigning 6 genes from module 4 to modules with higher KME.
      ..reassigning 4 genes from module 5 to modules with higher KME.
     ..merging modules that are too close..
         mergeCloseModules: Merging modules whose distance is less than 0.1
           Calculating new MEs...

``` r
#now we get to see our modules!
mergedColors2 <- labels2colors(net2$colors) 
plotDendroAndColors(net2$dendrograms[[1]], mergedColors2[net2$blockGenes[[1]]], "Module colors",
                                                             dendroLabels = FALSE, hang = 0.03,
                                                             addGuide = TRUE, guideHang = 0.05)
```

![](README_files/figure-commonmark/3.1-1.png)

``` r
#save modules
genes2<-colnames(dataExpr2) 
moduleColors2<-labels2colors(net2$colors) 
mymodules2<-cbind(genes2,moduleColors2) 
moduleColors2<-labels2colors(net2$colors) 
mymodules2<-data.frame(cbind(genes2,moduleColors2,paste("geneModule_",net2$colors)))


#get Eigengenes
MEs02 <- moduleEigengenes(dataExpr2, moduleColors2)$eigengenes
MEs2 <- orderMEs(MEs02)

#get pathological response status
phenotype_num2 <- neo_pdata$CR2

#correlate geneModule with response
moduleTraitCor2 <- cor(MEs2, phenotype_num2, use = "p") 
moduleTraitCor2
```

                       [,1]
    MEbrown     -0.05738071
    MEyellow     0.11071491
    MEgreen     -0.02629087
    MEblue      -0.10901833
    MEturquoise  0.03647357
    MEgrey       0.04114259

``` r
corPvalueStudent(moduleTraitCor2, nrow(dataExpr2))
```

                     [,1]
    MEbrown     0.5250275
    MEyellow    0.2190053
    MEgreen     0.7710222
    MEblue      0.2261915
    MEturquoise 0.6863418
    MEgrey      0.6487065

``` r
#Select the geneModule that has highest correlation
module_high_cor2 <- moduleTraitCor2[order(abs(moduleTraitCor2),decreasing = T),]
module_high_cor2
```

       MEyellow      MEblue     MEbrown      MEgrey MEturquoise     MEgreen 
     0.11071491 -0.10901833 -0.05738071  0.04114259  0.03647357 -0.02629087 

``` r
module_high_cor_name2 <- sub('..', '', names(module_high_cor2[1]))

biomarker <- mymodules2[mymodules2$moduleColors2 %in% module_high_cor_name2,"genes2"]
biomarker
```

      [1] "S100A2"   "KRT6A"    "FN1"      "ANXA1"    "POSTN"    "LUM"     
      [7] "COL1A1"   "S100A10"  "LAMC2"    "SFRP4"    "ITGA2"    "DCN"     
     [13] "THBS1"    "SERPINE1" "THBS2"    "CD44"     "KRT17"    "KRT5"    
     [19] "CHI3L1"   "KRT6B"    "COL1A2"   "COL6A3"   "DSP"      "TIMP1"   
     [25] "BGN"      "ACTG2"    "VCAN"     "ASPN"     "TMEM45A"  "S100A3"  
     [31] "SFN"      "TIMP2"    "CHST11"   "DCBLD2"   "LAMB3"    "SPARC"   
     [37] "ACTA2"    "SFRP2"    "LGALS1"   "COL3A1"   "INHBA"    "LDHA"    
     [43] "CCDC80"   "SULF1"    "MMP2"     "MSN"      "PLOD2"    "COL12A1" 
     [49] "GREM1"    "CTSK"     "VIM"      "CAV1"     "PLAT"     "CTSL"    
     [55] "EMP1"     "AEBP1"    "SULF2"    "ITGA3"    "AHNAK2"   "ANXA2"   
     [61] "EGR1"     "TNC"      "MYADM"    "MMP1"     "CAV2"     "LY6E"    
     [67] "ITGB4"    "PRNP"     "CASZ1"    "MT1X"     "SNAI2"    "ACTN1"   
     [73] "PMEPA1"   "EFEMP1"   "TGFBI"    "C1R"      "OSMR"     "PLK2"    
     [79] "MYL9"     "FLNA"     "GLIPR1"   "SOCS3"    "FBN1"     "F13A1"   
     [85] "SERPINE2" "CD248"    "IFITM3"   "RND3"     "RAB31"    "MXRA5"   
     [91] "CRYAB"    "FGFR1"    "CD99"     "RGS2"     "IGF2BP2"  "MMP14"   
     [97] "PLXDC2"   "TMSB10"   "ANTXR1"   "HEG1"     "PMP22"    "LPAR1"   
    [103] "BICC1"    "CNN2"     "OLFML2B"  "FKBP10"   "TREM1"    "CALD1"   
    [109] "SLCO2A1"  "NT5E"     "APBB2"    "ITGA5"    "PHLDA1"   "CDH11"   
    [115] "EMP3"     "SCARNA22" "DPT"      "CNN1"     "LOXL2"    "RNF217"  

``` r
# Determine if there is overlap between the DEGs and the module
DE_overlap2 <- gene.list %in% biomarker
print(sum(DE_overlap2))
```

    [1] 0

``` r
print(gene.list[DE_overlap2])
```

    character(0)

``` r
# write the genes from the highest correlated modules to a text file
write.table(biomarker, "../Results/biomarker.txt", row.names=F, col.names=F, quote=F)
```

GSEA overlap results (FDR \< 0.05):  
![GSEA overlap results - significant gene sets
table](https://github.com/PhilipPeffer/MIBC_transcriptomic_biomarkers/blob/main/Results/CR_WGCNA_overlap1.png?raw=true)  
![GSEA overlap results - gene/geneset overlap
matrix](https://github.com/PhilipPeffer/MIBC_transcriptomic_biomarkers/blob/main/Results/CR_WGCNA_overlap2.png?raw=true)  

# ssGSEA and survival analysis

``` r
##insall pakcages if needed
# BiocManager::install("singscore")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("GSEABase")

library(singscore)
library(SummarizedExperiment)
library(GSEABase)

#create a summarized and ranked object for singscore analysis
se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(neo_eset)))
ranked <- rankGenes(se0)

#create a gene set signature from biomarker genes (DEGs)
geneset_Sig <- GeneSet(as.character(gene.list), geneIdType=SymbolIdentifier())

#perform the singsocre single sample GSEA analysis
scoredf <- simpleScore(ranked, upSet = geneset_Sig, knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)

#order all samples based on the singscore
final_score <- scoredf[order(scoredf$TotalScore,decreasing = T),]

#get a number corrspond to 30% of patients 
no.patient <- nrow(final_score) * 0.35

# 35% patients that have the highest singscore
highBioMarkerSigGroup <- head(final_score,no.patient)
write.csv(highBioMarkerSigGroup,"../Results/highBioMarkerSigGroup.csv")

# 35% patients that have the lowest singscore
lowBioMarkerSigGroup <- tail(final_score,no.patient)
write.csv(lowBioMarkerSigGroup,"../Results/lowBioMarkerSigGroup.csv")
```

``` r
sig_group <- c()
for(i in 1:dim(final_score)[1]){
    if(i <= dim(final_score)[1]/2){
        sig_group[i] <- "high"
    }
    else{
        sig_group[i] <- "low"
    }
}

final_score$sig_group <- sig_group

ggplot(final_score, aes(x = dim(final_score)[1]:1, y = TotalScore, color = sig_group)) +
    geom_point() +
    labs(x = "Patients (increasing biomarker score)", y = "Biomarker activity score")
```

![](README_files/figure-commonmark/4.2-1.png)

``` r
# Heatmap of CSS DEGs w/ gene signature score annotations
sig_eset <- neo_eset[gene.list, rownames(final_score)]

scl.expr <- t(scale(t(sig_eset), center=T, scale=T))

set.seed(123)
ha <- HeatmapAnnotation(Biomarker_Score=final_score$TotalScore)

# Row hierarchical clustering
library(dendsort)
corr.dist=function(x) { as.dist(1-cor(t(x))) }
row_dend = dendsort(hclust(corr.dist(scl.expr), method="ward.D2"))

# Heatmap
par(mfrow=c(1,1))
Heatmap(scl.expr,
        cluster_rows = row_dend,
        cluster_columns = F,
        row_names_gp = gpar(fontsize = 8),
        show_column_names = FALSE,
        top_annotation = ha,
        heatmap_legend_param = list(title="expr"))
```

![](README_files/figure-commonmark/4.2-2.png)

``` r
# survival analysis
library(survival)

high_pdata <- neo_pdata[(neo_pdata$GEO.Sample.ID %in% rownames(highBioMarkerSigGroup)),]
high_pdata$Sig <- "high"
low_pdata <- neo_pdata[(neo_pdata$GEO.Sample.ID %in% rownames(lowBioMarkerSigGroup)),]
low_pdata$Sig <- "low"
surv_pdata <- rbind(high_pdata, low_pdata)
surv_eset <- neo_eset[,surv_pdata$GEO.Sample.ID]

t2d <- Surv(time=as.numeric(surv_pdata$OS.time..months.), event=as.numeric(surv_pdata$OS.event), type="right")

x <- coxph(t2d~Sig, data=surv_pdata)
summary(x)
```

    Call:
    coxph(formula = t2d ~ Sig, data = surv_pdata)

      n= 86, number of events= 39 
       (1 observation deleted due to missingness)

              coef exp(coef) se(coef)      z Pr(>|z|)    
    Siglow -1.1870    0.3051   0.3490 -3.401 0.000671 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

           exp(coef) exp(-coef) lower .95 upper .95
    Siglow    0.3051      3.277     0.154    0.6047

    Concordance= 0.657  (se = 0.035 )
    Likelihood ratio test= 12.76  on 1 df,   p=4e-04
    Wald test            = 11.57  on 1 df,   p=7e-04
    Score (logrank) test = 12.94  on 1 df,   p=3e-04

``` r
fit <- survfit(t2d~Sig, data=surv_pdata)

library(survminer)
ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = surv_pdata,             # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   xlab = "Time in months",   # customize X axis label.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
   risk.table.y.text.col = T, # colour risk table text annotations.
   risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
)
```

![](README_files/figure-commonmark/4.3-1.png)

# Random forest prediction model

``` r
# Make a table containing the independent variables (gene signature expression values) and dependent variable (response) for each patient
training_data <- as.data.frame(t(neo_eset[gene.list,]))
training_data$Response <- neo_pdata$CR2-1

library(randomForestSRC)

# Fit a Random Forest model
set.seed(123)
rf_model <- rfsrc(Response ~ ., 
                  data = training_data, 
                  ntree = 500,  # Number of trees
                  importance = TRUE)  # Enables variable importance calculation

# View model summary
print(rf_model)
```

                             Sample size: 125
                         Number of trees: 500
               Forest terminal node size: 5
           Average no. of terminal nodes: 12.488
    No. of variables tried at each split: 3
                  Total no. of variables: 7
           Resampling used to grow trees: swor
        Resample size used to grow trees: 79
                                Analysis: RF-R
                                  Family: regr
                          Splitting rule: mse *random*
           Number of random split points: 10
                         (OOB) R squared: 0.23955609
       (OOB) Requested performance error: 0.18809948

``` r
# View feature importance
importance <- rf_model$importance
importance <- importance[order(importance, decreasing=T)]
print(importance)
```

         PARP6   MGC34796      LCE3A     HOXD10       SPP1        APP        GRP 
    0.07744624 0.04702895 0.03421514 0.03245276 0.01769435 0.01694057 0.01250145 

``` r
# Save the model
save(rf_model, file = "../Results/rf_model.RData")

# Make predictions
pred <- predict(rf_model, newdata = training_data)
class_probs <- pred$predicted # predicted class probabilites

library(pROC)
library(caret)

# Compute the ROC curve and AUC
roc_curve <- roc(training_data$Response, class_probs)
auc_value <- auc(roc_curve)

# Plot the ROC curve
par(mfrow=c(1,1))
plot(roc_curve, main = paste("AUC =", round(auc_value, 3)))
```

![](README_files/figure-commonmark/5.1-1.png)

``` r
# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))
```

    [1] "Optimal Threshold = 0.576675396825397"

``` r
# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)
```

      threshold sensitivity specificity
    1 0.5766754   0.9577465   0.9814815

``` r
# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(training_data$Response))
print(conf_matrix)
```

    Confusion Matrix and Statistics

              Reference
    Prediction  0  1
             0 53  3
             1  1 68
                                              
                   Accuracy : 0.968           
                     95% CI : (0.9201, 0.9912)
        No Information Rate : 0.568           
        P-Value [Acc > NIR] : <2e-16          
                                              
                      Kappa : 0.9351          
                                              
     Mcnemar's Test P-Value : 0.6171          
                                              
                Sensitivity : 0.9815          
                Specificity : 0.9577          
             Pos Pred Value : 0.9464          
             Neg Pred Value : 0.9855          
                 Prevalence : 0.4320          
             Detection Rate : 0.4240          
       Detection Prevalence : 0.4480          
          Balanced Accuracy : 0.9696          
                                              
           'Positive' Class : 0               
                                              

``` r
# 5-fold cross-validation
set.seed(123)
folds <- createFolds(training_data$Response, k = 5)
auc_values <- c()
sens_values <- c()
spec_values <- c()
accuracy <- c()

par(mfrow=c(3,2))
for (i in seq_along(folds)) {
  # Split into training and test data
  train_idx <- setdiff(seq_len(nrow(training_data)), folds[[i]])
  train_data <- training_data[train_idx, ]
  test_data <- training_data[folds[[i]], ]
  
  # Train Random Forest model
  rf_model <- rfsrc(Response ~ ., data = train_data, ntree = 500)
  
  # Predict risk scores on test set
  class_probs <- predict(rf_model, newdata = test_data)$predicted
  
  # Compute AUC
  roc_curve <- roc(test_data$Response, class_probs)
  auc_values <- c(auc_values, auc(roc_curve))

  # Plot ROC curve
  plot(roc_curve, main = paste("Fold", i, "ROC Curve"))

  # Identify the best threshold
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

  # Extract sensitivity and specificity at the optimal threshold
  optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  sens_values <- c(sens_values, optimal_coords$sensitivity)
  spec_values <- c(spec_values, optimal_coords$specificity)

  # Predicted classes based on optimal threshold
  predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)
  acc <- sum(ifelse(predicted_classes == test_data$Response, 1, 0))/length(predicted_classes)
  accuracy <- c(accuracy, acc)
}

# Mean accuracy, sensitivity, specificity, and AUC across folds
mean_acc <- mean(accuracy)
mean_sens <- mean(sens_values)
mean_spec <- mean(spec_values)
mean_auc <- mean(auc_values)
paste("Mean accuracy =", mean_acc)
```

    [1] "Mean accuracy = 0.752"

``` r
paste("Mean sensitivity =", mean_sens)
```

    [1] "Mean sensitivity = 0.824615384615385"

``` r
paste("Mean specificity =", mean_spec)
```

    [1] "Mean specificity = 0.677863247863248"

``` r
paste("Mean AUC =", mean_auc)
```

    [1] "Mean AUC = 0.757910256410256"

![](README_files/figure-commonmark/5.1-2.png)

# Model validation on independent dataset

``` r
# load validation set expression data
val_eset <- read.csv("../Data/Validation/41467_2020_18640_MOESM7_ESM(Supplementary Data 4).csv", row.names=1)

# load NAC patient list
nac_pats <- read.table("../Data/Validation/validation_NAC_patients.txt", header=F)
nac_pats <- nac_pats$V1
val_eset <- val_eset[, colnames(val_eset) %in% nac_pats]

# load validation set phenotype data
val_pdata <- read.table("../Data/Validation/validation_pdata2.txt", header=T, sep="\t", row.names=1)
val_pdata <- val_pdata[colnames(val_eset),]

# Biomarker gene IDs not present in the validation dataset
gene.list[!(gene.list %in% rownames(val_eset))]
```

    [1] "MGC34796"

``` r
# Biomarker after removing gene IDs not present in the validation data
final.list <- gene.list[gene.list %in% rownames(val_eset)]

# Make a table containing the independent variables (gene signature expression values) and dependent variables (survival info) for each patient
val_data <- as.data.frame(t(val_eset[final.list,]))
val_data$Response <- as.numeric(factor(val_pdata$Response, levels=c("No response", "Response")))-1
head(val_data)
```

         HOXD10  APP PARP6 SPP1 LCE3A   GRP Response
    P043  -2.57 9.14  4.08 7.15 -2.57 -2.57        1
    P074   1.93 8.58  3.48 4.49 -2.57 -0.80        1
    P045   1.68 8.65  4.28 1.68 -2.57 -2.57        0
    P009   4.73 7.85  2.05 5.06 -2.57 -2.57        0
    P090   0.14 8.79  3.71 5.45 -2.57 -2.57        1
    P092   3.30 8.98  3.72 3.22 -2.57 -2.57        1

``` r
# Compare validation expression distributions to training data
par(mfrow=c(2,1))
boxplot(t(training_data[, final.list]))
boxplot(t(val_data[, final.list]))
```

![](README_files/figure-commonmark/6.1-1.png)

``` r
# Z-score transformation of validation data (RNA-seq)
val_data[, final.list] <- scale(val_data[, final.list])
head(val_data)
```

               HOXD10          APP       PARP6        SPP1 LCE3A        GRP
    P043 -2.137699492  0.561514811  0.59949429  1.19443561   NaN -0.3776944
    P074 -0.006764872  0.001428791 -0.05677227 -0.30998769   NaN  1.2086221
    P045 -0.125150129  0.071439543  0.81824981 -1.89924688   NaN -0.3776944
    P009  1.319150003 -0.728683342 -1.62087423  0.01238873   NaN -0.3776944
    P090 -0.854403310  0.211461048  0.19479658  0.23296207   NaN -0.3776944
    P092  0.641986335  0.401490234  0.20573436 -1.02826497   NaN -0.3776944
         Response
    P043        1
    P074        1
    P045        0
    P009        0
    P090        1
    P092        1

``` r
# Remove gene that could not be scaled from the biomarker
final.list <- final.list[-5]
val_data <- val_data[,c(final.list,"Response")]

# Print final biomarker genes
print(final.list)
```

    [1] "HOXD10" "APP"    "PARP6"  "SPP1"   "GRP"   

``` r
# Validation set does not contain all of the genes in the gene list
training_data <- as.data.frame(t(neo_eset[final.list,]))
training_data$Response <- neo_pdata$CR2-1

# Z-score transformation of training data (microarray)
training_data[, final.list] <- scale(training_data[, final.list])
par(mfrow=c(2,1))
boxplot(t(training_data[, final.list]))
boxplot(t(val_data[, final.list]))
```

![](README_files/figure-commonmark/6.2-1.png)

``` r
library(randomForestSRC)

# Fit a Random Forest model
set.seed(123)
rf_model <- rfsrc(Response ~ ., 
                  data = training_data, 
                  ntree = 500,  # Number of trees
                  importance = TRUE)  # Enables variable importance calculation

# View model summary
print(rf_model)
```

                             Sample size: 125
                         Number of trees: 500
               Forest terminal node size: 5
           Average no. of terminal nodes: 13.856
    No. of variables tried at each split: 2
                  Total no. of variables: 5
           Resampling used to grow trees: swor
        Resample size used to grow trees: 79
                                Analysis: RF-R
                                  Family: regr
                          Splitting rule: mse *random*
           Number of random split points: 10
                         (OOB) R squared: 0.14458667
       (OOB) Requested performance error: 0.21159063

``` r
# View feature importance
importance <- rf_model$importance
importance <- importance[order(importance, decreasing=T)]
print(importance)
```

         PARP6     HOXD10        APP       SPP1        GRP 
    0.07870175 0.03293762 0.02070231 0.01978755 0.01908221 

``` r
# Save the model
save(rf_model, file = "../Results/rf_model2.RData")

# Make predictions
pred <- predict(rf_model, newdata = training_data)
class_probs <- pred$predicted # predicted class probabilites

library(pROC)
library(caret)

# Compute the ROC curve and AUC
roc_curve <- roc(training_data$Response, class_probs)
auc_value <- auc(roc_curve)

# Plot the ROC curve
par(mfrow=c(1,1))
plot(roc_curve, main = paste("AUC =", round(auc_value, 3)))
```

![](README_files/figure-commonmark/6.2-2.png)

``` r
# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))
```

    [1] "Optimal Threshold = 0.534823015873016"

``` r
# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)
```

      threshold sensitivity specificity
    1  0.534823   0.9295775   0.9444444

``` r
# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(training_data$Response))
print(conf_matrix)
```

    Confusion Matrix and Statistics

              Reference
    Prediction  0  1
             0 51  5
             1  3 66
                                             
                   Accuracy : 0.936          
                     95% CI : (0.8778, 0.972)
        No Information Rate : 0.568          
        P-Value [Acc > NIR] : <2e-16         
                                             
                      Kappa : 0.8702         
                                             
     Mcnemar's Test P-Value : 0.7237         
                                             
                Sensitivity : 0.9444         
                Specificity : 0.9296         
             Pos Pred Value : 0.9107         
             Neg Pred Value : 0.9565         
                 Prevalence : 0.4320         
             Detection Rate : 0.4080         
       Detection Prevalence : 0.4480         
          Balanced Accuracy : 0.9370         
                                             
           'Positive' Class : 0              
                                             

``` r
# 5-fold cross-validation
set.seed(123)
folds <- createFolds(training_data$Response, k = 5)
auc_values <- c()
sens_values <- c()
spec_values <- c()
accuracy <- c()

par(mfrow=c(3,2))
for (i in seq_along(folds)) {
  # Split into training and test data
  train_idx <- setdiff(seq_len(nrow(training_data)), folds[[i]])
  train_data <- training_data[train_idx, ]
  test_data <- training_data[folds[[i]], ]
  
  # Train Random Forest model
  rf_model <- rfsrc(Response ~ ., data = train_data, ntree = 500)
  
  # Predict risk scores on test set
  class_probs <- predict(rf_model, newdata = test_data)$predicted
  
  # Compute AUC
  roc_curve <- roc(test_data$Response, class_probs)
  auc_values <- c(auc_values, auc(roc_curve))

  # Plot ROC curve
  plot(roc_curve, main = paste("Fold", i, "ROC Curve"))

  # Identify the best threshold
  optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

  # Extract sensitivity and specificity at the optimal threshold
  optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")
  sens_values <- c(sens_values, optimal_coords$sensitivity)
  spec_values <- c(spec_values, optimal_coords$specificity)

  # Predicted classes based on optimal threshold
  predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)
  acc <- sum(ifelse(predicted_classes == test_data$Response, 1, 0))/length(predicted_classes)
  accuracy <- c(accuracy, acc)
}

# Mean accuracy, sensitivity, specificity, and AUC across folds
mean_acc <- mean(accuracy)
mean_sens <- mean(sens_values)
mean_spec <- mean(spec_values)
mean_auc <- mean(auc_values)
paste("Mean accuracy =", mean_acc)
```

    [1] "Mean accuracy = 0.68"

``` r
paste("Mean sensitivity =", mean_sens)
```

    [1] "Mean sensitivity = 0.672564102564103"

``` r
paste("Mean specificity =", mean_spec)
```

    [1] "Mean specificity = 0.716581196581197"

``` r
paste("Mean AUC =", mean_auc)
```

    [1] "Mean AUC = 0.667752136752137"

![](README_files/figure-commonmark/6.2-3.png)

``` r
# Assess predictive performance on the validation dataset
# Load the model
load("../Results/rf_model2.RData")

# Make predictions
pred <- predict(rf_model, newdata = val_data)
class_probs <- pred$predicted # predicted class probabilites

library(pROC)
library(caret)

# Compute the ROC curve and AUC
roc_curve <- roc(val_data$Response, class_probs)
auc_value <- auc(roc_curve)

# Plot the ROC curve
par(mfrow=c(1,1))
plot(roc_curve, main = paste("AUC =", round(auc_value, 3)))
```

![](README_files/figure-commonmark/6.3-1.png)

``` r
# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))
```

    [1] "Optimal Threshold = 0.592842857142857"

``` r
# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)
```

      threshold sensitivity specificity
    1 0.5928429   0.7692308         0.5

``` r
# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(val_data$Response))
print(conf_matrix)
```

    Confusion Matrix and Statistics

              Reference
    Prediction  0  1
             0  4  3
             1  4 10
                                              
                   Accuracy : 0.6667          
                     95% CI : (0.4303, 0.8541)
        No Information Rate : 0.619           
        P-Value [Acc > NIR] : 0.4183          
                                              
                      Kappa : 0.2759          
                                              
     Mcnemar's Test P-Value : 1.0000          
                                              
                Sensitivity : 0.5000          
                Specificity : 0.7692          
             Pos Pred Value : 0.5714          
             Neg Pred Value : 0.7143          
                 Prevalence : 0.3810          
             Detection Rate : 0.1905          
       Detection Prevalence : 0.3333          
          Balanced Accuracy : 0.6346          
                                              
           'Positive' Class : 0               
                                              

``` r
##insall pakcages if needed
# BiocManager::install("singscore")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("GSEABase")

library(singscore)
library(SummarizedExperiment)
library(GSEABase)

#create a summarized and ranked object for singscore analysis
se0 <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(val_eset)))
ranked <- rankGenes(se0)

#create a gene set signature from biomarker genes
geneset_Sig <- GeneSet(as.character(final.list), geneIdType=SymbolIdentifier())

#perform the singsocre single sample GSEA analysis
scoredf <- simpleScore(ranked, upSet = geneset_Sig, knownDirection = TRUE, centerScore = TRUE, subSamples = NULL)

#order all samples based on the singscore
final_score <- scoredf[order(scoredf$TotalScore,decreasing = T),]

#get a number corrspond to 30% of patients 
no.patient <- nrow(final_score) * 0.35

# 30% patients that have the highest singscore
highBioMarkerSigGroup <- head(final_score,no.patient)

# 30% patients that have the lowest singscore
lowBioMarkerSigGroup <- tail(final_score,no.patient)


# survival analysis
library(survival)

high_pdata <- val_pdata[(rownames(val_pdata) %in% rownames(highBioMarkerSigGroup)),]
high_pdata$Sig <- "high"
low_pdata <- val_pdata[(rownames(val_pdata) %in% rownames(lowBioMarkerSigGroup)),]
low_pdata$Sig <- "low"
surv_pdata <- rbind(high_pdata, low_pdata)
surv_eset <- val_eset[,rownames(surv_pdata)]

t2d <- Surv(time=as.numeric(surv_pdata$OS), event=as.numeric(surv_pdata$Deceased), type="right")

x <- coxph(t2d~Sig, data=surv_pdata)
summary(x)
```

    Call:
    coxph(formula = t2d ~ Sig, data = surv_pdata)

      n= 15, number of events= 5 

               coef exp(coef) se(coef)     z Pr(>|z|)
    Siglow -0.01263   0.98745  1.22935 -0.01    0.992

           exp(coef) exp(-coef) lower .95 upper .95
    Siglow    0.9875      1.013   0.08873     10.99

    Concordance= 0.583  (se = 0.154 )
    Likelihood ratio test= 0  on 1 df,   p=1
    Wald test            = 0  on 1 df,   p=1
    Score (logrank) test = 0  on 1 df,   p=1

``` r
fit <- survfit(t2d~Sig, data=surv_pdata)

library(survminer)
ggsurvplot(
   fit,                     # survfit object with calculated statistics.
   data = surv_pdata,       # data used to fit survival curves.
   risk.table = TRUE,       # show risk table.
   pval = TRUE,             # show p-value of log-rank test.
   xlab = "Time in months", # customize X axis label.
   ggtheme = theme_light(), # customize plot and risk table with a theme.
 risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
                            # in legend of risk table
)
```

![](README_files/figure-commonmark/6.4-1.png)
