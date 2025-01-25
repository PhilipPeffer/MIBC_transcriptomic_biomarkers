# ---
# title: "RBIF114 Final Project"
# author: "Philip Peffer"


# In this project, I analyze a gene expression dataset from tumor tissues of muscle-invasive bladder cancer (MIBC) patients 
# that received platinum-based neoadjuvant chemotherapy (NAC) prior to radical cystectomy (RC). The goal was to identify biomarkers 
# that differentiate patients that achieve a pathologic response to NAC from those that are NAC resistant.





## Exploratory analysis


# read in normalized expression data
# data was RMA-normalized, adjusted for labeling batch using ComBat24 and merged by median probe set values per gene symbol.
norm_data <- read.table("../Data/Discovery/GSE169455_normalized_by_gene.txt", header=T, sep='\t', row.names=1)
norm_data[1:10,1:10]



# read in patient phenotype data
p_data <- read.table("../Data/Discovery/GSE169455_pdata.txt", header=T, sep='\t', row.names=NULL)
p_data[1:10,1:10]



# remove induction chemo cohort (only interested in neoadjuvant chemo)
induc <- p_data$Cohort=="induktion"
neo_pdata <- p_data[!induc,]
neo_eset <- norm_data[,!induc]
dim(neo_pdata)[1]



par(mfrow=c(1,1))
boxplot(neo_eset)



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





## Differential gene expression analysis


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
print(samfit)
plot(samfit)



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

# write significant gene list to file for gene set overlap analysis
write.table(gene.list, file="../Results/CR_DEG_list.txt", quote=F, sep="", row.names=F, col.names=F)



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




## GSEA overlap analysis


# utilized GSEA/MSigDB webtool: https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp

#GSEA overlap results (FDR < 0.05):
#GSEA overlap results - significant gene sets table: refer to 'CR_GS_overlap1.png' in results folder
#GSEA overlap results - gene/geneset overlap matrix: refer to 'CR_GS_overlap2.png' in results folder




## Gene Co-expression (WGCNA) Analysis


# install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel","matrixStats") ) 
# BiocManager::install("impute")
# install.packages("WGCNA")
library(matrixStats)
library(WGCNA)
allowWGCNAThreads() 



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

#now we get to see our modules!
mergedColors2 <- labels2colors(net2$colors) 
plotDendroAndColors(net2$dendrograms[[1]], mergedColors2[net2$blockGenes[[1]]], "Module colors",
                                                             dendroLabels = FALSE, hang = 0.03,
                                                             addGuide = TRUE, guideHang = 0.05)
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
corPvalueStudent(moduleTraitCor2, nrow(dataExpr2))

#Select the geneModule that has highest correlation
module_high_cor2 <- moduleTraitCor2[order(abs(moduleTraitCor2),decreasing = T),]
module_high_cor2
module_high_cor_name2 <- sub('..', '', names(module_high_cor2[1]))

biomarker <- mymodules2[mymodules2$moduleColors2 %in% module_high_cor_name2,"genes2"]
biomarker

# Determine if there is overlap between the DEGs and the module
DE_overlap2 <- gene.list %in% biomarker
print(sum(DE_overlap2))
print(gene.list[DE_overlap2])

# write the genes from the highest correlated modules to a text file
write.table(biomarker, "../Results/biomarker.txt", row.names=F, col.names=F, quote=F)


#GSEA overlap results (FDR < 0.05):\
# utilized GSEA/MSigDB webtool: https://www.gsea-msigdb.org/gsea/msigdb/human/annotate.jsp
# GSEA overlap results - significant gene sets table: refer to 'CR_WGCNA_overlap1.png' in Results folder
# GSEA overlap results - gene/geneset overlap matrix: refer to 'CR_WGCNA_overlap2.png' in Results folder




## ssGSEA and survival analysis


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




## Random forest prediction model


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

# View feature importance
importance <- rf_model$importance
importance <- importance[order(importance, decreasing=T)]
print(importance)

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

# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))

# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)

# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(training_data$Response))
print(conf_matrix)


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
paste("Mean sensitivity =", mean_sens)
paste("Mean specificity =", mean_spec)
paste("Mean AUC =", mean_auc)


# Validate model on independent dataset


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

# Biomarker after removing gene IDs not present in the validation data
final.list <- gene.list[gene.list %in% rownames(val_eset)]

# Make a table containing the independent variables (gene signature expression values) and dependent variables (survival info) for each patient
val_data <- as.data.frame(t(val_eset[final.list,]))
val_data$Response <- as.numeric(factor(val_pdata$Response, levels=c("No response", "Response")))-1
head(val_data)

# Compare validation expression distributions to training data
par(mfrow=c(1,1))
boxplot(t(training_data[, final.list]))
boxplot(t(val_data[, final.list]))

# Z-score transformation of validation data (RNA-seq)
val_data[, final.list] <- scale(val_data[, final.list])
head(val_data)

# Remove gene that could not be scaled from the biomarker
final.list <- final.list[-5]
val_data <- val_data[,c(final.list,"Response")]

# Print final biomarker genes
print(final.list)



# Validation set does not contain all of the genes in the gene list
training_data <- as.data.frame(t(neo_eset[final.list,]))
training_data$Response <- neo_pdata$CR2-1

# Z-score transformation of training data (microarray)
training_data[, final.list] <- scale(training_data[, final.list])
par(mfrow=c(1,1))
boxplot(t(training_data[, final.list]))
boxplot(t(val_data[, final.list]))

library(randomForestSRC)

# Fit a Random Forest model
set.seed(123)
rf_model <- rfsrc(Response ~ ., 
                  data = training_data, 
                  ntree = 500,  # Number of trees
                  importance = TRUE)  # Enables variable importance calculation

# View model summary
print(rf_model)

# View feature importance
importance <- rf_model$importance
importance <- importance[order(importance, decreasing=T)]
print(importance)

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

# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))

# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)

# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(training_data$Response))
print(conf_matrix)


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
paste("Mean sensitivity =", mean_sens)
paste("Mean specificity =", mean_spec)
paste("Mean AUC =", mean_auc)



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

# Identify the best threshold
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")

# Print the optimal threshold
print(paste("Optimal Threshold =", optimal_threshold))

# Extract sensitivity and specificity at the optimal threshold
optimal_coords <- coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"), best.method = "youden")

# Print sensitivity and specificity
print(optimal_coords)

# Predicted classes based on optimal threshold
predicted_classes <- ifelse(class_probs > as.numeric(optimal_threshold), 1, 0)

# Tabulate confusion matrix
conf_matrix <- confusionMatrix(as.factor(predicted_classes), as.factor(val_data$Response))
print(conf_matrix)



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
