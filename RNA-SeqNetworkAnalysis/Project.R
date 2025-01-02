rm(list=ls())
dev.off()
gc()

library(edgeR)
library(cluster)
library(mclust)
library(RColorBrewer)
library(ggplot2)
library(knitr)
library(limma)
library(reshape2)
library(RColorBrewer)
library(WGCNA)
library(gplots)
library(AnnotationDbi)
library(org.Dm.eg.db)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(100)

female_counts = read.csv("female_counts.csv", row.names = 1)
female_meta = read.csv("bio_inf_meta_female.csv")
colnames(female_meta)[1] = "Sample"

# Data preprocessing

# Filter out rows with < 12 low counts
y_f <- female_counts[rowSums(female_counts)>12,]

# DGEList
x_f <- DGEList(counts = y_f, group = female_meta$Treatment)
x_f

# Normalization
x_f <- calcNormFactors(x_f, method = "TMM")
x_f$samples

win.graph(700, 700)
plotMDS(x_f)

# Estimate dispersion
x_f <- estimateCommonDisp(x_f, verbose = TRUE)
x_f <- estimateTagwiseDisp(x_f)
x_f <- estimateTrendedDisp(x_f)

win.graph(700, 700)
plotBCV(x_f)

# Testing
f_test <- exactTest(x_f)
f_top <- topTags(f_test, n = 50, adjust.method = "fdr")
f_top

f_tests <- decideTests(f_test)
f_tests

# Differential express
sum(f_tests[,1] != 0) # Total differential express 875 genes

DE_f = which(f_tests[,1] != 0) # Index
allDE_f = f_tests[DE_f, ]
summary(allDE_f) # Down 175 Up 700

DE_f_genes <- x_f[rownames(allDE_f),]
head(DE_f_genes$count)

# PCA
log_DE_f = log2(x_f$counts + 1/6)
pca = prcomp(t(log_DE_f))
summary(pca)

win.graph(700, 700)
plot(pca$sdev, type = "l", ylab = "Standard Deviation", xlab = "Number of Components")

win.graph(700, 700)
plot(predict(pca), col = c("red", "red", "red", "red",
                        "blue", "blue", "blue", "blue"), pch = 16)
legend("topright", legend = c("HS", "LS"), col = c("red", "blue"), pch = 20,
       cex = 1)

# Clustering
win.graph(700, 700)
heatmap(log_DE_f, cexCol = 0.8,
        col = colorRampPalette(brewer.pal(8, "Blues"))(25), margins = c(5, 15))
legend(x = "topright", legend = c("Downregulated", "Upregulated"), 
       fill = colorRampPalette(brewer.pal(8, "Blues"))(3))

# Model based clustering
cluster.BIC <- mclustBIC(log_DE_f)
summary(cluster.BIC)
log_DE_f_clust = Mclust(log_DE_f, x = cluster.BIC)
log_DE_f_clust$G
log_DE_f_clust$modelName

clust_class_f <- log_DE_f_clust$classification

# Network analysis
datExpr = t(y_f)

# Clustering the samples for obvious outliers
sampleTree = hclust(dist(datExpr), method = "average")
par(cex = 0.6)
win.graph(700, 700)
plot(sampleTree, main = "Sample Clustering to Detect Outlier", sub = "",
     xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
enableWGCNAThreads()

powers = c(c(1:10), seq(from = 12, to = 30, by = 2)) # Range of power
# Tunning soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate
win.graph(1400, 700)
par(mfrow = c(1,2))

# Scale-free topology fit index as a function of the soft-thresholding power 
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold Power",
     ylab = "Scale Free Topology Model Fit Signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9,col = "red")
abline(h = 0.9,col = "red")

# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold Power",ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, 
     cex = 0.9,col = "red")
# No optimal solution

charData = apply(as.matrix(datExpr), 2, as.character)
numData = apply(datExpr,2, as.numeric)

cor <- WGCNA::cor # Conflict between WGCNA with other software packages
net = blockwiseModules(numData, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femalefflyTOM",
                       verbose = 3)
table(net$colors)
cor <- stats::cor

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
win.graph(1400, 700)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]

traitData = female_meta
dim(traitData)
names(traitData)

femaleSamples <- rownames(datExpr)
traitRows <- match(femaleSamples, traitData$Sample)
datTraits <- as.data.frame(traitData[traitRows, -1])
rownames(datTraits) <- traitData[traitRows, 1]
dim(datTraits)

Treatment <- traitData$Treatment
datTraits <- as.data.frame(model.matrix(~0+Treatment))
# Recalculate MEs with color labels
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
win.graph(1400, 700)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

highSug <- as.data.frame(datTraits$TreatmentHS)
names(highSug) = "High Sugar"
# Names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),
                                          nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(datExpr, highSug, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),
                                          nSamples))

names(geneTraitSignificance) = paste("GS.", names(highSug), sep = "")
names(GSPvalue) = paste("p.GS.", names(highSug), sep = "")

module = "salmon"
column = match(module, modNames)
moduleGenes = moduleColors == module
win.graph(1400, 700)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance from High Sugar Diet",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "darkred"
column = match(module, modNames)
moduleGenes = moduleColors == module

win.graph(1400, 700)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance from High Sugar Diet",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Similarity matrix by Pearson correlation and Euclidean Distance
similarity <- function(data) {
        cor_matrix  <- cor(t(data)) # Pearson correlation
        dist_matrix <- as.matrix(dist(data, diag = TRUE, upper = TRUE))
        # Eucliden Distance
        dist_matrix <- log1p(dist_matrix) # Sigmoid Transformation
        dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
        
        sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix) / 2)
}

sim_matrix <- similarity(DE_f_genes$count)

win.graph(700, 700)
heatmap.2(t(sim_matrix), col = redgreen(75), labRow = NA, labCol = NA,
          trace = 'none', dendrogram = 'row', xlab = 'Gene', ylab = 'Gene',
          main = 'Similarity Matrix', density.info='none', revC = TRUE)

# Tune Soft Thresholding Power
powers = c(seq(1,100, by = 1)) 
sft = pickSoftThreshold(sim_matrix, powerVector = powers, networkType = "signed",
                        verbose = 2, dataIsExpr = T)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold", ylab = "Scale Free Topology Model Fit", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,cex = 0.9,col = "red")

# Construct adjacency matrix
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power = 9, type = 'signed')

rm(sim_matrix)
gc()

gene_ids <- rownames(adj_matrix)

adj_matrix <- matrix(adj_matrix, nrow = nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids

win.graph(700, 700)
heatmap.2(t(adj_matrix), col = redgreen(75), labRow = NA, labCol = NA, trace = 'none',
          dendrogram = 'row', xlab = 'Gene', ylab = 'Gene', main = 'Adjacency Matrix',
          density.info = 'none', revC = TRUE)

distance = as.dist(1 - adj_matrix)
gene_tree <- hclust(distance, method = "average")
win.graph(700, 700)
plot(gene_tree)
# Break apart the Hierarchical Clustering into modules
module_labels <- cutreeDynamicTree(dendro = gene_tree, minModuleSize = 15,
                                   deepSplit = TRUE)

# Assign color to each module
module_colors <- labels2colors(module_labels)

source("export.R")

# Retrieve gene annotations
gene_info <- select(org.Dm.eg.db, keytype = 'FLYBASE', keys = rownames(DE_f_genes$count),
                    columns = c('GENENAME', 'ALIAS', 'SYMBOL', 'ENTREZID'))

# Grab the description for the first transcript
gene_info <- gene_info[!duplicated(gene_info$FLYBASE), ]
gene_info <- cbind(gene_info, module = module_colors)
gene_info$color_rgb <- col2hex(gene_info$module)

# A lot of low correlation
sum(adj_matrix >= 0.001)
g <- export_network_to_graphml(adj_matrix, filename = 'network.graphml',
                               threshold = 0.4, nodeAttrDataFrame = gene_info)
g