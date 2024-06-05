#load packages
library(WGCNA)
library(tibble)
library(ggplot2)
library(tidyverse)
library(dplyr)

#do not omit this step
options(stringsAsFactors = FALSE)

#Read in the dataset
data=read.delim("~/Desktop/norm_counts.txt", header=TRUE, row.names=1, sep="\t")
data <- data[,-c(17,25)]
data<-data[rowSums(data == 0) == 0, ]

#take a look in the dataset
dim(data)
names(data)

#make the dataframe
datExpr0 = data.frame(data)
#transpose the data
datExpr = as.data.frame(t(datExpr0))
#remove the X
rownames(datExpr) = sub("*\\X", "", rownames(datExpr))
#check the genes and samples for too many missing values
gsg = goodSamplesGenes(datExpr, verbose =3)
gsg$allOK

#cluster the samples to look for any outliers
sampleTree = hclust(dist(datExpr), method = "average")
#plot the sample tree
par(cex = 0.6)
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

datExpr1 = datExpr
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)
#datExpr1 now has the expression data ready for network analysis

#loading clinical data
traitData = read.csv("~/Desktop/pheno.csv", header=TRUE ,sep = ",")
traitData <- traitData[,-c(7:10)]
dim(traitData)
names(traitData)

#form a data frame analogous to expression data that will hold the clinical traits
hcjSamples = rownames(datExpr1)
#match the samples for which they were measured to the expression samples
traitRows = match(hcjSamples, traitData$DMB)
datTraits0 = traitData[traitRows,]
rownames(datTraits0) = traitData[traitRows, 1]

#remove the column-DMB
datTraits0<- datTraits0[,-1]
datTraits = datTraits0

collectGarbage()

#expression data is in the variable- datExpr1
#clinical traits is in the variable- datTraits
#visualize the clinical traits relating to the sample dendrogram
#recluster the samples
sampleTree2 = hclust(dist(datExpr1), method = "average")

#convert traits to color rep, white- low, red- high, grey- missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

#plot the sample dendrogram and the colors underneath

plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

#save the data
save(datExpr1, datTraits, file = "hcjInit.RData")

#network construction-> 1 step
#multithreading
allowWGCNAThreads()


#choose the soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#call the network topology analysis function
sft = pickSoftThreshold(datExpr1, networkType = "signed", powerVector = powers, verbose = 5)

sizeGrWindow(16,9)
par(mfrow = c(1,2))
cex1 = 0.90
#scale-free topology fit index as a function of the sft power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
#the following line will correspond to using an R^2 cut-off of h
abline(h=0.80, col="red")
#mean connectivity as a function of the sft power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0, col="red")

#one-step network construction
temp_cor <- cor       
cor <- WGCNA::cor  # Force it to use WGCNA cor function (fix a namespace conflict issue)
net = blockwiseModules(datExpr1,
                       power = 10,
                       networkType = "signed",
                       deepSplit = 2,
                       pamRespectsDendro = F,
                       minModuleSize = 30,
                       maxBlockSize = 4000,
                       reassignThreshold = 0,
                       mergeCutHeight = 0.25,
                       saveTOMs = T,
                       saveTOMFileBase = "ER",
                       numericLabels = T,
                       verbose = 3)

table(net$colors)
cor <- temp_cor     # Return cor function to original namespace

#convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

#plot the dendro and module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]


#relate modules to external clinical traits
#quantify the module-trait associations
#we have eigengene for each module
#correlate the eigengenes with external traits
#define numbers of genes and samples
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)
#recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#make the graphical representation
sizeGrWindow(10,6)
#will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = " ")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 17.5, 2.5, 1))

#display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.30,
               zlim = c(-1,1))

data[duplicated(row.names(data)), ]

#identifying the gene relationships in modules
#change module to desired module
#gene relationship to trait and important modules - calcium
ca_levels = as.data.frame(datTraits$Serum_Ca_mg_dL)
names(ca_levels) = "ca_levels"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr1, ca_levels, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(ca_levels), sep="");
names(GSPvalue) = paste("p.GS.", names(ca_levels), sep="")

#mm-gs for modules
#all the mm-gs done by replacing module with the respective colours
module = "lightyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for calcium levels",
                   main = paste("Module membership vs. gene Significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")


#gene relationship to trait and important modules - vitamin D
d_levels = as.data.frame(datTraits$Vit_D_ng_mL)
names(d_levels) = "d_levels"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr1, d_levels, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(d_levels), sep="");
names(GSPvalue) = paste("p.GS.", names(d_levels), sep="")

#mm-gs for modules
#all the mm-gs done by replacing module with the respective colours
module = "lightyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for vitamin D levels",
                   main = paste("Module membership vs. gene Significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")


#gene relationship to trait and important modules - eIS scores
is_levels = as.data.frame(datTraits$IS_CACTI_exA)
names(is_levels) = "is_levels"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr1, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr1, is_levels, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(is_levels), sep="");
names(GSPvalue) = paste("p.GS.", names(is_levels), sep="")

#mm-gs for modules
#all the mm-gs done by replacing module with the respective colours
module = "lightyellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;

par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CACTI levels",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#display the gene names inside the modules
#change moduleColors to desired module
dm = colnames(datExpr1)[moduleColors=="lightyellow"]
write.csv(dm, file = "~/Desktop/lyw_trim.csv")
