#SRR14294834 Single Cell RNA-Seq Data Anaysis

library(Seurat)
library(dplyr)

# Open data count matrix -- excluded by Linux processes
countmat = read.table('counts.txt', header=TRUE)
dim(countmat)
View(countmat)

# Detect first column as a row name
rownames(countmat) = countmat[,1]
View(countmat)
countmat = countmat[,-1]
View(countmat)

# 1) First filtering: Every genes must been in min 3 cells. Every cells must have min 200 genes
srobj = CreateSeuratObject(countmat, project = 'SRR14294834', min.cells = 3, min.features = 200) 
srobj
srobj[['RNA']]@counts

# 2) Second filtering: Removing cells with more than 1 gene subset
# Find Mitochondria genes percent in each cells
srobj[['MTpercent']] = PercentageFeatureSet(srobj, pattern = '^MT-')
srobj$MTpercent
# Violin Plot
# We can also save all the plots in a pdf etc., formats
pdf(file='C:/Users/Newsha/Downloads/SingleCellBehrouz/Plots/SRR14294834_plots.pdf', width = 12, title = "SRR14294834")
VlnPlot(srobj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent", ncol(3)))
# Scatter plot
plot1 = FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "MTpercent")
plot2 = FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Removing cells with below min 200 and above max 800 genes and cells with above 16 mitochondrial genes
srobj = subset(srobj, subset=nFeature_RNA>200 & nFeature_RNA<800 & MTpercent<16)
dim(srobj)
# Check again Plots
VlnPlot(srobj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent", ncol(3)),cols = 3)
plot1 = FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "MTpercent", cols = 3)
plot2 = FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = 3)
plot1 + plot2
dev.off()

# 3) Normalizing Data
srobj = NormalizeData(srobj, normalization.method = "LogNormalize", scale.factor = 10000)
srobj[["RNA"]]@data

# 4) Highest Variation Genes:
library(ggplot2)
srobj = FindVariableFeatures(srobj, selection.method = "vst", nfeatures = 2000)
topVariableGenes = head(VariableFeatures(srobj), 15)
plot1 = VariableFeaturePlot(srobj)
plot2 = LabelPoints(plot=plot1, points = topVariableGenes, size=3, ynudge = 0.06, xnudge = -0.08)
plot2

# 5) Centering and Scaling Data (Mean=0, Variance=1)
srobj = ScaleData(srobj, features = rownames(srobj))
dim(srobj[["RNA"]]@scale.data)
mean(srobj[["RNA"]]@scale.data[,1])
var(srobj[["RNA"]]@scale.data[,1])

# 6) Dimension Reduction
srobj = RunPCA(srobj, features = VariableFeatures(srobj))
print(srobj[['pca']], dims = 1:10, nfeatures = 5)
VizDimLoadings(srobj, dims = 1:2, reduction = "pca", ncol = 3)
DimPlot(srobj, reduction = "pca")
DimHeatmap(srobj, dims = 1, cells = 300)
ElbowPlot(srobj)

# 7) Determine Number of Dimensions
srobj = JackStraw(srobj, num.replicate = 100)
srobj = ScoreJackStraw(srobj, dims = 1:20)
JackStrawPlot(srobj, dims = 1:20)
ElbowPlot(srobj)

# 8) Clustering (K Nearest Neighbors):
srobj = FindNeighbors(srobj, dims = 1:13)
srob = FindClusters(srobj, resolution = 0.5)
Idents(srobj)
levels(Idents(srobj))
length(levels(Idents(srobj)))
View(srobj)

# 9) Non-Linear Dimension Reduction
srobj = RunTSNE(srobj, dims = 1:13)
DimPlot(srobj, reduction = "tsne")
# srobj = RunUMAP(srobj, dims: 1:13)
# DimPlot(srobj, reduction = "umap")