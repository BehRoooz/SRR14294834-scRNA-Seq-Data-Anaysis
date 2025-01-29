#!/usr/bin/env Rscript

#===============================================================================
# Single Cell RNA-Seq Data Analysis Pipeline
# Dataset: SRR14294834
#===============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)

#-------------------------------------------------------------------------------
# 1. Data Import and Initial Processing
#-------------------------------------------------------------------------------

# Import count matrix
countmat <- read.table('counts.txt', header = TRUE)
message("Initial dimensions: ", paste(dim(countmat), collapse = " x "))

# Set row names and remove redundant first column
rownames(countmat) <- countmat[,1]
countmat <- countmat[,-1]

#-------------------------------------------------------------------------------
# 2. Create Seurat Object with Initial Filtering
#-------------------------------------------------------------------------------

# Create Seurat object with minimal filtering criteria
# Minimum 3 cells per gene, minimum 200 genes per cell
srobj <- CreateSeuratObject(
    countmat,
    project = 'SRR14294834',
    min.cells = 3,
    min.features = 200
)

#-------------------------------------------------------------------------------
# 3. Quality Control
#-------------------------------------------------------------------------------

# Calculate mitochondrial gene percentage
srobj[['MTpercent']] <- PercentageFeatureSet(srobj, pattern = '^MT-')

# Generate QC plots
pdf('SRR14294834_QC_plots.pdf', width = 12)

# Pre-filtering QC plots
VlnPlot(srobj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)

# Feature correlation plots
plot1 <- FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "MTpercent")
plot2 <- FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

# Filter cells based on QC metrics
srobj <- subset(srobj,
    subset = nFeature_RNA > 200 &
            nFeature_RNA < 800 &
            MTpercent < 16
)

# Post-filtering QC plots
VlnPlot(srobj, features = c("nFeature_RNA", "nCount_RNA", "MTpercent"), ncol = 3)
plot1 <- FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "MTpercent")
plot2 <- FeatureScatter(srobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)

dev.off()

#-------------------------------------------------------------------------------
# 4. Normalization and Feature Selection
#-------------------------------------------------------------------------------

# Normalize data
srobj <- NormalizeData(
    srobj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

# Find variable features
srobj <- FindVariableFeatures(
    srobj,
    selection.method = "vst",
    nfeatures = 2000
)

# Plot variable features
pdf('SRR14294834_variable_features.pdf', width = 10)
top_variable_genes <- head(VariableFeatures(srobj), 15)
plot1 <- VariableFeaturePlot(srobj)
plot2 <- LabelPoints(
    plot = plot1,
    points = top_variable_genes,
    size = 3,
    ynudge = 0.06,
    xnudge = -0.08
)
print(plot2)
dev.off()

#-------------------------------------------------------------------------------
# 5. Data Scaling
#-------------------------------------------------------------------------------

# Scale data (mean = 0, variance = 1)
srobj <- ScaleData(srobj, features = rownames(srobj))

#-------------------------------------------------------------------------------
# 6. Dimension Reduction and Analysis
#-------------------------------------------------------------------------------

# Run PCA
srobj <- RunPCA(srobj, features = VariableFeatures(srobj))

# Generate PCA analysis plots
pdf('SRR14294834_PCA_analysis.pdf', width = 12)
print(srobj[['pca']], dims = 1:10, nfeatures = 5)
VizDimLoadings(srobj, dims = 1:2, reduction = "pca", ncol = 3)
DimPlot(srobj, reduction = "pca")
DimHeatmap(srobj, dims = 1, cells = 300)
ElbowPlot(srobj)

# Determine significant dimensions
srobj <- JackStraw(srobj, num.replicate = 100)
srobj <- ScoreJackStraw(srobj, dims = 1:20)
JackStrawPlot(srobj, dims = 1:20)
ElbowPlot(srobj)
dev.off()

#-------------------------------------------------------------------------------
# 7. Clustering Analysis
#-------------------------------------------------------------------------------

# Perform clustering
srobj <- FindNeighbors(srobj, dims = 1:13)
srobj <- FindClusters(srobj, resolution = 0.5)

# Generate t-SNE visualization
srobj <- RunTSNE(srobj, dims = 1:13)

pdf('SRR14294834_clustering.pdf', width = 10)
DimPlot(srobj, reduction = "tsne")
dev.off()

# Optional UMAP visualization
# srobj <- RunUMAP(srobj, dims = 1:13)
# DimPlot(srobj, reduction = "umap")

#-------------------------------------------------------------------------------
# Save processed object
#-------------------------------------------------------------------------------

saveRDS(srobj, file = "SRR14294834_processed.rds")

# Print final cluster information
message("Number of clusters identified: ", length(levels(Idents(srobj))))
