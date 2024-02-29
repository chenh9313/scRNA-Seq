#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(patchwork)
###Step1: Load the PBMC dataset
#The read in data under sampleSAMN15315807/outs/filtered_feature_bc_matrix
#they are: barcodes.tsv  genes.tsv  matrix.mtx
pbmc.data <- Read10X(data.dir = "/Desktop/SingleCell/Example_data_Seurat/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size

###Step2: QC and selecting cells for further analysis
#The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#FeatureScatter is typically used to visualize feature-feature relationships, but can be used
#for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#We filter cells that have unique feature counts over 2,500 or less than 200 and cells that have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
dim(pbmc)
#But what should be the min gene and max gene? 
#Louoe Browser suggests UMI counts between 2^10 and 2^14;
# Louoe Browserbarcodes with gene numbers between 2^9.5 and 2^13; 
#Louoe Browser Enter 0 and 15% for the percentage of mitochondrial genes with barcodes
test <- subset(pbmc, subset = nFeature_RNA > 724 & nFeature_RNA < 8192 & nCount_RNA > 1024 & nCount_RNA < 16384 & percent.mt < 5)
dim(test)

###Step3: Normalizing the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#pbmc <- NormalizeData(pbmc)
pbmc
#Identification of highly variable features (feature selection); downstream analysis, like PCA
#high cell-to-cell variation in the dataset(i.e, they are highly expressed in some cells, and lowly expressed in others)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

###Step4: Scaling the data
#apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#remove unwanted sources of variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt")#strongly suggest to use this one
#pbmc

###Step5: Perform linear dimensional reduction
#perform PCA on the scaled data
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:10, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:3, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
#change dims=1-10 to check heterogeneity in dataset and decide which PCs to include for further downstream analyses
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

###Step6: Determine the ‘dimensionality’ of the dataset
#Q:how many components should we choose to include? 10? 20? 100
#identify ‘significant’ PCs as those who have a strong enrichment of low p-value features by randomly permute a subset of the data (1% by default) and rerun PCA, 
#constructing a ‘null distribution’ of feature scores, and repeat this procedure.
#NOTE:This process can take a long time for big datasets, comment out for expediency. More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
#more PC will help revel the rare cell, such as, immune subsets
pbmc <- JackStraw(pbmc, num.replicate = 100) #takes a long time to run
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15) #comparing the distribution of p-values for each PC
ElbowPlot(pbmc)

###Step7: Cluster the cells
#K-nearest neighbor (KNN) graph
#setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

###Step8:Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "/Users/chenhua9/Desktop/SingleCell/Example_data_Seurat/pbmc_tutorial.rds")

###Step9: Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2; identifies positive and negative markers of a single cluster 
#ident.1 means which cluster(1,2,3...10); 
#min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
tail(cluster2.markers, n = 5)
dim(cluster2.markers)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25, only.pos = TRUE)
head(cluster5.markers, n = 5)
dim(cluster5.markers)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

#visualizing marker expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#Feature marker expression
FeaturePlot(pbmc, features = c("MS4A1"))
#FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
#expression heatmap for given cells and features
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

###Step10: Assigning cell type identity to clusters
#you have to know the marker genes of each group to name the clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "/Users/chenhua9/Desktop/SingleCell/Example_data_Seurat/pbmc3k_final.rds")


