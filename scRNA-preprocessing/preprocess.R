set.seed(42) #ensures same results for any randomization steps. helps to keep reproducibility

library(tidyverse)
library(Seurat)


nsclc_sm <- Read10X_h5(filename="./40k_NSCLC_DTC_3p_HT_nextgem_donor_1_count_sample_feature_bc_matrix.h5")

#Genome matrix has multiple modalities, returning a list of matrices for this genome
#1 has gene expression data, 1 has antibody capture, 1 has multiplexing capture

cts <- nsclc_sm$'Gene Expression' #sparse matrix. 


#seurat object. use count matrix, give it a name, min # of cells, and min # of features
#cells with at least 3 features (3 genes), genes that are at leadt expressed in 200 cells
nsclc_seu <- CreateSeuratObject(counts = cts, project = "NSCLS Preprocess", min.cells = 3, min.features =  200 )

#QC and Filtering 
#Remove cells with low # of genes detected
#Remove cells with very high gene count
#Remove cells with high mt genes 

nsclc_seu[['precent_mt']] <- PercentageFeatureSet(nsclc_seu, pattern= "^MT-") #humans is uppercase, mouse is lowercase

View(nsclc_seu@meta.data) #nCount_RNA- # of RNA molecules in the cell. nFeatures_RNA is # of unique genes detected in the cell
#Low nFeauture_RNA for a cell indicates that it may be dead or empty droplet. High nCount or nFeature suggests doublets or multiplets.  

VlnPlot(nsclc_seu, features = c("nCount_RNA", "nFeature_RNA", "precent_mt"), ncol=3 )


#Look at # of molecules per sample vs number of genes. Should follow a linear pattern. The lower right corner cells have high count and low features. Top left would mean we have high # of genes but not deeply sequenced enough
FeatureScatter(nsclc_seu, feature1="nCount_RNA", feature2="nFeature_RNA") + geom_smooth(method='lm') #add line

#filter out too low and too high features and counts

nsclc_seu <- subset(nsclc_seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & precent_mt < 5)

#Normalize data so to remove systematic variation, reduce noise, and improve the interpretability of your data
#global scaling that normalize feature expression for each cell by total expression and then multiply by 10000 factor and then log transforms the results. These are the default params
nsclc_seu <- NormalizeData(nsclc_seu, normalization.method = "LogNormalize", scale.factor = 10000)


str(nsclc_seu) #we can see the normalized data



#Identify features that have high variability between cells. Expressed in some and not other. Distinguish subset of cells from others. Focus on these genes helps highlight biological signal of scRNA data analysis

nsclc_seu <- FindVariableFeatures(nsclc_seu, selection.method = "vst", nfeatures = 2000) #defualt is 2000. Variance Stabilizing Transformation

top10 <- head(VariableFeatures(nsclc_seu),10)

top10_plot <- VariableFeaturePlot(nsclc_seu )
LabelPoints(top10_plot, points=top10, repel=TRUE) #black less variable

#Scale data - remove unwanted sources of variation. Technical variation sunch as batch effect, cell cycle stage, mt contamination,etc...

all_genes <- rownames(nsclc_seu)
nsclc_seu <- ScaleData(nsclc_seu, features = all_genes)

nsclc_seu <- RunPCA(nsclc_seu, features= VariableFeatures(nsclc_seu))

#Each PC will have loadings for each gene which tells you how much that gene explains the variation in the dataset PC1 explains more variation than PC2 and so on...

print(nsclc_seu[['pca']], dims =1:3 , nfeatures =5) # first 3 pcs. first +/- 5 features


DimHeatmap(nsclc_seu, dims=1,cells = 500, balanced = TRUE) #pc1, 500 cells.  PC1 already groups cells according to principal component. We look into a group of genes.


DimPlot(nsclc_seu, reduction = "pca") #each dot is a cell. Cells group into clusters.

#Elbow plot allows us to determine # of PCs we need to explain enough variation

ElbowPlot(nsclc_seu) #about 15-20 is enough

#Clustering: Identify groups or clusters of similar cells based on their gene expression profiles.
nsclc_seu <- FindNeighbors(nsclc_seu, dims = 1:20)

# Louvain
nsclc_seu <- FindClusters(nsclc_seu, resolution = c(.1, .3, .5, .7, 1)) #different resolutions. The lower the resolution the less clusters. High resolutions leads to more clusters


DimPlot(nsclc_seu, group.by = "RNA_snn_res.0.1", label=TRUE)


#UMAP: Visualize the high-dimensional data in 2D or 3D while preserving local and global relationships between cells.


nsclc_seu <- RunUMAP(nsclc_seu, dims = 1:20)
DimPlot(nsclc_seu, reduction = 'umap')

#next we can look into genetic markers to identify each cluster. 


#save seurat object with saveRDS
#saveRDS("nscls_seu")


#plot with cluster labels. Clustered close based on gene expression profile
DimPlot(nsclc_seu, reduction = 'umap', label=TRUE)

#Find markers for cluster 2

cluster2.markers <- FindAllMarkers(nsclc_seu)

#avg_log2FC = 2.69 → This means NKG7 is upregulated in cluster 0 compared to other clusters. A log2 fold change of 2.69 suggests that its expression is significantly higher in this cluster.
#95.3% of the cells in cluster 0 express NKG7
FeaturePlot(nsclc_seu, features = "NKG7") + 
  DimPlot(nsclc_seu, reduction = "umap", label = TRUE, repel = TRUE)


FeaturePlot(nsclc_seu, features = "NKG7", label = TRUE)  
VlnPlot(nsclc_seu, features = "NKG7", group.by = "seurat_clusters")  

#we see that cells 0,1,3 highly express NKG7 which may serve as a reason for why they cluster together



