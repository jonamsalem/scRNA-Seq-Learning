# Harmony ensures that cells in each clusters come from as many batches as possible. Diversify batches and calculates correction factor. Iterates until convergence

# Harmony vs seurat CCA. CCA uses anchor and returns an integrated object. Harmony only computes corrected dimensionality reduction values and does not return integrated object.


# CCA - anchors between datasets and integrates into new expression matrix. Harmony corrects PCA embeddings does not change seurat data. 
# Harmony is gaster and works well for datasets with multiole batches or conditions. CCA is good for strong batch effects but is slow
# Harmony is better for large datasets and complex batch effects where as CCA is useful when strong biological caraitions needs alignment across datasets



# These methods are not only used for batch correction. They can be used for inter-sample varaition due to differneces in seq tech or experimental condition

#We will use harmony to integrate data from differenet conditions

set.seed(1234)
library(Seurat)
library(harmony)
library(SeuratData)
library(tidyverse)
library(ggplot2)

AvailableData() #see availalbe dataset. We use ifnb


InstallData("ifnb") #install the data


LoadData("ifnb") #already a seurat object


str(ifnb) #we see the data is not processed bc not scaling, reduction, features

ifnb = UpdateSeuratObject(object = ifnb) #update object


#Process data using standard qc
ifnb$mito <- PercentageFeatureSet(ifnb, pattern = '^MT-')

View(ifnb@meta.data)

#filtering
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 & nFeature_RNA > 200 & mito < 5 ) 

#Normalize, Feature selection, scale, and run PCA
ifnb.filtered <- NormalizeData(ifnb.filtered) #comparable gene expression values. Log normalization

ifnb.filtered <- FindVariableFeatures(ifnb.filtered)

ifnb.filtered <- ScaleData(ifnb.filtered) # All genes contribute equally to PCA. Scale after find variable features bc they carry the most important biological signals

ifnb.filtered <- RunPCA(ifnb.filtered)

ElbowPlot(ifnb.filtered) #use about 15-20

ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim'  ) #cells cluster based on condition. Separation might be too strong. You want cell types to drive clusterings.

before_cell <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'seurat_annotations')

#run harmony
ifnb.harmony<- ifnb.filtered %>% #tidyverse syntax. take the filtered object and passes it to run harmony. Correct stim varaible. Do not plot plot convergence
  RunHarmony(group.by.vars='stim', plot_convergence= FALSE)

ifnb.harmony@reductions #harmony embeddings saved in a subslot called harmony

ifnb.harmony.harmony.embed <- Embeddings(ifnb.harmony, "harmony") #embedding values for each cell. Embeddings are the coordinates of each cell in the dimensional reduction space.


ifnb.harmony <-ifnb.harmony %>%
  RunUMAP( reduction= 'harmony', dims = 1:20) %>% #use harmony instead of PCA using the new dimensinality reduction
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = .5) #Assign cluster labels 


after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim') #cells now do not cluster by condition
after_cell  <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'seurat_annotations') #cells now do not cluster by condition

before | after



#clusters match the cell type
cluster <-DimPlot(ifnb.harmony, group.by = 'RNA_snn_res.0.5') 
cluster2 <- DimPlot(ifnb.harmony, group.by = 'seurat_annotations') 

before_cell | after_cell #harmony correction

cluster | cluster2 #cells cluster by cell type
