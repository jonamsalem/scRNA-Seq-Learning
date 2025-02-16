# Here’s the corrected version:  
#   
#   **When to integrate?** When you have multiple RNA-seq datasets from different samples, treatments, etc. You can also use integration as a reference with a queried dataset. It can also be used to merge modalities.  
# 
# **We have different patients → integrate and correct for batch effects.**  
#   
#   - **Horizontal integration** – Same modality measured across different cells (e.g., RNA-seq performed on the same tissue but different patients).  
# - **Vertical integration** – Different modalities from the same group of cells (e.g., methylation data, scATAC-seq data).  
# - **Diagonal integration** – Different modalities from different cells.  
# 
# We will perform **horizontal integration** since we have the same modality but different cells.  
# 
# **What are batch effects?**  
#   Technical variations between datasets that need to be corrected for accurate comparisons.  
# 
# **Methods to correct:** Seurat v3, MNN, Harmony,etc



library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)


#get data location

dirs <- list.dirs(path="scRNA-Integration/feature_data/", recursive = F, full.names = F) #7 directories


#create a for loop to create count matrix and create a seurat object that contains patient info and tissue type

#extract the count matrix, barcodes, and features for each foler, create a count matrix, and then create a seurat object
for(x in dirs){
  name <- gsub('filtered_feature_bc_matrix', "" ,x)
  
  cts <- ReadMtx(mtx= paste0('scRNA-Integration/feature_data/', x, "/matrix.mtx.gz"), features = paste0("scRNA-Integration/feature_data/", x, "/features.tsv.gz"),cells = paste0("scRNA-Integration/feature_data/", x, "/barcodes.tsv.gz") )
  
  assign(name, CreateSeuratObject(counts = cts))
}


#merge all objects so we can perform QC and filtering simoultaneously

merged_seurat <- merge(
  HB17_background_, 
  y = c(HB17_PDX_, HB17_tumor_, HB30_PDX_, HB30_tumor_, HB53_background_, HB53_tumor_),
  add.cell.ids = c("HB17_background", "HB17_PDX", "HB17_tumor", "HB30_PDX", "HB30_tumor", "HB53_background", "HB53_tumor"),
  project = "HB"
)

View(merged_seurat@meta.data)

#add column sample with sample name to the metadata
merged_seurat$sample <- rownames(merged_seurat@meta.data)

#split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col='sample', into=c("Patient", "Type", "Barcode"))


#ensure we have data from all patients and tissue type
unique(merged_seurat$Patient)
unique(merged_seurat$Type)

#QC and filtering

merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')


#filtering cells with atl eadt 800 tracripts and at least 500 genes and mitoPercent less than 10%
merged_seurat <- subset(merged_seurat, subset = nCount_RNA > 800 & nFeature_RNA > 500 & mitoPercent < 10)



#check if integration is required. Explore data in low dim space

merged_seurat <- NormalizeData(object = merged_seurat ) #adjust for sequencing depth between cells
merged_seurat <- FindVariableFeatures(object = merged_seurat) #Identify most variable genes across all cells. Default is 2000
merged_seurat <- ScaleData(merged_seurat) #standardize the genes to remove noise from varying gene expression level and ensures that each gene contributes equally to wownstream analyis.

merged_seurat <- RunPCA(merged_seurat) #Reduce dimensionality of data for visualization. Assumes linear-relationship between varaibles.


ElbowPlot(merged_seurat) #Use all 20

merged_seurat <- FindNeighbors(merged_seurat, dims=1:20) #Identify cells similar to each other based on gene expression profile
merged_seurat<- FindClusters( merged_seurat) #Group cells into clusters based on values from FindNeighbors. 
merged_seurat <- RunUMAP(merged_seurat, dims= 1:20) #Reduce dimensionality to visualize clusters. Does not assume linearity. 



#plot  cells by patient
p1<-DimPlot(merged_seurat, reduction='umap', group.by='Patient') #Cells from same patient cluster together ---> batch effect (cluster by non-biological effect)

#plot cells by tissue type
p2<-DimPlot(merged_seurat, reduction='umap', group.by='Type', cols=c('red', 'green', 'blue'))
            
grid.arrange(p1,p2, ncol=2, nrow=2)


#plot on left shows us the clustering we see is due to technical varaition and not biological differences. cells from different cells cluster differently and masking biological variation clustering. We need to corrent for batch effects 

onb_list <- SplitObject(merged_seurat, split.by = 'Patient') #split by patient



#run normalization and select features for each object

for (i in 1:length(onb_list)){
  onb_list[[i]] <- NormalizeData(onb_list[[i]]) #normalize again bc  patient set is not the same
  onb_list[[i]] <- FindVariableFeatures(onb_list[[i]]) #find variable features again bc most variable genes may change when using a subset
}


#select integration features

features <- SelectIntegrationFeatures(onb_list) #Identify genes across patient subset.

#find integration anchors (CCA). Can try MNN instead.
anchors <- FindIntegrationAnchors(onb_list, anchor.features=features) #Anchors are pairs of cells one from each dataset/subset that are similar in terms of gene expression. Helps to align dataset by shared biological features. 

#integrate data
seurat.integrated <- IntegrateData(anchorset=anchors) #Integration done based on integration pairs. Full data is still includes.

seurat.integrated <- ScaleData(seurat.integrated) #Scale again after integration bc integration could change underlying structure of the data. 

seurat.integrated <- RunPCA(seurat.integrated) 

seurat.integrated <- RunUMAP(seurat.integrated, dims=1:50) #Visualize after integration

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type')

grid.arrange(p3, p4, ncol =2)


