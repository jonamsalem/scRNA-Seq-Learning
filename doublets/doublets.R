# Doublets in Single-Cell RNA Sequencing  
# A doublet is a droplet that contains two cells instead of one.  
# Doublets can lead to miscalculations or misinterpretations of the data,  
# so they must be removed.  

# Types of Doublets  
# 1. Homotypic Doublets – Derived from transcriptionally similar cells.  
# 2. Heterotypic Doublets – Derived from transcriptionally distinct cells.  

# DoubletFinder is more sensitive to detecting heterotypic doublets.  

# DoubletFinder Parameters  
# DoubletFinder requires three key parameters:  
# pN  = Proportion of artificial doublets introduced (default: 25%)  
# pK  = Neighborhood size used to compute the number of artificial nearest neighbors  
# Exp = Expected number of real doublets  

# DoubletFinder Workflow  
# 1. Simulate artificial doublets from existing data and introduce a proportion (pN) of artificial doublets  
# 2. Merge artificial doublets with the real dataset and perform preprocessing  
# 3. Perform PCA on the merged data – artificial doublets co-localize with real cells  
# 4. Detect nearest neighbors and compute the proportion of artificial nearest neighbors (depends on pK)  
# 5. Use the proportion of artificial nearest neighbors to predict real doublets based on the expected number of doublets (Exp)  

# Finding the Best pK Value  
# 1. If ground truth labels (known doublets) are available, use them to determine the best pK  
# 2. If no ground truth is available, estimate pK using heuristic methods  

# Determining Exp (Expected Doublet Rate)  
# - Exp can be estimated from the 10x Genomics doublet rate table  
# - The percentage of doublets is based on the number of loaded and recovered cells  

# Best Practices for Using DoubletFinder  
# - Do not aggregate or integrate data before running DoubletFinder  
# - DoubletFinder is a quality control (QC) step and should be applied early  
# - Do not run DoubletFinder on merged datasets, as different samples may contain varying proportions of cells  
# - Applying DoubletFinder to large merged datasets can also crash R  
# - Run DoubletFinder separately on each sample after performing initial QC  


# install.packages("remotes")
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)


#Create count matrix 
cts <- ReadMtx(mtx= 'doublets/raw_feature_bc_matrix/matrix.mtx.gz', features = 'doublets/raw_feature_bc_matrix/features.tsv.gz', cells= 'doublets/raw_feature_bc_matrix/barcodes.tsv.gz' )

#Create a Seurat Object
pbmc.seurat <- CreateSeuratObject(counts = cts, project = "doublets")

#QC

pbmc.seurat$mitoPercent <-  PercentageFeatureSet(pbmc.seurat, pattern = "^MT-")

pbmc.seurat <- subset(pbmc.seurat, subset = nCount_RNA > 800 & nFeature_RNA> 500 & mitoPercent < 10)


View(pbmc.seurat@meta.data)



#pre-process standard workflow

pbmc.seurat <- NormalizeData(object = pbmc.seurat)
pbmc.seurat <- FindVariableFeatures(pbmc.seurat)
pbmc.seurat <- ScaleData(pbmc.seurat)
pbmc.seurat <- RunPCA(object = pbmc.seurat)

ElbowPlot(object = pbmc.seurat)


pbmc.seurat <- FindNeighbors(object = pbmc.seurat, dims=1:20)
pbmc.seurat <- FindClusters(object = pbmc.seurat)

pbmc.seurat <- RunUMAP(object = pbmc.seurat, dims= 1:20)

DimPlot(pbmc.seurat, reduction = 'umap', label=TRUE)

# pN default of .25 is standard
#exp comes from 10x genomic table

#pK value (no ground-truth)

sweep.res.list.nscls <- paramSweep(pbmc.seurat, PCs = 1:20, sct= FALSE) #merge artificial doublets in different proprtions. Pre-process data and calculate proportion of artifical nearest neghbors. Provides a list of the proprtion of artificla nearest nehgbors

swwep.stat.nscls <- summarizeSweep(sweep.list = sweep.res.list.nscls, GT=FALSE) #summary for each combination of pN and pK

bcmvn_nscls <- find.pK(swwep.stat.nscls) #find optimal pK value. Highest meanVaraince. The. highest value of that corresponds to the optimal pK value. Generates a plot and creates a table.

#pK value corresponding to max BC metric which is .21
ggplot(bcmvn_nscls , aes(pK, BCmetric, group=1)) + 
  geom_point() +
  geom_line()

#store optimal pK value.
pK<- bcmvn_nscls  %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)

#since it is stored as an array we extract it
pK <- as.numeric(as.character(pK[1]))


#Homotypic Doublet Proportion Estimate

annotations <- pbmc.seurat@meta.data$seurat_clusters


homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filtered <- doubletFinder(pbmc.seurat, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_28_691")


# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_28_691)
