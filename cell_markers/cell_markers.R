#FindMarkers
# Purpose: Identifies differentially expressed genes between two specific groups or clusters.
#Use Case: When you want to compare a particular cluster against another specific cluster or condition.

#Example: Comparing Cluster A vs. Cluster B to find genes that are upregulated in Cluster A compared to Cluster B.


# FindAllMarkers
#Purpose: Identifies differentially expressed genes for each cluster compared to all other clusters.

#Use Case: When you want to find marker genes for each cluster without specifying comparisons.

#Example: Finding genes that are upregulated in Cluster A compared to all other clusters, and similarly for other clusters.



#FindConservedMarkers
#Purpose: Identifies genes that are differentially expressed within each condition and are conserved across multiple conditions.

#Use Case: When you have data from multiple conditions (e.g., control vs. treatment) and want to find genes that are consistently differentially expressed in a specific cluster across these conditions.
#Example: Finding genes that are upregulated in Cluster A in both control and treatment conditions.

library(SeuratData)
library(tidyverse)
library(ggplot2)

ifnb <- readRDS("./ifnb_harmony.rds")

str(ifnb)
View(ifnb@meta.data)


#visualize data
clusters <- DimPlot(ifnb, reduction = 'umap', group= 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb, reduction = 'umap', group= 'stim')

clusters | condition #We can see that the data is integrated and separated by cluster. We want to identify what cells from what cluster (even though we have it).

#Identify cell types that form each cluster

#FindConservedmarkers ---> data from two conditons and cells in both


DefaultAssay(ifnb) #returns RNA bc we used harmony. No other assay. We need it to be RNA for FindAllMarkers
#DefaultAssay(ifnb) <- 'RNA' #in case not rna

#FindAllMakers - data from one condition. We can test it out
FindAllMarkers(ifnb, logfc.threshold = 0.25, #default is .25
               min.pct=.1 #only detect genes that are at 50% freq in other of the two clusters we compare. Default is .1
                only.pos=TRUE #only return markers that are upregulated
                test.use = 'DESeq2', slot = 'counts' ) #many different ones

#We will use FindConservedMarkers. Indentify cluster 3. No Ident.2 bc we are not comparing 2 clusters. Group by condition
markers_cluster_3 <- FindConservedMarkers(ifnb, ident.1 = 3, grouping.var = 'stim') #separate cells by the condition. for each condition it compared cluster 3 to all other conditions

#First gene has a p-value of 0 and fold change of 4.099783 Detected in 97% of cluster 3 and 20% of all other genes
head(markers_cluster_3)

FeaturePlot(ifnb, features =c ('FCGR3A'), min.cutoff = 'q10') #We see this gene is highly expressed in cluster 3 . Using q10 means that for cells that cells with this gene expressed under q10 will be gray. The gene expression values about this will have cells with color


#Min cuttoff
seq(1,5)

#Find the quantile 50
SetQuantile('q50', seq(1,5)) #median

SetQuantile('q10', seq(1,5))



#rename cluster 3 ident to CD16 Monocyte
Idents(ifnb) #cell idemtificaions are the cluster number. Replacce to name of cell
ifnb = RenameIdents(ifnb, '3' = 'CD16 Monocytes')


DimPlot(ifnb, reduction = 'umap', label=TRUE)


Idents(ifnb) <- ifnb@meta.data$seurat_annotations #use annotations for indets


# now we want to look at gene expression levles before and after treatment. Split cluster 3 into control and stimulated

#create cell type andc ondition column in data

ifnb$celltype.md <- paste0(ifnb$seurat_annotations, '-', ifnb$stim)

Idents(ifnb) <- ifnb$celltype.md

DimPlot(ifnb, reduction = 'umap', label=TRUE) #each cluster has multiple labels bc cells belong to two conditions

#plotting conserved features vs DE fatures between conditions. FindMarkers to compare CD 16 cells of the two conditions

interferon_response <- FindMarkers(ifnb, ident.1 = 'CD16 Mono-STIM', ident.2 = 'CD16 Mono-CTRL')

head(interferon_response) #found genes that are differnetly expressed across the same cell but differnet condition. EX: IFIT1 is expressed more in stim than control


#Plot expression of markers of conserved vs DE fearures between conditinos

head(markers_cluster_3)

# two makrers from conserved, one from findmarkers

FeaturePlot(ifnb, features = c('FCGR3A','VMO1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10') #first two markers are expressed in both control and stim. Last marker expressed in all cells in stim but not at all in control 
