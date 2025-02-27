## This repo is used to centralize scRNA Seq Data analysis.

`Preoprocessing` - Retrieve h5 file that contains feature expression data, convert to seurat, perform normalizatiom, filtering, sclaing, PCA, clustering, and UMAP. Visualize umap by genes to see differential expression by clusters.

`Integration` - Check and Corrects batch effects by identifying shared features across patients, aligning datasets using canonical correlation analysis (CCA), and generating an integrated Seurat object for unbiased clustering and visualization.

`Doublets` -  detect and remove doublets in single-cell RNA sequencing data. Doublets occur when two cells are captured within the same droplet, leading to mixed gene expression profiles.

`Harmony` - use Harmony to correct batch effects and integrate data from different conditions.

`Cell-Markers` - find markers using both FindConservedMarkers(across clusters but persistnet per condition), and FindMarkers (between conditions) in order to find differential expressed markers and to annotate cells. 
