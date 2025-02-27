# Batch Correction and Batch Effect

Single-cell RNA sequencing (scRNA-seq) allows us to study gene expression at the single-cell level in any tissue.
Workflow:

    Isolate the tissue of interest
    Extract single cells
    Capture RNA
    Sequence
    Obtain gene expression readouts and generate a gene expression matrix

Major Challenge: Batch Effect

Batch effects are systematic technical differences introduced during sample handling, library preparation, sequencing platforms, or operator variations. These effects can obscure true biological signals.
Example:

In a UMAP plot, different samples of the same cell type may cluster separately due to batch effects.
After correction, the cells should overlap because they belong to the same cell type.
Methods for Batch Effect Correction
1. Limma and ComBat

    Limma: Assumes batch effects are the same across all cells.
    ComBat: Uses empirical Bayes shrinkage, making it more flexible. However, it still assumes batch effects are uniform across all cells.
    Limitation: Both approaches are linear, meaning they may not work well if the data is non-linear.

2. Mutual Nearest Neighbors (MNN) – A Better Approach

    Identifies overlapping cell populations to correct batch effects.
    Finds pairs of similar cells across different batches.
    Removes batch effects only in overlapping subpopulations, making it more adaptable.

How MNN Works

    Identify cells based on their gene expression profiles. Each cell’s gene expression is represented as a vector.
    Find similar cells across batches that have similar expression vectors.
    Define Mutual Nearest Neighbor (MNN) pairs:
        Each cell finds its k nearest neighbors.
        If two cells from different batches identify each other as their nearest neighbor, they form an MNN pair.

This approach is effective because it aligns cells based on shared biological signals rather than assuming uniform batch effects.



# Markers
Markers: differentially expressed genes that let us identify cells
FindAllMarkers - looks at all population to find markers 
FindMarkers - look at 2 populations to find markers
Wilcox - default statistical test. 
min cells per group : how many cells per a cluster size.
min pct: min % of cells expressing that gene in either groups


# Imputation in Single-Cell RNA Sequencing (scRNA-seq)

In single-cell RNA sequencing (scRNA-seq), a common challenge is the presence of zero values in the data. These zeros can arise due to a variety of reasons:

    Lack of gene expression: The gene may not be expressed in that particular cell.
    Technical issues: Such as sequencing errors or inefficiencies in data capture.
    Sampling errors: Due to limitations in detecting lowly expressed genes.

Zero Inflation and Dropout:

    Dropout: A phenomenon where genes are truly expressed but are not detected due to technical limitations, leading to zeros in the data.
    Zero Inflation: This occurs when zeros dominate in the data, possibly because we start with zero expression, and any amplification during sequencing will just produce more zeros.

Using UMIs (Unique Molecular Identifiers) helps reduce zero inflation because it tags each RNA molecule before amplification. This ensures that we don’t overcount or misrepresent the presence of RNA molecules due to PCR amplification bias, thus minimizing technical zeros.
Imputation Methods:

Imputation techniques are applied to distinguish between technical zeros (arising from dropout or sampling issues) and biological zeros (genuine absence of gene expression). Imputation aims to predict or fill in the missing values based on patterns in the data, essentially reconstructing the data to more accurately reflect the biological process.
Types of Imputation:

    Model-based imputation: Utilizes mathematical models to infer missing values. This helps identify technical zeros rather than biological zeros.

    Data Smoothing: This technique adjusts the expression values for each cell based on the expression values of similar or neighboring cells. Smoothing helps reduce noise and can provide a more accurate picture of gene expression.

    Data Reconstruction: This involves leveraging external datasets or pretrained models to fill in missing or zero values, based on knowledge from other similar datasets or biological contexts.


Single R: reference based annotations

Seurat reference assembly integrstion