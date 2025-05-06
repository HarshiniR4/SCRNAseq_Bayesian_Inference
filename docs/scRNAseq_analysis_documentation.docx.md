# Single-cell RNA‑seq Analysis (GSE98969): Technical Documentation

***Notebook source: scRNAseq\_analysis.ipynb***

# **1 Objective**

Provide a complete, reproducible workflow for analysing single‑cell RNA‑seq data from the GEO dataset GSE98969, with emphasis on microglial transcriptional states in Alzheimer’s disease (AD). Specific goals include:  
• Extract raw 10x matrices and construct a gene‑by‑cell count matrix.  
• Perform rigorous quality control and normalisation.  
• Identify discrete and transitional microglial states via unsupervised clustering.  
• Discover marker genes and enriched pathways for each state.  
• Infer activation trajectories and quantify cell‑type frequency shifts between AD and control.

# **2 Procedure**

## **Data Ingestion**

Untar GEO archive (GSE98969\_RAW.tar), decompress sample‑wise 10x matrices, concatenate into a single AnnData object and append sample metadata.

## **Quality Control**

Filter low‑quality cells (n\_genes \< 200, pct\_mito \> 5 %) and low‑information genes, followed by library‑size normalisation and log transformation.

![image1](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture1.png)

## **Feature Selection & PCA**

Select \~1 800 highly variable genes (HVGs) and compute 50 principal components to capture dominant sources of variance.

![image2](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture2.png)

## **Neighbour Graph & Clustering**

Construct a k‑nearest‑neighbour graph (k = 15) on top PCs and apply Leiden community detection (resolution = 0.5) to define cell clusters.

![image3](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture3.png)

The UMAP plot you generated shows the clustering of your RNA data after dimensionality reduction and clustering. Here’s what the plot tells us:

**Key Observations:** 1\. **Clusters (Colors)**: \- The plot displays distinct clusters, each represented by a different color. \- Each color corresponds to a different Leiden cluster ID, as indicated by the color bar.

2\. **Separation of Clusters**:

• Some clusters are well-separated, indicating distinct groups of cells with different gene  expression profiles.

• Overlapping clusters suggest potential biological similarities or gradual transitions between cell states.

3\. **Density and Spread**:

• The density of points within a cluster varies. Dense regions could indicate cell types with a high degree of similarity.

• More spread-out clusters might represent heterogeneous populations or transitional states.

4\. **Cluster Sizes**:

• Some clusters are larger than others, suggesting that certain cell types or states are more prevalent in your dataset.

**Biological Implications:** \- **Distinct Cell Types or States**: \- The separation into clusters likely corresponds to different cell types, subtypes, or distinct cellular states.

• **Heterogeneous Populations**:

**–** Overlapping clusters or transitional regions may represent cells transitioning from one state to another or mixed populations.

• **Potential Outliers**:

**–** Any small, isolated clusters could be rare cell types or potentially outliers depending with any of these follow-up analyses?

**Identifying Marker Genes for Each Cluster** Marker genes help determine the biological identity of each cluster by highlighting genes that are highly expressed in a specific cluster compared to others.

Cell-Type Marker Genes:

The figure provides a list of validated marker genes for neurons, astrocytes, and microglia: Neurons (NeuN+): Markers: Rbfox3, Grin1, Cx3cl1, Nos1 Astrocytes (GFAP+): Markers: Gfap, Aldh1l1, Aqp4, Slc1a2 Microglia (CD11b+): Markers: Itgam (CD11b), Aif1, Cx3cr1, Ptprc (CD45) These genes can be directly used to annotate clusters in your dataset. For example: If a cluster shows high expression of Rbfox3 or Grin1, it likely represents neurons. If a cluster expresses Gfap or Aldh1l1, it may represent astrocytes. If Itgam or Aif1 is highly expressed, it could correspond to microglia.

## **Visualisation**

Embed cells in 2‑D with UMAP for qualitative inspection of cluster relationships.

![image4](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture4.png)

The UMAP plot shows distinct clusters of cells, indicating diverse cell types or states. Clusters’ proximity and separation suggest biological relationships or transitions between these cell populations

## **Marker Gene Discovery**

Rank genes per cluster using Wilcoxon tests (logFC \> 0.25), generating cluster‑specific signatures.

![image5](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture5.png)

## **Pathway Enrichment**

Run GO/KEGG over‑representation analysis on marker sets to contextualise functional programmes.

## **Trajectory Inference**

Build a PAGA connectivity graph and compute diffusion pseudotime to model activation paths from homeostatic to disease‑associated microglia.

![image6](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture6.png)

**Pseudotime Projection on UMAP**

• **Color Gradient**: Represents pseudotime, where:

*  Dark purple areas indicate earlier pseudotime (starting cells in the trajectory).

**–** 	Yellow/green areas indicate later pseudotime (more differentiated or terminal cells).

**Key Insights**: 

1\. **Developmental Trajectories**: \- Cells in the dark purple regions likely represent early progenitors or naive cells. \- Cells transitioning to yellow/green regions show progression toward differentiated states, such as monocytes or dendritic cells.

2\. **Cluster Overlap**:

• The smooth gradient from one region to another implies continuous cellular transitions between certain cell types.

• For example, some cells transition from progenitor-like states (low pseudotime) to monocytes or dendritic cells (high pseudotime).

3\. **Dynamic Processes**:

• Pseudotime reveals dynamic biological processes like ie assistance with further analysis or validation steps?

**Biological Implications:** 

1\. **Immune Cell Development**: \- Early pseudotime clusters could represent stem or progenitor cells, while later stages show differentiated immune cells such as monocytes or dendritic cells.

2\. **Cell Fate Decisions**:

• Pseudotime gradients suggest how cells transition from one type to another, providing insights into developmental pathways or immune responses.

3\. **Potential for Rare Cell Types**:

• Investigating the “NA” cluster and pseudotime extremes could help identify rare or novel cell types.

## **Differential Expression & Composition**

Aggregate counts into pseudobulk profiles per sample×cell‑type, fit negative‑binomial GLMs for AD vs control contrasts, and model compositional shifts with scCODA.

# **3 Key Outputs**

* adata\_qc.h5ad – Filtered single‑modality AnnData (25 032 cells × 24 379 genes).  
* mdata\_processed.h5mu – Multimodal MuData with joint embeddings and cluster annotations.  
* figures/\*.png – UMAPs, QC plots, PAGA graph, volcano and enrichment plots.  
* de\_results/\*.csv – Differential‑expression tables (gene, logFC, p‑value, FDR) per cell type.  
* pseudobulk/counts.csv – Gene counts aggregated by sample and cell type.

# **4 Interpretations**

Unsupervised clustering resolved distinct microglial states, including a disease‑associated microglia (DAM/ARM) cluster enriched for lipid‑handling and phagocytosis genes (Apoe, Trem2, Cst7), and an interferon‑responsive cluster characterised by Ifit3 and Irf7. Trajectory analysis suggests a continuum from homeostatic microglia towards the DAM phenotype, supporting progressive activation in the AD brain. Pathway enrichment highlights lipid catabolic processes and cytokine signalling as dominant themes in activated states, while composition modelling indicates an increased proportion of DAM cells in 5XFAD samples relative to wild type.

## **4.1 Representative Cluster Signatures**

| Cluster ID | Top Markers | Putative Identity | Functional Notes |
| :---- | :---- | :---- | :---- |
| 0 | P2ry12, Tmem119 | Homeostatic microglia | Baseline surveillance |
| 3 | Apoe, Trem2, Cst7 | DAM/ARM | Phagocytosis, lipid metabolism in AD |
| 5 | Ifit3, Ifitm3, Irf7 | Interferon‑responsive | Type‑I IFN  |

# **5 Limitations**

Batch correction is not applied; integrating additional donors or modalities may require Harmony or scVI. Pseudobulk differential expression assumes ≥3 replicates per condition. ADT and ATAC modalities are stubbed and await data insertion.

![image7](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture7.png)

![image8](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture8.png)

![image9](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture9.png)

![image10](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture10.png)

![image11](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture11.png)

![image12](https://github.com/HarshiniR4/SCRNAseq_Bayesian_Inference/blob/main/images/Picture12.png)

This volcano plot represents the results of a differential expression analysis for T\_Cells. Here’s what it implies:

1\. **Axes**:

• The x-axis shows the **log fold change (logFC)** in gene expression between conditions.

**–** Positive values indicate upregulation in T\_Cells.

**–** Negative values indicate downregulation.

• The y-axis shows **\-log10(FDR-adjusted p-values)**.

**–** Higher values indicate greater statistical significance.

2\. **Thresholds**:

• The vertical blue dashed lines represent the **log fold change threshold** (±1.5).

**–** Genes outside these lines have substantial changes in expression.

• The horizontal green dashed line represents the **FDR threshold (0.01)**.

**–** Genes above this line are statistically significant.

3\. **Points**:

• **Red points**: Genes that meet both criteria for differential expression (significant FDR and logFC \> 1.5 or \< \-1.5). These are the **significantly differentially expressed genes** (DEGs).

• **Gray points**: Genes that do not meet the criteria. These are either not statistically significant or do not show a strong enough fold change.

**Implications:** 

1\. **Upregulated Genes**: \- Genes with high positive logFC (on the right of the plot) and significant FDR are upregulated in T\_Cells under the condition being tested.

2\. **Downregulated Genes**:

• Genes with high negative logFC (on the left of the plot) and significant FDR are downregulated in T\_Cells.

3\. **Biological Insight**:

• The plot highlights genes whose expression levels are strongly and significantly altered in T\_Cells compared to other conditions or groups.

• These genes may be of biological interest for further investigation, as they could be

critical markers, pathways, or regulators involved in the processbased on specific research

goals.

