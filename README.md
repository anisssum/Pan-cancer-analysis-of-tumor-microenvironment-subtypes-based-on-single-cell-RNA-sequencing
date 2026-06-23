# Pan-cancer-analysis-of-tumor-microenvironment-subtypes-based-on-single-cell-RNA-sequencing
## Overview

This repository contains the computational workflow used to identify tumor microenvironment (TME) subtypes across multiple cancer types using single-cell RNA sequencing (scRNA-seq) data. The pipeline integrates datasets from several public cancer atlases, harmonizes cell-type annotations, quantifies TME composition, performs consensus clustering, derives gene-expression signatures associated with T-cell infiltration, and investigates cell-cell communication programs.

The study includes both a pan-cancer analysis and an independent colorectal cancer (CRC) validation cohort.

# Pan-cancer Analysis

Folder: Pan-cancer-analysis-of-tumor-microenvironment-subtypes-based-on-single-cell-RNA-sequencing/Pan-cancer/

---

## Dataset-Specific Processing

### TISCH

#### Pseudobulk_tisch.py

Generates malignant-cell pseudobulk profiles from TISCH datasets.

#### Liana_tisch.py

Performs ligand-receptor analysis within TISCH cohorts.

---

### Korean Cohort

#### Pseudobulk_korean.py

Generates malignant-cell pseudobulk profiles from Korean cancer cohorts.

#### Liana_koren.py

Performs ligand-receptor analysis for Korean datasets.

---

### 3CA Atlas

#### Pseudobulk_3ca.py

Generates pseudobulk expression matrices from the 3CA atlas.

#### Liana_3ca.py

Performs cell-cell communication analysis using LIANA.

---

## Liana_combined.py

Integrates ligand-receptor interaction results from multiple datasets.

### Main tasks

* Dataset harmonization
* Aggregation of communication scores
* Cross-cohort comparison

### Output

* Combined interaction tables
* Consensus communication networks

---

## Pancancer_non_malignant_celltypist_composition.py

Processes harmonized scRNA-seq datasets and calculates sample-level cellular composition profiles.

### Main tasks

* Standardization of cell-type annotations
* Grouping related immune and stromal populations
* Quantification of cell counts per sample
* Generation of cell-type abundance matrices

### Output

* Sample-by-cell-type count matrix
* Cell-type composition profiles used for clustering

---

## Pseudobulk.R

Generates pseudobulk expression profiles from malignant cells.

### Main tasks

* Aggregation of single-cell counts
* Sample-level expression matrix construction
* Preparation of expression data for signature development

### Output

* Pseudobulk expression matrix

---

## ConsensusClustering_CLR_nonzero.R

Identifies recurrent TME subtypes based on cell-type composition.

### Main tasks

* Conversion of cell counts to proportions
* Zero replacement using pseudocounts
* Centered log-ratio (CLR) transformation
* Euclidean distance calculation in CLR space
* Consensus clustering using PAM

### Output

* Consensus matrices
* Cluster stability metrics
* Optimal cluster number estimation

---

## GeneStats_malignant.py

Computes gene-level statistics across malignant-cell pseudobulk samples.

### Main tasks

* Expression filtering
* Gene ranking
* Differential expression metrics
* Summary statistics generation

### Output

* Gene-level statistical summaries

---

## Cluster_analysis.R

Performs downstream characterization of TME clusters.

### Main tasks

* Visualization of consensus clustering results
* Heatmap generation
* UMAP embedding
* Cluster assignment
* Comparison of cellular compositions between clusters

### Output

* Cluster labels
* Heatmaps
* UMAP projections
* Statistical comparisons
* 
---

## ElasticNet.R

Builds a pan-cancer gene-expression signature associated with CD8⁺ T-cell infiltration.

### Main tasks

* CLR transformation of composition data
* Construction of CD8 infiltration phenotype
* Elastic Net regression
* Gene selection and coefficient estimation
* Functional enrichment analysis

### Output

* Pan-cancer signature genes
* Gene coefficients
* Signature scores
* Enrichment results

---

## Liana_signature.R

Infers ligand-receptor interactions associated with the pan-cancer signature.

### Main tasks

* Signature stratification
* Cell-cell communication analysis
* Pathway-level interpretation

### Output

* Signature-associated communication networks

---

# CRC Analysis

Folder: Pan-cancer-analysis-of-tumor-microenvironment-subtypes-based-on-single-cell-RNA-sequencing/CRC/

## CRC_analysis.py

Processes colorectal cancer datasets and constructs CRC-specific composition matrices.

### Main tasks

* Cell-type harmonization
* Sample filtering
* Composition calculation
* Preparation for clustering and signature modeling

### Output

* CRC cell-type composition matrix

---

## ConsensusClustering_CLR_nonzero.R

Performs compositional clustering of CRC tumors.

### Main tasks

* CLR transformation
* Consensus clustering
* Cluster stability evaluation

### Output

* CRC TME subtypes

---

## Cluster_analysis.R

Characterizes CRC TME clusters.

### Main tasks

* Heatmap generation
* UMAP visualization
* Cluster assignment
* Cluster-level biological interpretation

### Output

* CRC subtype labels
* Visualization figures

---

## ElasticNet.R

Builds the CRC-specific CD8⁺ T-cell infiltration signature.

### Main tasks

* Pseudobulk preprocessing
* Elastic Net modeling
* Signature derivation
* Functional enrichment analysis

### Output

* CRC signature genes
* Gene coefficients
* Functional enrichment results

---

# Methodological Workflow

1. Collect and harmonize scRNA-seq datasets.
2. Standardize cell-type annotations.
3. Quantify TME composition per sample.
4. Apply CLR transformation to compositional data.
5. Identify TME subtypes via consensus clustering.
6. Generate malignant-cell pseudobulk profiles.
7. Train Elastic Net models to derive CD8-associated signatures.
8. Perform pathway enrichment analysis.
9. Investigate cell-cell communication programs using LIANA.
10. Validate findings in an independent CRC cohort.

---

# Software Requirements

## Python

* scanpy
* anndata
* pandas
* numpy
* scipy
* celltypist
* liana

## R

* ConsensusClustering
* compositions
* zCompositions
* glmnet
* limma
* enrichR
* tidyverse
* pheatmap
* uwot
* ggplot2
* patchwork

---
