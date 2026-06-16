# GCLink: Genetic Causality and Shared Molecular Mechanisms Linking Complex Diseases
GCLink is an integrated analytical pipeline designed to uncover the genetic causal relationships and shared molecular mechanisms between complex diseases, such as anxiety disorder (ANX) and allergic diseases. By integrating large-scale epidemiological data, genome-wide association study (GWAS) summary statistics, quantitative trait loci (QTL), bulk and single cell RNA sequencing data, GCLink provides a comprehensive framework for dissecting disease comorbidity from genetic causality to molecular mechanisms.

## Overview
Understanding whether two complex diseases are causally related and identifying the molecular mechanisms underlying their comorbidity remain major challenges in human genetics.
GCLink addresses these challenges through a two-phase analytical framework:
### Phase 1: Observational Studies and Genetic Causal Inferences
(1) Observational studies; (2) Genetic correlation analyses; (3) Bidirectional Mendelian Randomization (MR).
### Phase 2: Molecular Mechanisms Underlying Genetic Correlations
(1) Significant mediators; (2) Shared genetic variants; (3) Tissue enrichment; (4) Shared functional genes; (5) Specific cell types; (6) Perturbed pathways.
The pipeline can be readily adapted to investigate the genetic relationships between any pair of complex diseases with available GWAS summary statistics.

## Data Preparation
### Required GWAS Summary Statistics Format
Before running GCLink, all GWAS summary statistics should be harmonized into the following format: chr (Chromosome number), pos (Genomic position), rsid (SNP identifier), A1 (Effect allele), A2 (Reference allele), beta (Effect size), se (Standard error), p (P-value), eaf (Effect allele frequency), N (Sample size).

## Data Resources
Preprocessed GWAS summary statistics in this study and reference panels used by GCLink are available through Google Drive: https://drive.google.com/drive/my-drive?dmr=1&ec=wgc-drive-%5Bmodule%5D-goto

## Citation
If you use GCLink in your research, please cite:
1. Rao S., Chen X., Ou O. Y., et al., A Positive Causal Effect of Shrimp Allergy on Major Depressive Disorder Mediated by Allergy- and Immune-Related Pathways in the East Asian Population. Nutrients, 2023. 16(1): 79.
2. Chen T. B., Jiang J. W., Guo H. Y., et al., Causal relationship between hepatic function indicators and thrombocytopenia risk in early-stage hepatitis B virus infection: evidence from clinical observational studies and mendelian randomization analyses. Front Immunol, 2025. 16: p. 1440317.





