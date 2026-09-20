# RNA-Seq Analysis of Paclitaxel-Resistant MCF-7 Breast Cancer Cells

## Overview
Paclitaxel resistance is a major challenge in breast cancer treatment.
This project investigates transcriptomic alterations associated with paclitaxel resistance in MCF-7 breast cancer cells using RNA-Seq analysis.

## Objectives
- Identify differentially expressed genes (DEGs)
- Perform GO enrichment analysis
- Conduct Gene Set Enrichment Analysis (GSEA)
- Identify hub genes using Cytoscape and MCODE
- Explore potential therapeutic targets

## Dataset
Source: Gene Expression Omnibus (GEO)
Dataset: GSE113685
Samples taken: MCF-7 Paclitaxel-resistant: SRR32240819, SRR32240820, SRR32240821
               MCF-7 Parental(Paclitaxel-sensitive): SRR32240825, SRR32240826, SRR32240827

## Workflow
RNA-Seq Data
↓
FastQC
↓
Kallisto Quantification
↓
Differential Expression Analysis
↓
Volcano Plot Generation
↓
GO Enrichment Analysis
↓
GSEA
↓
Hub Gene Identification
↓
Target Discovery

## Tools Used
- R
- Kallisto
- clusterProfiler
- Cytoscape
- MCODE

## Key Skills Demonstrated
- RNA-Seq Analysis
- Differential Expression Analysis
- Functional Enrichment
- Network Biology
- Bioinformatics Data Interpretation
