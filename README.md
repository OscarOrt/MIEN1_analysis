# MIEN1_analysis

Analysis code for the MIEN1 promoter ablation provides new evidence for colorectal cancer genome editing-based therapeutics project

This repository is intended to provide the code used for the bioinformatic analysis for the project MIEN1 promoter ablation provides new evidence for colorectal cancer genome editing-based therapeutics.

The code is implemented in R (v3.4.4) and the following libraries:
- edgeR v3.20.9
- org.Hs.eg.db v3.5.0
- pheatmap v1.0.10
- ggplot2 3.0.0
- dendsort v0.3.3
- RColorBrewer v1.1-2

#Components of the analysis

- Differential expression
This code uses edgeR was to identify over and under expressed genes.

- Gene ontology analysis and Plot gene ontologies
This code identifies gene ontology terms enriched in the differentially expressed genes and plot the results. It uses the function goana()

- Heatmap
This code uses the package pheatmap to plot selected differentially expressed genes and the gene ontology terms they belong.

- Molecular signature MSigDB
This code uses The Molecular Signatures Database (MSigDB) to identify molecular signatures in the genes differentially expressed.

