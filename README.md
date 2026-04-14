# Glioma-Cell-Type-Specific-Causal-Genes

This repository contains code supporting "Single-cell multi-omic integration analysis prioritizes druggable genes and reveals cell-type-specific causal effects in glioblastomagenesis" by Huang et al. (2025). You can access the preprint at medRxiv.
cis-QTL reference panel weights generated for this study used in TWAS and PWAS are available at https://figshare.com/articles/dataset/_b_TWAS_PWAS_gene_expression_prediction_models_b_/28080035. 


MungeSumstats	v1.10.1	

Seurat 	v4.4.0	

inferCNV 	 v1.18.1	

PoPS	v0.2	

MAGMA 	v1.10	

Mapgen	v0.5.8	

FUSION	-	

ACAT-O	v0.91	

GCTA-GREML	v1.94.2	

SMR	v1.3.1	

coloc	v5.2.3	

STRING	v12.0	

TissueEnrich	v1.24.1	

clusterProfiler	v4.10.0	

irGSEA	v2.1.5	

limma 	v3.58.1	

edgeR 	v4.0.3	

CellChat 	v2.1.2	

CellPhoneDB	v5.0.1	

scPagwas	v1.3.0	

LDSC	v1.0.1	

CytoTRACE2	v1.0.0	

VECTOR	v0.0.4	

PAGA	-	

monocle3	v1.4.18	

cisMRcML	v0.0.09	

EWCE	v1.10.2	

CTFM	-	

CARMA	v1.0	

susieR	v0.12.35	

ldscr	v0.1.0	

TwoSampleMR	v0.5.8	

RadialMR	v1.1	

MRPRESSO	v1.0	

MRcML	v0.0.0.9	

HDL	v1.4.0	

LCV	-	

ViSNP	v0.1.0	


## R package scaffold for GWAS-post workflows

A new package scaffold is available under `GWASpost/`.

- Main implemented function: `run_coloc()`
- Purpose: run the full colocalization workflow directly from one GWAS summary file and a local directory of eQTL summary files.

See `GWASpost/README.md` for required columns and usage.
