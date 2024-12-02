setwd("G:/Database/CGGA")

pacman::p_load(vroom, dplyr, limma, edgeR, readr)
### Process CGGA data ----
A <- vroom("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt")
A <- as.data.frame(A)
rownames(A) <- A$gene_name
B <- A[,-1]
saveRDS(B, file = "CGGA_693_ALL.rds.gz")
clinical_data <- vroom("CGGA.mRNAseq_693_clinical.20200506.txt")
colnames(clinical_data)
GBM_samples <- clinical_data %>% 
  filter(Histology == "GBM") %>% 
  pull(CGGA_ID)
other_samples <- clinical_data %>% 
  filter(Histology != "GBM") %>% 
  pull(CGGA_ID)
A_GBM <- A %>% select(all_of(GBM_samples))
A_other <- A %>% select(all_of(other_samples))
saveRDS(A_GBM, file = "CGGA_693_GBM.rds.gz")
saveRDS(A_other, file = "CGGA_693_noGBM.rds.gz")
A_normal <- vroom("CGGA.normal_20.Read_Counts-genes.20230104.txt")
A_normal <- as.data.frame(A_normal)
rownames(A_normal) <- A_normal$gene_name
A_normal <- A_normal[,-1]
saveRDS(A_normal, file = "CGGA_20_NORMAL.rds.gz")
rm(list = ls())
gc()
B_tumour <- readRDS("CGGA_693_ALL.rds.gz")
B_normal <- readRDS("CGGA_20_NORMAL.rds.gz")
B_tumour <- B_tumour[rowMeans(B_tumour) > 1, , drop = FALSE]
B_normal <- B_normal[rowMeans(B_normal) > 1, , drop = FALSE]
B_tumour$rowname <- rownames(B_tumour)
B_normal$rowname <- rownames(B_normal)
all.data <- merge(B_tumour, B_normal, by = "rowname")
rownames(all.data) <- all.data$rowname
all.data$rowname <- NULL
saveRDS(all.data, file = 'CGGA_ALL_693_normal.rds.gz')
rm(list = ls())
gc()
dlbc.exp <- readRDS("CGGA_ALL_693_normal.rds.gz")
dlbc.exp <- dlbc.exp[rowMeans(dlbc.exp) > 0, ]
Sample_infor <- data.frame(sample_name = colnames(dlbc.exp))
Sample_infor$group <- ifelse(grepl("^N", Sample_infor$sample_name), "normal", "tumour")
###DEG Analyses----

dge <- DGEList(counts = dlbc.exp)
dge <- calcNormFactors(dge)
design <- model.matrix(~ factor(Sample_infor$group))
colnames(design) <- levels(factor(Sample_infor$group))
logCPM <- voom(dge, design, normalize = "quantile")
fit <- lmFit(logCPM, design)
fit <- eBayes(fit)
output <- topTable(fit, coef = 2, n = Inf)
output$GENE <- rownames(output)
A <- subset(output, abs(logFC) > 1 & adj.P.Val < 0.05)
write_tsv(output, file = "pan_CGGA_693_DEG_limma_all.txt")
write_tsv(A, file = "pan_CGGA_693_DEG_limma_significant.txt")
rm(list = ls())
gc()
