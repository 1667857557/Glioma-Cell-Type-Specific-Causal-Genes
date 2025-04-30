pacman::p_load("vroom","sva","stringr","RColorBrewer","FactoMineR","factoextra","limma","rjson","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse")
setwd("G:/Database/TCGA")
B_tumour <- readRDS("tcga_GBM_count.rds.gz")
B_normal <- readRDS("tcga_norm_count.gz")
B_tumour$rowname <- rownames(B_tumour)
B_normal$rowname <- rownames(B_normal)
setwd("G:/Database/TCGA/GTEX")
gtex<-readRDS("gtex_GBM_norm.gz")
gtex$rowname <- rownames(gtex)
all.data <- merge(B_tumour, B_normal, by = "rowname")
all.data <- merge(all.data, gtex, by = "rowname")
rownames(all.data) <- all.data$rowname
all.data$rowname <- NULL
all.data[,] <- lapply(all.data[,], as.numeric)
all.data <- all.data[rowMeans(all.data) > 1, , drop = FALSE]
#saveRDS(all.data, file = 'TCGA_GBM_normal_GTEX.rds.gz')
rm(list = ls())
gc()
dlbc.exp <- readRDS("TCGA_GBM_normal_GTEX.rds.gz")
Sample_infor <- data.frame(sample_name = colnames(dlbc.exp))
Sample_infor$group <- ifelse(grepl("TCGA", Sample_infor$sample_name), "tumour",
                             ifelse(grepl("GTEX", Sample_infor$sample_name) | grepl("^data", Sample_infor$sample_name), "normal", "normal"))

Sample_infor$Sample_origination <- ifelse(grepl("GTEX", Sample_infor$sample_name), "GTEX",
                                          ifelse(grepl("TCGA", Sample_infor$sample_name) | grepl("^data", Sample_infor$sample_name), "TCGA", "TCGA"))

str(dlbc.exp[1:3, 1:3])
PCA.plot = function(dat,col){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups")
}
gc()
library(sva)
library(edgeR)
expr_count_combat <- ComBat_seq(counts = as.matrix(dlbc.exp), 
                                batch = Sample_infor$Sample_origination,
                                group = Sample_infor$group)
dge <- DGEList(counts = expr_count_combat)
dge <- calcNormFactors(dge, method = "TMM")
expr_count_combat_cpm <-  edgeR::cpm(dge, log = TRUE, prior.count = 1)
p4 <- PCA.plot(expr_count_combat_cpm,paste0(Sample_infor$Sample_origination,"-",Sample_infor$group))
pdf("PCA_ComBat_seq_TCGA_GTEx_GBM.pdf", width = 8, height = 6)
p4
dev.off()
saveRDS(expr_count_combat_cpm ,file = "ComBat_seq_TCGA_GTEx_GBM.rds")
library(limma)
library(edgeR)
group <- factor(Sample_infor$group, levels = c("normal", "tumour"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_count_combat_cpm, design)
contrast_matrix <- makeContrasts(tumour - normal, levels = design)
fit_contrast <- contrasts.fit(fit, contrast_matrix)
fit_ebayes <- eBayes(fit_contrast)

de_genes <- topTable(fit_ebayes, 
                     number = Inf, 
                     adjust.method = "BH",
                     p.value = 0.05,
                     lfc = 1)
de_genes$GENE<-rownames(de_genes)
write.csv(de_genes, "DE_TCGA_GTEx_GBM.csv")
