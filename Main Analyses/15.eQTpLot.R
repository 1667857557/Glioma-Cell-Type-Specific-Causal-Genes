pacman::p_load("biomaRt", "vroom", "dplyr", "eQTpLot", "GenomicRanges", "ggnewscale", "ggplot2", 
               "ggplotify", "ggpubr", "gridExtra", "Gviz", "LDheatmap", "patchwork", 
               "ieugwasr", "plinkbinr")
gc()
setwd("G:/A/target/GWAS")
gwas_fn<-vroom("GBM_GLIOMA.txt.gz",col_select=c("CHR","POS","SNP","P","BETA","SE"))
gwas_fn$PHE<-"GBM"
colnames(gwas_fn)[colnames(gwas_fn) == "POS"] <- "BP"
eqtl_fn<-vroom("G:/ENSG00000146648.txt",col_select=c("SNP","GENE","P","BETA","N"))
eqtl_fn$Tissue<-"Astrocytes_eQTL(Fujita et al.)"
gc()
merged = merge(gwas_fn, eqtl_fn, by = "SNP", suffixes = c("1", "2"), all = FALSE)
colnames(eqtl_fn)<-c("SNP.Id","Gene.Symbol","P.Value","NES","N","Tissue")
eqtl_fn$Gene.Symbol<-"EGFR"
merged <- merged[!duplicated(merged$SNP), ]
merged<-dplyr::select(merged,c(SNP,CHR,BP))

LD<- ld_matrix_local(merged$SNP,
                     bfile = "G:/Database/1000G/1kg.v3/EUR",
                     plink_bin = "G:/plink_Windows.exe", 
                     with_alleles = FALSE)
SNP_pairs <- combn(dimnames(LD)[[1]], 2)
R2_values <- as.vector(LD)[lower.tri(LD)]

LD_data <- data.frame(
  SNP_A = SNP_pairs[1, ],
  SNP_B = SNP_pairs[2, ],
  R2 = R2_values
)
A<-vroom("G:/A/target/GWAS/GBM_GLIOMA.txt.gz",col_select=c("SNP","POS"))
B<-A
colnames(A)[colnames(A) == "SNP"] <- "SNP_A"
colnames(A)[colnames(A) == "POS"] <- "BP_A"

colnames(B)[colnames(B) == "SNP"] <- "SNP_B"
colnames(B)[colnames(B) == "POS"] <- "BP_B"
LD_data<-left_join(LD_data,A,by = "SNP_A",relationship = "many-to-many")
LD_data<-left_join(LD_data,B,by = "SNP_B",relationship = "many-to-many")

LD_data$BP_A <- as.integer(as.character(LD_data$BP_A))
LD_data$BP_B <- as.integer(as.character(LD_data$BP_B))
setwd("G:/")
eQTpLot(GWAS.df = gwas_fn, eQTL.df = eqtl_fn, gene = "EGFR", res = 600,tissue =  c("Astrocytes_eQTL(Fujita et al.)"),
        gbuild = "hg19", trait = "GBM",LD.df =LD_data, R2min = 0.25, LDmin = 100, congruence = TRUE)
