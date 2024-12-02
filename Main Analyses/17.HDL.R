pacman::p_load("vroom","data.table","readr","HDL")
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"# 修改为d盘路径
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
set.seed(1000)
### HDL-----
LD.path <- "G:/Database/UKB_imputed_SVD_eigen99_extraction"
gwas1.df <- vroom("G:/A/target/GWAS/GBM_GLIOMA.txt.gz")
colnames(gwas1.df)[colnames(gwas1.df) == "effect_allele"] <- "A1"
colnames(gwas1.df)[colnames(gwas1.df) == "other_allele"] <- "A2"
gwas1.df<-dplyr::select(gwas1.df,c(SNP,A1,A2,N,Z))
colnames(gwas1.df)<-c("SNP","A1","A2","N","Z")
# gwas1.df <- gwas1.df[!duplicated(gwas1.df$SNP), ]
gc()

gwas1.df<-as.data.frame(gwas1.df)
gwas1.df$N<-as.integer(gwas1.df$N)
results_df<-fread("G:/A/phenotypes/HDL_GBM.txt")

results_df <- data.frame()

setwd("G:/Database/phenotyes/buildGRCh37")
file_list <- list.files(pattern = "\\.gz$")
file_list <- file_list[!file_list %in% results_df$phenotypes]

for (i in seq_along(file_list)) {
  results <- data.frame()
  setwd("G:/Database/phenotyes/buildGRCh37")
  gwas2.df <- vroom(file_list[i],col_select = c("SNP","A1","A2","N","BETA","SE"))
  colnames(gwas2.df)[colnames(gwas2.df) == "A2"] <- "effect_allele"
  colnames(gwas2.df)[colnames(gwas2.df) == "A1"] <- "A2"
  colnames(gwas2.df)[colnames(gwas2.df) == "effect_allele"] <- "A1"
  gwas2.df$Z=gwas2.df$BETA/gwas2.df$SE
  gwas2.df<-dplyr::select(gwas2.df,c(SNP,A1,A2,N,Z))
  gwas2.df <- gwas2.df %>% drop_na(Z)
  gwas2.df$N<-as.integer(gwas2.df$N)
  
  set.seed(1000)
  colnames(gwas2.df)<-c("SNP","A1","A2","N","Z")
  gwas2.df <- gwas2.df[!duplicated(gwas2.df$SNP), ]
  gc()
  gwas2.df<-as.data.frame(gwas2.df)
  
  res.HDL <- HDL.rg.parallel(gwas1.df, gwas2.df, LD.path, numCores = 10)
  res.HDL

  select <- cbind(res.HDL$rg, res.HDL$rg.se, res.HDL$P)

  results <- rbind(results, select)
  colnames(results) <- c("rg","rg.se","P")
  results$phenotypes <- file_list[i]
  results_df <- rbind(results_df, results)
  
  rm(gwas2.df,res.HDL)
  gc()
}
write_tsv(results_df,file = "G:/A/phenotypes/HDL_GBM.txt")
