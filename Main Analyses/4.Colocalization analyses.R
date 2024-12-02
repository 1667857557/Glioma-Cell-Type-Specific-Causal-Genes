### COLOC Analyses-----
pacman::p_load("coloc","GWAS.utils","vroom","LDlinkR","data.table","readr","tidyr","dplyr","plyr","devtools","ggplot2")
tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
gc()
type1 <- "quant"
type2 <- "cc"   
s1 <- 1       
s2 <- 12496/30686 
setwd("G:/A/target/coloc")
A<-fread("GLIOMA.txt")
data_files <- A$ID
setwd("G:/A/target/GWAS")
gwas2<- read_tsv("GBM_GLIOMA.txt.gz",col_select = c("SNP","CHR","POS","eaf","N","SE","BETA","P","other_allele","effect_allele"))
sample_size2 <- max(gwas2$N) 
gwas2$eaf<-eaf2maf(eaf = gwas2$eaf)
colnames(gwas2)<-c("SNP","chr","pos","maf","n","se","beta","pval","oa","ea")
results_df <- data.frame()
for (file in data_files) {
  gwas1<-vroom(file,col_select = c("SNP","CHR","BP","FREQ","N","SE","BETA","P","other_allele","effect_allele"))
  gwas1$FREQ<-eaf2maf(eaf = gwas1$FREQ)
  colnames(gwas1)<-c("SNP","chr","pos","maf","N","se","beta","pval","oa","ea")
  sample_size1 <- max(gwas1$N) 
  gwas1 <- gwas1[is.finite(gwas1$se), ]
  top <- gwas1
  tab1 <- gwas1
  tab2 <- subset(gwas2, chr == top$chr[1] &
                   pos > (min(top$pos)) & 
                   pos < (max(top$pos))) %>% 
    subset(., !duplicated(SNP))
  
  commonsnps <- tab1$SNP[tab1$SNP %in% tab2$SNP]
  tab1 <- tab1[tab1$SNP %in% commonsnps, ] %>% dplyr::arrange(SNP)
  tab2 <- tab2[tab2$SNP %in% commonsnps, ] %>% dplyr::arrange(SNP)
  if (!setequal(tab1$SNP, tab2$SNP)) {
    stop("SNPs in tab1 and tab2 do not match")
  }
  order1 <- order(tab1$SNP)
  order2 <- order(tab2$SNP)
  tab1 <- tab1[order1, ]
  tab2 <- tab2[order2, ]
  tab1 <- tab1 %>% format_data(type = "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "maf", effect_allele_col = "ea", other_allele_col = "oa", pval_col = "pval")
  tab2 <- tab2 %>% format_data(type = "outcome", snps = tab1$SNP, snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "maf", effect_allele_col = "ea", other_allele_col = "oa", pval_col = "pval")
  dat <- harmonise_data(tab1, tab2, action = 2)
  colnames(dat)
  if (nrow(dat) > 0) {
    
    tab1<-dat[,c("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","beta.exposure",
                 "se.exposure","pval.exposure","eaf.exposure")]
    
    colnames(tab1)<-gsub(pattern = ".exposure", replacement = "", x = colnames(tab1))
    
    tab2<-dat[,c("SNP","chr.outcome","pos.outcome","effect_allele.outcome","other_allele.outcome","beta.outcome",
                 "se.outcome","pval.outcome","eaf.outcome")]
    colnames(tab2)<-gsub(pattern = ".outcome", replacement = "", x = colnames(tab2))
    
    index <- as.character(tab1$effect_allele) == as.character(tab2$effect_allele) & 
      as.character(tab1$other_allele) == as.character(tab2$other_allele) & as.character(tab1$SNP) == 
      as.character(tab2$SNP) & tab1$pos == tab2$pos
    if (sum(index) > 0) {
      tab1$sample_size1 <- sample_size1
      tab2$sample_size2 <- sample_size2
      tab1 <- tab1[index, ] %>% {
        list(pvalues = .$pval, N = .$sample_size1, MAF = .$eaf, beta = .$beta, 
             varbeta = .$se^2, type = type1, snp = .$SNP, z = .$beta/.$se, 
             chr = .$chr, pos = .$pos)
      }
      tab2 <- tab2[index, ] %>% {
        list(pvalues = .$pval, N = .$sample_size2, MAF = .$eaf, beta = .$beta, 
             varbeta = .$se^2, type = type2, snp = .$SNP, z = .$beta/.$se, 
             chr = .$chr, pos = .$pos)
      }
      if (type1 == "cc") {
        tab1$s <- s1 
      }
      if (type2 == "cc") {
        tab2$s <- s2
      }
      
      out <-list(dataset1 = tab1, dataset2 = tab2)
      res <- coloc::coloc.abf(out[[1]], out[[2]])
      my.res.df1<-t(as.data.frame(res[["priors"]]))
      my.res.df2<-t(as.data.frame(res[["summary"]]))
      my.res.df<- as.data.frame(cbind(my.res.df1,my.res.df2))
      my.res.df$exposure<-file
      results_df<-rbind(results_df,my.res.df)
      gc()
    }
    else 
    {next}
  }
  else 
  {next}
}
setwd("G:/A/target/coloc")
write.csv(results_df,file = "GLIOMA_significant.csv",row.names = FALSE)
gc()