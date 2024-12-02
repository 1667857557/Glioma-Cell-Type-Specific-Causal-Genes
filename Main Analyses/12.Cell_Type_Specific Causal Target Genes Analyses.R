pacman::p_load(
  data.table, biomaRt, vroom, GWAS.utils, ieugwasr, IRanges, readr, tidyr, dplyr, devtools, 
  ggplot2, MungeSumstats, tidyverse, mapgen, MVMR, LDlinkR, coloc, TwoSampleMR, plyr, cisMRcML
)
#####Training Cell-Type-SpecificTWAS Model------
cd /mnt/g/Linux/eQTL_to_TWAS/
  find ./Exc/ -name '*.txt' -exec basename {} \; > file_list.txt
find ./W/ -name '*.RDat' -exec basename {} \; | sed 's/\.RDat$//' > RDat_list.txt
grep -vFf RDat_list.txt file_list.txt > result.txt

process_line() {
  line="$1"
  Rscript ./compute_weights.R \
  --extract ./ldsc/w_hm3.snplist \
  --sumstats ./Exc/"$line" \
  --gcta ./gcta/gcta-1.94.1 \
  --gctb ./gctb/gctb \
  --gctb_ref ./gctb/ukbEURu_hm3_shrunk_sparse/ukbEURu_hm3_chr \
  --plink_ref_chr ./Trans-Phar/1000G_EUR_Phase3_plink/1000G.EUR.QC. \
  --plink_ref_keep ./EUR_1KG_phase3_samples.tsv \
  --ld_blocks ./ldetect/EUR \
  --rscript Rscript \
  --dbslmm ./DBSLMM/software \
  --plink ./plink/plink \
  --PRScs_path ./PRScs/PRScs.py \
  --PRScs_ref_path ./PRScs/ldblk_ukbb_eur \
  --ldpred2_ref_dir ./ldref \
  --output ./W/"$line"
}

export -f process_line
cat result.txt | parallel -j 1 process_line


###Cell-Type-SpecificTWAS Weight Processing----
rm (list=ls ())
gc()
rdat_list<-list.files(path='D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes', pattern='.RDat')
pos<-data.frame(PANEL='Astrocytes',
                WGT=paste0('Astrocytes/',rdat_list),
                ID=gsub('\\..*','',rdat_list))

Genes<-vroom("D:/GENE.txt")
# ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# Genes<-getBM(attributes=c('ensembl_gene_id_version','chromosome_name','start_position','end_position'), mart = ensembl)
# Genes$ensembl_gene_id_version<-gsub('\\..*','',Genes$ensembl_gene_id_version)
# names(Genes)<-c('ID','CHR','P0','P1')
Genes<-Genes[complete.cases(Genes),]
Genes<-Genes[!duplicated(Genes),]

pos<-merge(pos, Genes, all.x=T, by='ID')
for(i in 1:length(rdat_list)){
  load(paste0('D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes/',rdat_list[i]))
  pos$N[i]<-N.tot
  pos$CHR_pos[i]<-snps$V1[1]
  pos$P0_pos[i]<-min(snps$V4)
  pos$P1_pos[i]<-max(snps$V4)
  
}
rdat_list[i]
pos$CHR<-pos$CHR_pos
pos$P0[is.na(pos$P0)]<-pos$P0_pos[is.na(pos$P0)]
pos$P1[is.na(pos$P1)]<-pos$P0_pos[is.na(pos$P1)]

pos<-pos[,c('PANEL','WGT','ID','CHR','P0','P1','N')]
setwd("D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/")
write.table(pos, 'D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.pos', col.names=T, row.names=F, quote=F)

pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils")

eqtl_pos<-fread('D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.pos')
eqtl_h2<-NULL
for(i in 1:nrow(eqtl_pos)){
  load(paste0('D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/',eqtl_pos$WGT[i]))
  eqtl_h2<-rbind(eqtl_h2, data.frame(row=i,
                                     ID=eqtl_pos$ID[i],
                                     h2=hsq[1],
                                     se=hsq[2],
                                     p=hsq.pv))
}

fwrite(eqtl_h2, 'D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.hsq.txt', sep=' ', na='NA', quote=F)
rm (list=ls ())
gc()
eqtl_h2<-vroom("D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.hsq.txt")
eqtl_pos<-fread('D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.pos')
eqtl_h2<- subset(eqtl_h2, p<0.05)
eqtl_pos <- filter(eqtl_pos, ID %in% eqtl_h2$ID)
setwd("D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/")
write.table(eqtl_pos, 'D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes.pos', col.names=T, row.names=F, quote=F)
rm (list=ls ())
gc()

###Cell-Type-SpecificTWAS Weight Processing----
cd /mnt/g/linux/fusion_twas-master
for CHR in {1..22}
do
Rscript FUSION.assoc_test.edit.R \
--sumstats ./A/finngen_R9_C3_GBM.${CHR}.sumstats \
--weights ./brain/cell_specify/Exc_ROSMA.pos \
--weights_dir ./brain/cell_specify/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr ${CHR} \
--out finngen_R9_C3_GBM_Exc_ROSMA${CHR}.dat
done
touch finngen_R9_C3_GBM_Exc_ROSMA_TWAS.txt
for i in {1..22}
do
cat "finngen_R9_C3_GBM_Exc_ROSMA$i.dat" >> finngen_R9_C3_GBM_Exc_ROSMA_TWAS.txt
done

###ACAT-O-----
ACATO <- function(p){
  if (all(is.na(p))) return(NA)
  p <- p[!is.na(p)]
  p[p == 1] <- 1 - 1e-16
  #### check if there are very small non-zero p values
  is.small <- (p < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(tan((0.5 - p) * pi))/length(p)
  } else {
    cct.stat <- sum((1 / p[is.small]) / pi)
    cct.stat <- cct.stat + sum(tan((0.5 - p[!is.small]) * pi))
    cct.stat <- cct.stat/length(p)
  }
  #### check if the test statistic is very large.
  if (cct.stat > 1e+15){
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  pval
}
setwd("G:/Linux/fusion_twas-master/all_result/ALL")
data_files <- list.files(pattern = ".txt")
for (file in data_files) {
  twas<-vroom(file)
  twas$TWAS.P<-as.numeric(twas$TWAS.P)
  twas <- twas[complete.cases(twas$TWAS.P), ]
  twas$P_onesided<-twas$TWAS.P/2
  twas$P_onesided[twas$TWAS.Z < 0]<- 1-twas$P_onesided[twas$TWAS.Z < 0]
  twas_gene<-data.table(ensembl_gene_id=sort(unique(twas$ID)),
                        ACATO.P=aggregate(twas$TWAS.P, list(twas$ID), FUN=ACATO)$x,
                        ACATO.P_one=aggregate(twas$P_onesided, list(twas$ID), FUN=ACATO)$x,
                        mean_Z=aggregate(twas$TWAS.Z, list(twas$ID), FUN=mean)$x)
  twas_gene$ACATO.P_two<-twas_gene$ACATO.P_one
  twas_gene$ACATO.P_two[twas_gene$ACATO.P_one > 0.5]<- 1-twas_gene$ACATO.P_one[twas_gene$ACATO.P_one > 0.5]
  twas_gene$ACATO.P_two<-2*twas_gene$ACATO.P_two
  twas_gene$ACATO.Z<--qnorm(twas_gene$ACATO.P/2)
  twas_gene$ACATO.Z_one<--qnorm(twas_gene$ACATO.P_one)
  output_file <- paste0(sub(".txt", "", file), "_ACATO.txt")
  write_tsv(twas_gene,file =output_file)
}

setwd("G:/Linux/fusion_twas-master/all_result/ALL")
results_df <- data.frame()
data_files <- list.files(pattern = "_ACATO.txt")
for (file in data_files) {
  twas<-vroom(file)
  twas<- subset(twas, ACATO.P<0.05)
  twas$names<-file
  results_df<-rbind(results_df,twas)
}
write_tsv(results_df,file = "GBM_result.txt")
gc()

###SMR Analyses-----
cd /mnt/g/smr-1.3.1-linux-x86_64
export PATH=$PATH:/mnt/g/smr-1.3.1-linux-x86_64

for i in {1..22}
do
smr-1.3.1 --bfile EUR --gwas-summary nonGBM_SMR.ma --beqtl-summary ./BrainMeta_cis_sqtl_summary/BrainMeta_cis_sQTL_chr$i --diff-freq 0.2 --diff-freq-prop 0.5 --out nonGBM_BrainMeta_chr$i --thread-num 12
done
touch nonGBM_BrainMeta_sqtl.smr
for i in {1..22}
do
cat "GBM_BrainMeta_chr$i.smr" >> nonGBM_BrainMeta_sqtl.smr
done
###Cell-Type-Specific Colocalization -----

type1 <- "quant"         
type2 <- "cc"            
s1 <- 1         
s2 <- 6191/24381 
setwd("G:/A/cell_specify/COLOC")
A<-fread("GBM.txt")
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
setwd("G:/A/cell_specify/COLOC")

write.csv(results_df,file = "GBM_significant.csv",row.names = FALSE)
gc()


###cisMRcML----
setwd("G:/A/target/GWAS")

gwas <- read_tsv("GBM_GLIOMA.txt.gz", col_select = c("SNP", "CHR", "POS", "eaf", "N", "SE", "BETA", "P", "other_allele", "effect_allele")) %>%
  select(SNP = SNP, CHR = CHR, A1 = effect_allele, A2 = other_allele, EAF = eaf, BETA = BETA, SE = SE, P = P, N = N)

setwd("G:/Database/QTL_summary/coloc")
A <- read_tsv("cell_type.txt")
data_files <- A$A

for (file in data_files) {
  output_file <- gsub("G:/Database/QTL_summary/coloc/", ".txt", file)
  output_file <- gsub("/", "_", output_file)
  output_file <- gsub(".txt", "", output_file)
  
  s <- output_file
  
  t <- vroom(file, col_select = c("SNP", "CHR", "BP", "FREQ", "N", "SE", "BETA", "P", "other_allele", "effect_allele"))
  chr <- t$CHR[1]
  pqtl_df <- data.frame(SNP = t$SNP, A1 = t$effect_allele, A2 = t$other_allele, freq = t$FREQ, b = t$BETA, se = t$SE, p = t$P, N = t$N)
  
  gwas_df <- gwas %>% filter(SNP %in% pqtl_df$SNP)
  pqtl_df <- pqtl_df %>% filter(SNP %in% gwas_df$SNP)
  gwas_pqtl_df <- merge(gwas_df, pqtl_df, by = 'SNP') 
  ##https://github.com/ZhaotongL/GraphMRcML/blob/60f7dacebd9156238f3d993616964fa3cb16f268/real_data/allele_qc.R
  source("G:/allele_qc.R")
  remove_flip <- allele.qc(gwas_pqtl_df$A1.x, gwas_pqtl_df$A2.x, gwas_pqtl_df$A1.y, gwas_pqtl_df$A2.y)
  flip_snp <- gwas_pqtl_df$SNP[remove_flip$flip]
  gwas_pqtl_df$BETA[which(remove_flip$flip)] <- -gwas_pqtl_df$BETA[which(remove_flip$flip)]
  gwas_pqtl_df$EAF[which(remove_flip$flip)] <- 1 - gwas_pqtl_df$EAF[which(remove_flip$flip)]
  gwas_pqtl_df <- gwas_pqtl_df[remove_flip$keep,]
  
  gwas_df <- gwas_pqtl_df %>% select(SNP, CHR, A1 = A1.y, A2 = A2.y, freq = EAF, b = BETA, se = SE, p = P, N = N.x)
  pqtl_df <- gwas_pqtl_df %>% select(SNP, CHR, A1 = A1.y, A2 = A2.y, freq, b, se, p, N = N.y)
  
  pqtl_fn <- paste0('G:/Linux/GCTA/pqtl/', s, '.txt')
  fwrite(pqtl_df, pqtl_fn, sep = '\t')
  gwas_fn <- paste0('G:/Linux/GCTA/gwas/', s, '.txt')
  fwrite(gwas_df, gwas_fn, sep = '\t')
  
}

setwd("/mnt/g/A/target/GWAS")

gwas <- read_tsv("GBM_GLIOMA.txt.gz", col_select = c("SNP", "CHR", "POS", "eaf", "N", "SE", "BETA", "P", "other_allele", "effect_allele")) %>%
  select(SNP = SNP, CHR = CHR, A1 = effect_allele, A2 = other_allele, EAF = eaf, BETA = BETA, SE = SE, P = P, N = N)
setwd("/mnt/g/Database/QTL_summary/coloc")
A <- read_tsv("cell_type.txt")
data_files <- A$A
for (file in data_files) {
  output_file <- gsub("/mnt/g/Database/QTL_summary/coloc/", ".txt", file)
  output_file <- gsub("/", "_", output_file)
  output_file <- gsub(".txt", "", output_file)
  t <- vroom(file, col_select = c("SNP", "CHR", "BP", "FREQ", "N", "SE", "BETA", "P", "other_allele", "effect_allele"))
  s <- output_file
  chr <- max(t$CHR) 
  
  pqtl_df <- fread(paste0('/mnt/g/Linux/GCTA/pqtl/', s, '.txt'))
  pqtl_fn <- paste0('/mnt/g/Linux/GCTA/pqtl/cojo/', s, '.txt')
  fwrite(pqtl_df[,-2], pqtl_fn, sep = '\t')
  
  cojo_pqtl_cm <- paste0("/mnt/g/Linux/GCTA/gcta64 --bfile /mnt/g/Linux/GCTA/UK_1KG_impu_hrc_2016_LIFT_chr",chr, " --diff-freq 1 --maf 0.01 --cojo-wind 10 --chr ", chr, " --cojo-file ", pqtl_fn, " --cojo-slct --cojo-p 5e-6 --out /mnt/g/Linux/GCTA/pqtl/cojo/", s)
  system(cojo_pqtl_cm, intern = TRUE) 
  if (file.exists(paste0('/mnt/g/Linux/GCTA/pqtl/cojo/',s,'.jma.cojo'))) {
    
    pqtl_cojo_res <- fread(paste0('/mnt/g/Linux/GCTA/pqtl/cojo/', s, '.jma.cojo'))
    if (nrow(pqtl_cojo_res) > 2) {
      write.table(pqtl_cojo_res$SNP, paste0('/mnt/g/Linux/GCTA/pqtl/cojo/', s, '.cojo.snp'), quote = FALSE, row.names = FALSE, col.names = FALSE)
      gwas_df <- fread(paste0('/mnt/g/Linux/GCTA/gwas/', s, '.txt'))
      gwas_fn <- paste0('/mnt/g/Linux/GCTA/gwas/cojo/', s, '.txt')
      gwas_df1<-gwas_df
      fwrite(gwas_df1 %>% mutate(b = log(b), se = b * se) %>% select(-CHR), gwas_fn, sep = '\t')
      cojo_gwas_cm <- paste0("/mnt/g/Linux/GCTA/gcta64 --bfile /mnt/g/Linux/GCTA/UK_1KG_impu_hrc_2016_LIFT_chr", chr, " --diff-freq 1 --maf 0.01 --cojo-wind 10 --chr ", chr, " --cojo-file ", gwas_fn, " --cojo-slct --cojo-p 5e-6 --out /mnt/g/Linux/GCTA/gwas/cojo/", s)
      system(cojo_gwas_cm, intern = TRUE)  
      
      if (file.exists(paste0('/mnt/g/Linux/GCTA/gwas/cojo/', s, '.jma.cojo'))) {
        gwas_cojo_res <- fread(paste0('/mnt/g/Linux/GCTA/gwas/cojo/', s, '.jma.cojo'))
        IV <- union(pqtl_cojo_res$SNP, gwas_cojo_res$SNP)
      } else {
        IV <- pqtl_cojo_res$SNP
      }
      
      write.table(IV, paste0('/mnt/g/Linux/GCTA/pqtl/', s, '.IV'), quote = FALSE, row.names = FALSE, col.names = FALSE)
      joint_pqtl_cm <- paste0("/mnt/g/Linux/GCTA/gcta64 --bfile /mnt/g/Linux/GCTA/UK_1KG_impu_hrc_2016_LIFT_chr", chr, " --diff-freq 1 --maf 0.01 --cojo-p 5e-6 --cojo-wind 10 --chr ", chr, " --cojo-file ", pqtl_fn, " --extract /mnt/g/Linux/GCTA/pqtl/", s, ".IV --cojo-joint --out /mnt/g/Linux/GCTA/pqtl/cojo/", s)
      p <- system(joint_pqtl_cm, intern = TRUE)
      LD_mat = fread(paste0('/mnt/g/Linux/GCTA/pqtl/cojo/', s,'.ldr.cojo'))
      LD_mat = LD_mat[, 2:(ncol(LD_mat) - 1)]
      LD_mat = as.matrix(LD_mat)
      rownames(LD_mat) = colnames(LD_mat)
      pqtl_df <- pqtl_df %>% filter(SNP %in% IV)
      pqtl_df <- pqtl_df[match(colnames(LD_mat), pqtl_df$SNP), ]
      pqtl_df$p <- as.numeric(pqtl_df$p)
      pqtl_df$z <- pqtl_df$b / pqtl_df$se
      pqtl_cor <- pqtl_df$b / sqrt(pqtl_df$b^2 + (pqtl_df$N - 2) * pqtl_df$se^2)
      pqtl_df$cor <- pqtl_cor
      pqtl_df$se_cor <- pqtl_df$se * pqtl_cor / pqtl_df$b
      pqtl_df$bJ <- solve(LD_mat) %*% pqtl_cor
      if (nrow(pqtl_df) > 2) {
        gwas_df <- gwas_df %>% filter(SNP %in% IV)
        gwas_df <- gwas_df[match(colnames(LD_mat), gwas_df$SNP), ]
        gwas_df$p <- as.numeric(gwas_df$p)
        gwas_df$z <- gwas_df$b / gwas_df$se
        
        gwas_cor <- gwas_df$b / sqrt(gwas_df$b^2 + (gwas_df$N - 2) * gwas_df$se^2)
        gwas_df$cor <- gwas_cor
        gwas_df$se_cor <- gwas_df$se * gwas_cor / gwas_df$b
        
        gwas_df$bJ <- solve(LD_mat) %*% gwas_cor
        mr_dat = list(b_exp=pqtl_df$bJ, b_out=gwas_df$bJ, N1=median(pqtl_df$N), N2=median(gwas_df$N), LD_mat = LD_mat,
                      exp_df = pqtl_df, out_df = gwas_df,
                      exp_IV=pqtl_cojo_res$SNP, out_IV=setdiff(IV,pqtl_cojo_res$SNP))
        save(mr_dat, file=paste0('/mnt/g/Linux/GCTA/MRdat/',s,'.RData'))
      }
    }
  } 
}  


setwd("G:/Linux/GCTA/MRdat")
data_files <- list.files(pattern = "*.RData")
results_df <- data.frame()

for (file in data_files) {
  load(file)
  b_exp_cond<-mr_dat$exp_df$bJ
  b_out_cond<-mr_dat$out_df$bJ 
  Sig_exp1 <- solve(mr_dat$LD_mat) %*% (mr_dat$exp_df$se_cor %o% mr_dat$exp_df$se_cor * mr_dat$LD_mat) %*%
    solve(mr_dat$LD_mat)
  Sig_out1 <- solve(mr_dat$LD_mat) %*% (mr_dat$out_df$se_cor %o% mr_dat$out_df$se_cor * mr_dat$LD_mat) %*%
    solve(mr_dat$LD_mat)
  
  ciscML_res = cismr_cML_DP(b_exp=b_exp_cond,b_out=b_out_cond,
                            Sig_exp_inv=solve(Sig_exp1),Sig_out_inv=solve(Sig_out1),maxit=200,
                            n = mr_dat$N1,random_start = 5,
                            min_theta_range=-0.1,max_theta_range=0.1,
                            num_pert=100,random_start_pert=5,random_seed = 12345)
  
  my.res.df<-c(BIC_theta=ciscML_res$BIC_theta,BIC_se=ciscML_res$BIC_se,BIC_p=ciscML_res$BIC_p,BIC_DP_theta=ciscML_res$BIC_DP_theta,BIC_DP_se=ciscML_res$BIC_DP_se,BIC_DP_p=ciscML_res$BIC_DP_p)
  my.res.df<-t(as.data.frame(my.res.df))
  my.res.df<-as.data.frame(my.res.df)
  my.res.df$exposure<-file
  results_df<-rbind(results_df,my.res.df)
  colnames(results_df)<-colnames(my.res.df)
}
write_tsv(results_df,file = "GBM_cojo_10_cisMRcML.txt")