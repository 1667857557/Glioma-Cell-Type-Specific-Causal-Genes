##### Training PWAS/TWAS Model----
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


### PWAS/TWAS Weight Process------
library(data.table)
rm (list=ls ())
gc()
rdat_list<-list.files(path='D:/Work/Ten/brain/brain_cell_TWAS_Fujita et al/Astrocytes', pattern='.RDat')

pos<-data.frame(PANEL='Astrocytes',
                WGT=paste0('Astrocytes/',rdat_list),
                ID=gsub('\\..*','',rdat_list))

library(biomaRt)
library(vroom)
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

### TWAS or PWAS Analyses----
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

### ACAT-O-----
pacman::p_load("RadialMR","ieugwasr","vroom","TwoSampleMR","data.table","readr","tidyr","dplyr","devtools","ggplot2","MendelianRandomization","MungeSumstats","tidyverse","here")
pacman::p_load("bigsnpr","GWAS.utils","ieugwasr","vroom","IRanges","data.table","readr","tidyr","dplyr","devtools","ggplot2","MungeSumstats","tidyverse","mapgen","MVMR")

library(data.table)
ACATO <- function(p){
  if (all(is.na(p))) return(NA)
  p <- p[!is.na(p)]
  p[p == 1] <- 1 - 1e-16
  is.small <- (p < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(tan((0.5 - p) * pi))/length(p)
  } else {
    cct.stat <- sum((1 / p[is.small]) / pi)
    cct.stat <- cct.stat + sum(tan((0.5 - p[!is.small]) * pi))
    cct.stat <- cct.stat/length(p)
  }
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

### SMR Analyses-----
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