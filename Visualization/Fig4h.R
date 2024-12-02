pacman::p_load(
  vroom, readr, mapgen, susieR, tidyverse, ggplot2, GWAS.utils, bigsnpr, 
   data.table, tidyr, dplyr, devtools
)

setwd("C:/Users/16678/OneDrive - sxmu.edu.cn/文档/LocusZooms")
Unique.genes <- read.delim("Gencode_GRCh37_Genes_UniqueList2021.txt", stringsAsFactors = FALSE, header = TRUE)

source("functions/locus_zoom.R")
setwd("G:/CARMA/GBM_genome_wide")
linear<-vroom("CARMA_782.txt")
colnames(linear)
linear<-select(linear,c(CHR,SNP,POS,Pval))
colnames(linear)<-c("CHR","SNP","BP","P")
linear1<-vroom("CARMA_782.txt")

A<-read_rds("CARMA_782.rds.gz")
B<-vroom("CARMA_782.txt")
B$ID<-rownames(B)
D<-A[[1]][["Credible set"]][[2]]
C<-A[[1]][["Credible model"]][1]
B$Credible_model <- ifelse(B$ID %in% unlist(C[[1]]), "causal SNP", "no causal SNP")
B$Credible_set <- ifelse(B$ID %in% unlist(D[[1]]), "causal SNP", "no causal SNP")

setwd("G:/A/target/GWAS")
gwas.sumstats<-vroom("GBM_GLIOMA.txt.gz")
LD_blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
n <- as.numeric(names(sort(table(gwas.sumstats$N),decreasing=TRUE)[1]))
gc()
gwas.sumstats <- process_gwas_sumstats(gwas.sumstats, 
                                       chr='CHR', 
                                       pos='POS', 
                                       beta='BETA', 
                                       se='SE',
                                       a0='other_allele', 
                                       a1='effect_allele', 
                                       snp='SNP', 
                                       pval='P',
                                       LD_Blocks=LD_blocks)

if(max(gwas.sumstats$pval) <= 1){
  gwas.sumstats <- gwas.sumstats %>% dplyr::mutate(pval = -log10(pval))
}
locus <- 782
gwas.sumstats.locus <- gwas.sumstats[gwas.sumstats$locus == locus, ]
LDREF <- load_UKBB_LDREF(LD_blocks, 
                         locus = locus, 
                         LDREF.dir = "G:/Database/fine_mapping", 
                         prefix = "ukb_b37_0.1")
matched.sumstat.LD <- match_gwas_LDREF(gwas.sumstats.locus, LDREF$R, LDREF$var_info)
sumstats.locus <- matched.sumstat.LD$sumstats
R.locus <- matched.sumstat.LD$R
LD_matrices <- list(R.locus)
names(LD_matrices) <- locus

susie.locus.res <- run_finemapping(sumstats.locus, LD_matrices = LD_matrices, priortype = 'uniform', n = n, L = 10)
susie_plot(susie.locus.res[[1]], y='PIP')
susie.locus.sumstats <- merge_susie_sumstats(susie.locus.res, sumstats.locus)
condz <- LD_diagnosis_susie_rss(sumstats.locus$zscore, R = R.locus, n = n)
condz$plot
detected_index <- which(condz$conditional_dist$logLR > 2 & abs(condz$conditional_dist$z) > 2)
cat(sprintf("Detected %d variants with possible allele flipping", length(detected_index)), "\n")
condz$conditional_dist$detected <- 0
condz$conditional_dist$detected[detected_index] <- 1
sumstats.locus.filtered <- sumstats.locus[-detected_index, ]
matched.sumstat.LD <- match_gwas_LDREF(sumstats.locus.filtered, LDREF$R, LDREF$var_info)
sumstats.locus <- matched.sumstat.LD$sumstats
R.locus <- matched.sumstat.LD$R
LD_matrices <- list(R.locus)
names(LD_matrices) <- locus
R<-as.data.frame(matched.sumstat.LD$R)
colnames(R)<-matched.sumstat.LD$sumstats$snp
rownames(R)<-matched.sumstat.LD$sumstats$snp
R<-select(R,c(rs6964933))
R$snp<-rownames(R)
R<-left_join(R,matched.sumstat.LD$sumstats,by="snp")
LD<-vroom("D:/Work/Ten/figure/locus_782_rs6964933.txt")
linear$P <- 10^(-linear$P)
setwd("G:/")
locus.zoom(data = linear,                                   
           region = c(7, 55050293,55300000),                         
           offset_bp = 0,gene = "EGFR",                              
           ld.file = LD,ignore.lead = TRUE,snp="rs6964933",         
           genes.data = Unique.genes,			                 
           plot.title = "Association of rs6964933 in GBM",    
           file.name = "rs6964933.svg",                                
           secondary.snp = c("rs723526","rs1344307","rs1861007","rs6960438","rs759168","rs1024750","rs1024749","rs759166"),         
           secondary.label = F,plot.type = "svg" , nominal = 5) 


