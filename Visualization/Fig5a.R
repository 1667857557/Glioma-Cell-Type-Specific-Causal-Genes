pacman::p_load(tidyverse, Seurat, dplyr, locusplotr)

setwd("D:/Work/Ten/figure/figure4")
A<-fread("regional_plot.txt")
A<-na.omit(A)
desired_order <- c("Whole_brain_eQTL", "Cortex_eQTL", "Cerebellum_eQTL",
                   "Astrocytes_eQTL(Fujita et al.)","Astrocytes_eQTL(Bryois et al.)","Brain_pQTL(Wingo et al.)")
A$NAME <- factor(A$NAME, levels = desired_order)
unique(A$NAME)
  pdf("EGFR_locusplot.pdf", width = 10, height = 17)
  
  gg_locusplot(A,genome_build = "GRCh37",plot_pvalue_threshold = 5e-8,
               lead_snp = "rs74504435",
               rsid = SNP,trait=NAME,
               chrom = CHR,
               # plot_genes = T,
               effect = BETA,
               std_err = SE,
               pos = BP,
               ref = effect_allele,
               alt = other_allele,
               # plot_recombination = TRUE,
               p_value = P,population = "EUR",path = "D:/Work/Ten/figure/figure4"
  )
  dev.off()
gc()
