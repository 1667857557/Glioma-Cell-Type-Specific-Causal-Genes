###GBM_seperate_scPagwas------
pacman::p_load("vroom","cowplot","data.table",
               "readr","tidyr","dplyr","devtools","ggplot2","biomaRt","tidyverse"
               ,"Seurat","scater","scPagwas","parallel")
Single_data <-readRDS("GBM_core.rds")
Single_data<-NormalizeData(Single_data,normalization.method = "LogNormalize",scale.factor = 10000)
Idents(Single_data) <- "cell_type"

n_split=4
Split_index <- rep(1:n_split, time = ceiling(ncol(Single_data)/n_split), length = ncol(Single_data))

for (i in 1:n_split) {
  Example_splice <- Single_data[,Split_index==i]
  saveRDS(Example_splice,file = paste0("splice",i,".rds"))
}
Pagwas<-list()
gwas_data <- bigreadr::fread2("scPagwas_GBM")
Pagwas <- GWAS_summary_input(
  Pagwas = Pagwas,
  gwas_data = gwas_data
)
Pagwas$snp_gene_df <- SnpToGene(gwas_data = Pagwas$gwas_data, 
                                block_annotation = block_annotation, marg = 10000)

for (i in 1:n_split) {
  scPagwas_main(Pagwas =Pagwas,
                gwas_data =NULL,
                Single_data = paste0("splice",i,".rds"),
                output.prefix=i,
                output.dirs="GBM_cell_type_GBmap",
                Pathway_list=Genes_by_pathway_kegg,
                run_split=TRUE,
                min_clustercells=10,
                assay="RNA",
                block_annotation = block_annotation,
                chrom_ld = chrom_ld)
  gc()
}
output.dirs="GBM_cell_type_GBmap"
oriDir <- paste0("./",output.dirs)
files <- list.files(oriDir, pattern="*_singlecell_scPagwas.gPAS.score.Result.csv")
scPagwas.gPAS.score<-unlist(lapply( files,function(file){
  gs<-read.csv(file=paste0(oriDir,"/",file))
  ga <- gs$scPagwas.gPAS.score
  names(ga) <- gs$cellnames
  return(ga)
}))
Single_data <- Single_data[,names(scPagwas.gPAS.score)]
Single_data$scPagwas.gPAS.score<-scPagwas.gPAS.score
data_mat <- GetAssayData(Single_data, slot = "data", assay = "RNA")
PCC<-scPagwas::Corr_Random(data_mat,
                           scPagwas.gPAS.score,
                           seed=1234,
                           random=T,
                           Nrandom=10,
                           Nselect=200000
)

mean_gpas<-mean(scPagwas.gPAS.score)
a1<-which(scPagwas.gPAS.score >= mean_gpas)
a2<-which(scPagwas.gPAS.score < mean_gpas)

PCC_up <- scPagwas::Corr_Random(scPagwas.gPAS.score=scPagwas.gPAS.score[a1],data_mat=data_mat[,a1])
PCC_down <- scPagwas::Corr_Random(scPagwas.gPAS.score=scPagwas.gPAS.score[a2],data_mat=data_mat[,a2])

scPagwas_topgenes <- names(PCC[order(PCC, decreasing = T)])[1:500]
scPagwas_upgenes <- names(PCC_up)[order(PCC_up, decreasing = T)[1:500]]
scPagwas_downgenes <- names(PCC_down)[order(PCC_down, decreasing = F)[1:500]]
Single_data <- Seurat::AddModuleScore(Single_data, assay = 'RNA', list(scPagwas_topgenes,scPagwas_upgenes,scPagwas_downgenes), name = c("scPagwas.TRS.Score","scPagwas.upTRS.Score","scPagwas.downTRS.Score"))
topgenes<-as.data.frame(scPagwas_topgenes)
upgenes<-as.data.frame(scPagwas_upgenes)
downgenes<-as.data.frame(scPagwas_downgenes)
rm(data_mat)
gc()
correct_pdf <- scPagwas::Get_CorrectBg_p(Single_data=Single_data,
                                         scPagwas.TRS.Score=Single_data$scPagwas.TRS.Score1,
                                         iters_singlecell=100,
                                         n_topgenes=500,
                                         scPagwas_topgenes=scPagwas_topgenes)
gc()
Pagwas$scPagwas.TRS.Score = Single_data$scPagwas.TRS.Score1
Pagwas$Random_Correct_BG_pdf <- correct_pdf
Idents(object = Single_data) <- "cell_type"
Pagwas<-list()
Pagwas$Merged_celltype_pvalue<-scPagwas::Merge_celltype_p(single_p=correct_pdf$pooled_p,
                                                          celltype=Idents(Single_data))
Idents(object = Single_data) <- "annotation_level_4"
Pagwas$Merged_celltype_pvalue_annotation_level_4<-scPagwas::Merge_celltype_p(single_p=correct_pdf$pooled_p,
                                                                             celltype=Idents(Single_data))
Pagwas$topgenes<-scPagwas_topgenes
Pagwas$Merged_celltype_pvalue_annotation_level_4
saveRDS(object = Pagwas, file = "GBmap_core1_Pagwas.rds")
saveRDS(object = Single_data, file = "GBmap_core1_Single_data_Pagwas.rds")
rm (list=ls ())
gc()