pacman::p_load(sceasy, reticulate, Seurat, vroom, Matrix, preprocessCore, EWCE)

ewce_plot<-function (total_res, mtc_method = "bonferroni", q_threshold = 0.05, 
                     ctd = NULL, annotLevel = 1, heights = c(0.3, 1), make_dendro = FALSE, 
                     verbose = TRUE) 
{
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  check_mtc_method(mtc_method = mtc_method)
  multiList <- TRUE
  if (is.null(total_res$list)) 
    multiList <- FALSE
  if (isTRUE(make_dendro)) {
    if (is.null(ctd)) {
      messager("Warning: Can only add the dendrogram when ctd is provided.", 
               "Setting make_dendro=FALSE.", v = verbose)
      make_dendro <- FALSE
    }
    else {
      if (length(ctd[[annotLevel]]$plotting) > 0) {
        annotLevel <- which(unlist(lapply(ctd, FUN = cells_in_ctd, 
                                          cells = as.character(total_res$CellType))) == 
                              1)
        err_msg2 <- paste0("All of the cells within total_res should come", 
                           " from a single annotation layer of the CTD")
        if (length(annotLevel) == 0) {
          stop(err_msg2)
        }
        cell_ordr <- ctd[[annotLevel]]$plotting$cell_ordering
      }
      else {
        ctdIN <- prep_dendro(ctdIN = ctd[[annotLevel]], 
                             expand = c(0, 0.66))
        cell_ordr <- ctdIN$plotting$cell_ordering
      }
      total_res$CellType <- factor(x = fix_celltype_names(total_res$CellType), 
                                   levels = fix_celltype_names(cell_ordr), ordered = TRUE)
    }
  }
  if (!"q" %in% colnames(total_res)) {
    total_res$q <- stats::p.adjust(total_res$p, method = mtc_method)
  }
  ast_q <- rep("", dim(total_res)[1])
  ast_q[total_res$q < q_threshold] <- "*"
  total_res$ast_q <- ast_q
  total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
  graph_theme <- ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") + 
    ggplot2::theme(text = ggplot2::element_text(size = 14), 
                   axis.title.y = ggplot2::element_text(vjust = 0.6), 
                   strip.background = ggplot2::element_rect(fill = "white"), 
                   strip.text = ggplot2::element_text(color = "black"))
  upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
  total_res$y_ast <- total_res$sd_from_mean * 1.05
  total_res$abs_sd <- abs(total_res$sd_from_mean)
  if ("Direction" %in% colnames(total_res)) {
    the_plot <- ggplot2::ggplot(total_res) + ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", 
                                                                                   y = "abs_sd", fill = "Direction"), position = "dodge", 
                                                               stat = "identity") + graph_theme
  }
  else {
    the_plot <- ggplot2::ggplot(total_res) + ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", 
                                                                                   y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
      ggplot2::scale_fill_gradient(low = "#008B4599", high = "#DF8F44FF") + 
      graph_theme + ggplot2::theme(legend.position = "none")
  }
  the_plot <- the_plot + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 
                                                                      0, 0, 0), "mm"), axis.text.x = ggplot2::element_text(angle = 55, 
                                                                                                                           hjust = 1)) + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", 
                                                                                                                                                                                             fill = NA, linewidth = 1)) + ggplot2::xlab("Cell type") + 
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) + 
    ggplot2::ylab("Std.Devs. from the mean")
  the_plot <- the_plot + ggplot2::scale_y_continuous(breaks = c(0, 
                                                                ceiling(upperLim * 0.66)), expand = c(0, 1.1)) + ggplot2::geom_text(ggplot2::aes_string(label = "ast_q", 
                                                                                                                                                        x = "CellType", y = "y_ast"), size = 10)
  if (isTRUE(multiList)) {
    the_plot <- the_plot + ggplot2::facet_grid("list ~ .", 
                                               scales = "free", space = "free_x")
  }
  output <- list()
  output$plain <- the_plot
  if (isTRUE(make_dendro)) {
    if (length(ctd[[annotLevel]]$plotting) > 0) {
      ctdIN <- prep_dendro(ctdIN = ctd[[annotLevel]], expand = c(0, 
                                                                 0.66))
    }
    output$withDendro <- patchwork::wrap_plots(ctdIN$plotting$ggdendro_horizontal, 
                                               the_plot, heights = heights, ncol = 1)
  }
  return(output)
}
check_mtc_method <- function(mtc_method){
  err_msg <- paste0(
    "ERROR: Invalid mtc_method argument. Please see",
    " '?p.adjust' for valid methods."
  )
  if (!mtc_method %in% c(
    stats::p.adjust.methods
  )) {
    stop(err_msg)
  }
}
###Built EWCE Reference Panel-----
sc_dataset <- readRDS("GBM_core.rds")
Idents(object = sc_dataset) <- "cell_type"
gc()
data.input = sc_dataset@assays$RNA@data
meta = sc_dataset@meta.data
unique(meta$cell_type)
meta <- meta %>%
  mutate(new_cell_type = case_when(
    cell_type %in% c("neuron") ~ "Neuron",
    cell_type %in% c("oligodendrocyte") ~ "Oligodendrocyte",
    cell_type %in% c("oligodendrocyte precursor cell") ~ "OPC",
    cell_type %in% c("astrocyte") ~ "Astrocyte",
    cell_type %in% c("endothelial cell") ~ "Endothelial Cell",
    cell_type %in% c("macrophage") ~ "Other cell",
    cell_type %in% c("microglial cell") ~ "Microglial cell",
    cell_type %in% c("neutrophil") ~ "Other cell",
    cell_type %in% c("radial glial cell") ~ "Other cell",
    cell_type %in% c("mural cell") ~ "Other cell",
    cell_type %in% c("monocyte") ~ "Other cell",
    cell_type %in% c("mast cell") ~ "Other cell",
    cell_type %in% c("B cell") ~ "Other cell",
    cell_type %in% c("natural killer cell") ~ "Other cell",
    cell_type %in% c("plasma cell") ~ "Other cell",
    cell_type %in% c("mature T cell") ~ "Other cell",
    cell_type %in% c("dendritic cell") ~ "Other cell",
    TRUE ~ "Other cell"
  ))
sc_dataset <- CreateSeuratObject(counts = data.input,project = "GBM")
sc_dataset <- AddMetaData(sc_dataset,metadata = meta)
Idents(object = sc_dataset) <- "new_cell_type"
unique(sc_dataset$new_cell_type)
sc_dataset<-NormalizeData(sc_dataset,normalization.method = "LogNormalize",scale.factor = 10000)
sc_dataset<-Seurat::as.SingleCellExperiment(sc_dataset)
gc()
options(future.globals.maxSize = 50 * 1024^3)
exp_CortexOnly_DROPPED <- drop_uninformative_genes(exp=sc_dataset,mtc_method = "BH",
                                                   drop_nonhuman_genes = T,no_cores = 10,
                                                   input_species = "human",method = "gprofiler",output_species = "human",
                                                   level2annot=sc_dataset$new_cell_type)
annotLevels <- list(cell_type=sc_dataset$new_cell_type,
                    annotation_level_4=sc_dataset$annotation_level_4)
ctd <- EWCE::generate_celltype_data(
  exp = exp_CortexOnly_DROPPED,return_ctd=TRUE,
  annotLevels = annotLevels,no_cores = 1,savePath= "D:/",
  groupName = "gbm_select_bonferroni_ENSG") 

annotLevels <- list(cell_type=sc_dataset$new_cell_type,
                    annotation_level_4=sc_dataset$annotation_level_4)

rm(list = setdiff(ls(), c("annotLevels")))
gc()
load("D:/ctd_gbm_select_bonferroni_ENSG.rda")
###Significant Cell-Type-Specific Causal Genes----
hits <- c("SLC4A8", "GALNT6","HEATR3","SOX10", "TMEM184B", "CSNK1E", "BAIAP2L2", 
          "CTA-228A9.3")
reps <- 10000
annotLevel <- 1
unique(ctd[[1]][["annot"]])
unconditional_results <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                         sctSpecies = "human",
                                                         genelistSpecies = "human",
                                                         hits = hits, method = "gprofiler",
                                                         reps = reps,
                                                         annotLevel = annotLevel)
conditional_results_oli <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = hits,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = reps,method = "gprofiler",
  annotLevel = annotLevel,
  controlledCT = "Oligodendrocyte")
conditional_results_opc <- EWCE::bootstrap_enrichment_test(
  sct_data = ctd,
  hits = hits,
  sctSpecies = "human",
  genelistSpecies = "human",
  reps = reps,method = "gprofiler",
  annotLevel = annotLevel,
  controlledCT = "OPC")
merged_results <- rbind(
  data.frame(unconditional_results$results,
             list="Oligodendrocyte-Specific Gene Set (Unconditional)"),
  data.frame(conditional_results_opc$results,
             list="Oligodendrocyte-Specific Gene Set (OPC controlled)")
  # ,data.frame(conditional_results_oli$results,
  # list="Oligodendrocyte-Specific Gene Set (oligodendrocyte controlled)")
)
pdf("EWCE_casual_plots.pdf", width = 6, height = 4)
try({
  plot_list <- ewce_plot(total_res = merged_results, ctd = ctd,
                         mtc_method = "bonferroni") 
  print(plot_list$plain)
})
dev.off()

save(merged_results,file = "merged_results_ensg_10000_limma_fdr.rds.gz")
save(plot_list,file = "plot_list_merged_results_ensg_10000_limma_fdr.rds.gz")
gc()
###Potential Cell-Type-Specific Causal Genes----

gc()

reps <- 10000
annotLevel <- 1
unique(ctd[[1]][["annot"]])
hits <- genes <- c("EGFR", "NARF", "PWWP3A", "ARFRP1", "ZNF512B", "DCT", "ANXA5", "GTF2H2", "SLC1A4", "WWOX", "ZFYVE16", "TBC1D16")
#Astrocyte
unconditional_results_AST <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                             sctSpecies = "human",
                                                             genelistSpecies = "human",
                                                             hits = hits, method = "gprofiler",
                                                             reps = reps,
                                                             annotLevel = annotLevel)
hits <- genes <- c("CDKN2A", "TGFA", "PNPLA7", "EGFR", "TBC1D2", "PCM1", "MTAP", "ARL17B", "GTF2H2", "GRID2", "CEP192")
#OPC
unconditional_results_OPC <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                             sctSpecies = "human",
                                                             genelistSpecies = "human",
                                                             hits = hits, method = "gprofiler",
                                                             reps = reps,
                                                             annotLevel = annotLevel)
hits <- genes <- c("HEATR3", "JAK1", "RAVER2", "NARF", "PCMTD2", "LMF1", "CHRNA4", "CFAP61", "LINC02360", "AHI1", "GMEB2", "CNNM1", "HTR4", "TMEM163", "ARL17B", "SMARCD2", "GTF2H2", "TMX4", "HAR1A")#NEU
#Neuron
unconditional_results_NEU <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                             sctSpecies = "human",
                                                             genelistSpecies = "human",
                                                             hits = hits, method = "gprofiler",
                                                             reps = reps,
                                                             annotLevel = annotLevel)
#Oligodendrocyte
hits <- genes <- c("SLC4A8", "GALNT6", "HEATR3", "SOX10", "CSNK1E", "TMEM184B", "BAIAP2L2", "CTA-228A9.3", "PICK1", "POLR2F", "ZNF148", "PPFIA1", "NUMB", "LHFPL2", "ADAL", "ZNF577", "MTERF3", "PLOD1", "RSU1", "TAF1C", "GXYLT2", "SYNRG", "GTF2H2", "TMEM181")
unconditional_results_OLI <- EWCE::bootstrap_enrichment_test(sct_data = ctd,
                                                             sctSpecies = "human",
                                                             genelistSpecies = "human",
                                                             hits = hits, method = "gprofiler",
                                                             reps = reps,
                                                             annotLevel = annotLevel)
merged_results <- rbind(
  data.frame(unconditional_results_AST$results,
             list="Potential Astrocyte-Specific Gene Set"),
  data.frame(unconditional_results_OPC$results,
             list="Potential OPC-Specific Gene Set"),
  data.frame(unconditional_results_NEU$results,
             list="Potential Neuron-Specific Gene Set"),
  data.frame(unconditional_results_OLI$results,
             list="Potential Oligodendrocyte-Specific Gene Set")
)
dev.off()
pdf("EWCE_potential_plots.pdf", width = 6, height = 8)
pdf("EWCE_potential_plots.pdf", width = 6, height = 6)

try({
  plot_list <- ewce_plot(total_res = merged_results, ctd = ctd,
                         mtc_method = "bonferroni") 
  print(plot_list$plain)
})
dev.off()
merged_results<-readRDS("merged_results_ctd_gbm_fdr_ENSG_release.rds.gz")
setwd("G:/A/cell_specify/EWCE")
saveRDS(merged_results,file = "merged_results_ctd_gbm_fdr_ENSG_release.rds.gz")
saveRDS(plot_list,file = "plot_list_merged_results_ctd_gbm_fdr_ENSG_release.rds.gz")
gc()