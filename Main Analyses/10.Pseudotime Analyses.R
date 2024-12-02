pacman::p_load(CytoTRACE2, Seurat, ggplot2, gridExtra, tidyverse, patchwork, ggpubr, sceasy, reticulate, monocle3)
### CytoTRACE2----
setwd("D:/single_eqtl")
tempdir()
tempfile()
tempdir <- function() "D:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
A<-readRDS("GBM_core.rds")
A <- cytotrace2(A, 
                is_seurat = TRUE, 
                slot_type = "counts", 
                species = 'human',
                seed = 1234)
annotation <- data.frame(phenotype = A@meta.data$annotation_level_4) %>% 
  set_rownames(., colnames(A))
gc()
A$Phenotype <- A$annotation_level_4

message("Creating plots.")
plot_list <- list()
labels <- c("Differentiated", "Unipotent", "Oligopotent", 
            "Multipotent", "Pluripotent", "Totipotent")
colors <- c("#9E0142", "#F46D43", "#FEE08B", "#E6F598", 
            "#66C2A5", "#5E4FA2")
x_limits <- range(A@reductions$umap@cell.embeddings[, 
                                                    1], na.rm = TRUE)
y_limits <- range(A@reductions$umap@cell.embeddings[, 
                                                    2], na.rm = TRUE)
A@meta.data[["CytoTRACE2_Score_clipped"]] <- 5.5 - 
  6 * A@meta.data[["CytoTRACE2_Score"]]
A@meta.data[["CytoTRACE2_Score_clipped"]] <- -pmax(pmin(A@meta.data[["CytoTRACE2_Score_clipped"]], 
                                                        5), 0)
potency_score_umap <- FeaturePlot(A, "CytoTRACE2_Score_clipped") + 
  scale_colour_gradientn(colours = rev(colors), 
                         na.value = "transparent", labels = c(labels), 
                         limits = c(-5, 0), name = "Potency score \n", 
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black")) + xlab("UMAP1") + 
  ylab("UMAP2") + ggtitle("CytoTRACE 2") + theme(legend.text = element_text(size = 10), 
                                                 legend.title = element_text(size = 12), axis.text = element_text(size = 12), 
                                                 axis.title = element_text(size = 12), plot.title = element_text(size = 12, 
                                                                                                                 face = "bold", hjust = 0.5, margin = margin(b = 20))) + 
  theme(aspect.ratio = 1) + coord_cartesian(xlim = x_limits, 
                                            ylim = y_limits)
plot_list <- c(plot_list, setNames(list(potency_score_umap), 
                                   paste("CytoTRACE2_UMAP")))
potency_category_umap <- DimPlot(A, reduction = "umap", 
                                 group.by = "CytoTRACE2_Potency", label = FALSE) + 
  scale_color_manual(values = colors, name = "Potency category", 
                     breaks = rev(c("Differentiated", "Unipotent", 
                                    "Oligopotent", "Multipotent", "Pluripotent", 
                                    "Totipotent"))) + xlab("UMAP1") + ylab("UMAP2") + 
  ggtitle("CytoTRACE 2") + theme(legend.text = element_text(size = 10), 
                                 legend.title = element_text(size = 12), axis.text = element_text(size = 12), 
                                 axis.title = element_text(size = 12), plot.title = element_text(size = 12, 
                                                                                                 face = "bold", hjust = 0.5, margin = margin(b = 20))) + 
  theme(aspect.ratio = 1) + coord_cartesian(xlim = x_limits, 
                                            ylim = y_limits)
plot_list <- c(plot_list, setNames(list(potency_category_umap), 
                                   paste("CytoTRACE2_Potency_UMAP")))
rel_order_umap <- FeaturePlot(A, "CytoTRACE2_Relative") + 
  scale_colour_gradientn(colours = (c("#000004FF", 
                                      "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", 
                                      "#FCFDBFFF")), na.value = "transparent", limits = c(0, 
                                                                                          1), breaks = seq(0, 1, by = 0.2), labels = c("0.0 (More diff.)", 
                                                                                                                                       "0.2", "0.4", "0.6", "0.8", "1.0 (Less diff.)"), 
                         name = "Relative\norder \n", guide = guide_colorbar(frame.colour = "black", 
                                                                             ticks.colour = "black")) + ggtitle("CytoTRACE 2") + 
  xlab("UMAP1") + ylab("UMAP2") + theme(legend.text = element_text(size = 10), 
                                        legend.title = element_text(size = 12), axis.text = element_text(size = 12), 
                                        axis.title = element_text(size = 12), plot.title = element_text(size = 12, 
                                                                                                        face = "bold", hjust = 0.5, margin = margin(b = 20))) + 
  theme(aspect.ratio = 1) + coord_cartesian(xlim = x_limits, 
                                            ylim = y_limits)
plot_list <- c(plot_list, setNames(list(rel_order_umap), 
                                   paste("CytoTRACE2_Relative_UMAP")))

phenotype_umap <- DimPlot(A, reduction = "umap", 
                          group.by = "Phenotype", label = FALSE) + xlab("UMAP1") + 
  ylab("UMAP2") + ggtitle("Phenotypes") + theme(legend.text = element_text(size = 8), 
                                                legend.title = element_text(size = 12), axis.text = element_text(size = 10), 
                                                axis.title = element_text(size = 10), plot.title = element_text(size = 12, 
                                                                                                                face = "bold", hjust = 0.5, margin = margin(b = 20))) + 
  theme(aspect.ratio = 1) + coord_cartesian(xlim = x_limits, 
                                            ylim = y_limits)
plot_list <- c(plot_list, setNames(list(phenotype_umap), 
                                   paste("Phenotype_UMAP")))
mtd <- A@meta.data[c("Phenotype", "CytoTRACE2_Score")]
medians <- mtd %>% group_by(Phenotype) %>% summarise(CytoTRACE2_median_per_pheno = median(CytoTRACE2_Score, 
                                                                                          na.rm = TRUE)) %>% arrange(desc(CytoTRACE2_median_per_pheno))
phenotypes <- unique(mtd$Phenotype) 
medians <- data.frame(Phenotype = character(),
                      CytoTRACE2_median_per_pheno = numeric(),
                      stringsAsFactors = FALSE)

for (phenotype in phenotypes) {
  subset_mtd <- mtd[mtd$Phenotype == phenotype, ]
  median_score <- median(subset_mtd$CytoTRACE2_Score, na.rm = TRUE)
  
  medians <- rbind(medians, data.frame(Phenotype = phenotype,
                                       CytoTRACE2_median_per_pheno = median_score))
}

medians <- medians[order(-medians$CytoTRACE2_median_per_pheno), ]


mtd <- mtd %>% inner_join(medians, by = "Phenotype")

mtd$Phenotype <- factor(mtd$Phenotype, levels = medians$Phenotype)


potencyBoxplot_byPheno <- ggplot(mtd[!is.na(mtd$Phenotype), ], aes(x = Phenotype, y = CytoTRACE2_Score)) +
  geom_boxplot(aes(fill = CytoTRACE2_Score), width = 0.8, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = CytoTRACE2_Score), width = 0.05, height = 0, alpha = 0.5, shape = 21, stroke = 0.1, size = 1) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_fill_gradientn(colors = rev(colors)) +
  scale_color_gradientn(colors = rev(colors)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  labs(x = "Phenotype", y = "Potency score") +
  ggtitle("Developmental potential by phenotype") +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 20)),
        aspect.ratio = 0.8)      
dev.off()   

potencyBoxplot_byPheno <- ggplot(mtd[!is.na(mtd$Phenotype), 
], aes(x = Phenotype, y = CytoTRACE2_Score)) + 
  geom_boxplot(aes(fill = CytoTRACE2_median_per_pheno), 
               width = 0.8, alpha = 0.5, outlier.shape = NA) + 
  geom_jitter(aes(fill = CytoTRACE2_median_per_pheno), 
              width = 0.05, height = 0, alpha = 0.5, shape = 21, 
              stroke = 0.1, size = 1) + theme_classic() + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), 
                     limits = c(0, 1), sec.axis = sec_axis(trans = ~., 
                                                           breaks = seq(0, 1, by = 1/12), labels = c("", 
                                                                                                     "Differentiated", "", "Unipotent", "", 
                                                                                                     "Oligopotent", "", "Multipotent", "", 
                                                                                                     "Pluripotent", "", "Totipotent", ""))) + 
  scale_fill_gradientn(colors = rev(colors), 
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 
                                                                        1), labels = c(labels)) + scale_color_gradientn(colors = rev(colors), 
                                                                                                                        breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 
                                                                                                                                                                         1), labels = c(labels)) + scale_x_discrete(labels = function(x) str_wrap(x, 
                                                                                                                                                                                                                                                  width = 10)) + labs(x = "Phenotype", y = "Potency score") + 
  ggtitle("Developmental potential by phenotype") + 
  theme(legend.position = "None", axis.text = element_text(size = 8), 
        axis.title = element_text(size = 12), legend.text = element_text(size = 12), 
        plot.title = element_text(size = 12, face = "bold", 
                                  hjust = 0.5, margin = margin(b = 20)), 
        axis.ticks.y.right = element_line(color = c("black", 
                                                    NA, "black", NA, "black", NA, "black", 
                                                    NA, "black", NA, "black", NA, "black")), 
        aspect.ratio = 0.8, axis.ticks.length.y.right = unit(0.3, 
                                                             "cm"))
plot_list <- c(plot_list, setNames(list(potencyBoxplot_byPheno), 
                                   paste("CytoTRACE2_Boxplot_byPheno")))


p1 <- plot_list$CytoTRACE2_UMAP
p2 <- plot_list$CytoTRACE2_Potency_UMAP
p3 <- plot_list$CytoTRACE2_Relative_UMAP 
p4 <- plot_list$CytoTRACE2_Boxplot_byPheno
p5 <- plot_list$Phenotype_UMAP
p5
(p1+p2+p3+p5) + plot_layout(ncol = 2)

gc()
p4       
p1 <- ggboxplot(A@meta.data, x = "cell", y = "CytoTRACE2_Score", 
                width = 0.6, color = "black", fill = "annotation_level_4",
                palette = "npg", xlab = F, bxp.errorbar = T, bxp.errorbar.width = 0.5,
                size = 1, outlier.shape = NA, legend = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p1    

### Vector----
cds<-readRDS("GBM_core.rds")
VEC = cds@int_colData$reducedDims$UMAP
colnames(VEC) = c('UMAP_1','UMAP_2')
pbmc <- CreateSeuratObject(counts = as.matrix (cds@assays@data[["counts"]]), project = "GBM", min.cells = 0, min.features = 0)
rm(list = setdiff(ls(), c("VEC","pbmc")))
gc()
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc),npcs = 150)
PCA = pbmc@reductions$pca@cell.embeddings
source('Vector.R')
##https://github.com/jumphone/Vector
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
pdf("vecto_all.pdf", width = 10, height = 8)
vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
dev.off()
### PAGA----
A <- readRDS("GBM_core.rds")
Idents(object = A) <- "annotation_level_2"
A<-subset(x = A, idents = c("Glial-Neuronal","Differentiated-like","Stem-like","Vascular"))
A <- NormalizeData(A, normalization.method = "LogNormalize", scale.factor = 10000)
A <- FindVariableFeatures(A, selection.method = "vst", nfeatures = 2000)
A <- ScaleData(A, features = rownames(A))
A <- RunPCA(A, features = VariableFeatures(object = A))
A@assays$RNA@scale.data <- as.matrix(0)
gc()
sceasy::convertFormat(A, from="seurat", to="anndata",
                      outFile='GBM_core_select.h5ad')
sc <- import("scanpy")

adata_DS1 <- sc$read("GBM_core_select.h5ad")
obsm_keys <- adata_DS1$obsm_keys()
sc$pp$neighbors(adata_DS1, n_neighbors = 20L, use_rep = 'X_pca')

sc$tl$paga(adata_DS1, groups='Cell_type_level_3')
plt <- import("matplotlib")
plt$use("Agg", force = TRUE)
plt$rcParams[["figure.figsize"]] <- list(8, 8) 

save_dir <- "figures"
sc$pl$paga(adata_DS1,
           color = 'Cell_type_level_3',
           fontsize = 7,
           frameon = FALSE,
           save = "DS1_paga_SELECT_cell.pdf")

### OPC------
A <- readRDS("GBM_core.rds")
Idents(object = A) <- "cell_type"
rm(list = setdiff(ls(), c("A")))
gc()
A<-subset(x = A, idents = c("malignant cell","oligodendrocyte precursor cell"))
data <- GetAssayData(A, assay = 'RNA', slot = 'counts')
cell_metadata <- A@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
rm(list = setdiff(ls(), c("cds","gene_annotation","A")))

gc()
cds <- preprocess_cds(cds, num_dim = 12)     
plot_pc_variance_explained(cds) 
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds <- cluster_cells(cds)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(A, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed 
plot_cells(cds, reduction_method="UMAP", color_cells_by="Cell_type_level_3")      
dev.off()
cds <- learn_graph(cds)     
head(colData(cds))
plot_cells(cds,
           color_cells_by = "Cell_type_level_3",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)
dev.off()
cds = order_cells(cds) 
pdf("GBM_core_monocle3_plot_OPC.pdf", width = 10, height = 8)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
dev.off()

saveRDS(cds,file = "GBM_core_monocle3_pseudotime.rds.gz")
rm(list = setdiff(ls(), c("cds")))
gc()
### AC-----
A <- readRDS("GBM_core.rds")
Idents(object = A) <- "cell_type"
rm(list = setdiff(ls(), c("A")))
gc()
A<-subset(x = A, idents = c("malignant cell","astrocyte"))
data <- GetAssayData(A, assay = 'RNA', slot = 'counts')
cell_metadata <- A@meta.data

gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
rm(list = setdiff(ls(), c("cds","gene_annotation","A")))

gc()
cds <- preprocess_cds(cds, num_dim = 50)     
plot_pc_variance_explained(cds) 
cds <- reduce_dimension(cds,preprocess_method = "PCA") 
cds <- reduce_dimension(cds, reduction_method="tSNE")
cds <- cluster_cells(cds)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(A, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  
plot_cells(cds, reduction_method="UMAP", color_cells_by="Cell_type_level_3")      
dev.off()

cds <- learn_graph(cds)     
head(colData(cds))
plot_cells(cds,
           color_cells_by = "Cell_type_level_3",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1.5)
dev.off()
cds = order_cells(cds)  
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
dev.off()

saveRDS(cds,file = "GBM_core_monocle3_Astrocyte.rds.gz")
rm(list = setdiff(ls(), c("cds")))
gc()

Track_genes_sig<-c("EGFR")
plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
plot_genes_in_pseudotime(cds[Track_genes_sig,],
                         color_cells_by="Cell_type_level_3",
                         min_expr=0.5, ncol= 2,cell_size=1.5) + 
  scale_color_manual(values = c("#5CB85C","#337AB7","#FFFFCC","#D9534F","#F0AD4E")) 
