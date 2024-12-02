###Prepare CellChat Database----
pacman::p_load(devtools, sceasy, reticulate, CytoTRACE2, Seurat, vroom, dplyr, data.table, tidyverse, CellChat, NMF, ggalluvial, patchwork, ggplot2, svglite, ComplexHeatmap)

#R
options(stringsAsFactors = FALSE)
setwd("D:/")
interaction_input <- read.csv(file = './cellphonedb-data-master/data/interaction_input.csv')
complex_input <- read.csv(file = './cellphonedb-data-master/data/complex_input.csv', row.names = 1)
geneInfo <- read.csv(file = './cellphonedb-data-master/data/gene_input.csv')
geneInfo$Symbol <- geneInfo$hgnc_symbol
geneInfo <- select(geneInfo, -c("ensembl"))
geneInfo <- unique(geneInfo)
idx_partnerA <- match(interaction_input$partner_a, geneInfo$uniprot)
idx.use <- !is.na(idx_partnerA)
interaction_input$ligand <- interaction_input$partner_a
interaction_input$ligand[idx.use] <- geneInfo$hgnc_symbol[idx_partnerA[idx.use]]
idx_partnerB <- match(interaction_input$partner_b, geneInfo$uniprot)
idx.use <- !is.na(idx_partnerB)
interaction_input$receptor <- interaction_input$partner_b
interaction_input$receptor[idx.use] <- geneInfo$hgnc_symbol[idx_partnerB[idx.use]]
interaction_input$interaction_name <- interaction_input$interactors
interaction_input$interaction_name_2 <- interaction_input$interaction_name
interaction_input$pathway_name <- interaction_input$classification
interaction_input$pathway_name <- gsub(".*by ", "", interaction_input$pathway_name)
interaction_input$annotation <- interaction_input$directionality
interaction_input <- select(interaction_input, -c("partner_a","partner_b","protein_name_a","protein_name_b","interactors","classification","directionality"))
complexsubunits <- dplyr::select(complex_input, starts_with("uniprot"))
for (i in 1:ncol(complexsubunits)) {
  idx_complex <- match(complex_input[,paste0("uniprot_",i)], geneInfo$uniprot)
  idx.use <- !is.na(idx_complex)
  complex_input[idx.use,paste0("uniprot_",i)] <- geneInfo$hgnc_symbol[idx_complex[idx.use]]
}
colnames(complex_input)[1:ncol(complexsubunits)] <- paste0("subunit_",seq_len(ncol(complexsubunits)))
complex_input <- dplyr::select(complex_input, starts_with("subunit"))
other_info <- list(complex = complex_input)
db.new <- updateCellChatDB(db = interaction_input, gene_info = geneInfo, other_info = other_info, trim.pathway = T,merged = T,species_target = "human")
cellchat@DB <- db.new
save(db.new, file = "CellChatDB_cellphonedb.human_user.rda")
###Data Processing-----
#process healthy data
#python
import scanpy as sc
import pandas as pd
from scipy import sparse
import scanpy as sc
import multiprocessing as mp
import numpy as np
import gc

#neuron
adata1 = sc.read_h5ad("D:/0bb62fec-0cf1-46e1-9d10-de65a6d4f814.h5ad")
#non-neuron
adata2 = sc.read_h5ad("D:/cc9bfb86-96ed-4ecd-bcc9-464120fc8628.h5ad")
adata = anndata.concat([adata1, adata2], join='outer')
del adata1
del adata2
gc.collect() 
frontal_lobe_regions = [
  "Cerebral cortex (Cx) - Inferior frontal gyrus (IFG) - Ventrolateral prefrontal cortex - A44-A45",
  "Cerebral cortex (Cx) - Middle frontal gyrus (MFG) - A46",
  "Cerebral cortex (Cx) - Precentral gyrus (PrCG) - Primary motor cortex - M1C",
  "Cerebral cortex (Cx) - Frontal agranular insular cortex - FI",
  "Cerebral cortex (Cx) - Rostral gyrus (RoG) - Dorsal division of MFC - A32",
  "Cerebral cortex (Cx) - Gyrus rectus (ReG) - Medial orbitofrontal cortex - A14",
  "Cerebral cortex (Cx) - Cingulate gyrus, rostral (CgGr) - Ventral division of MFC - A24",
  "Cerebral cortex (Cx) - Subcallosal Gyrus (SCG) - Subgenual (subcallosal) division of MFC - A25"
]

temporal_lobe_regions = [
  "Cerebral cortex (Cx) - Anterior parahippocampal gyrus, posterior part (APH) - Medial entorhinal cortex - MEC",
  "Cerebral cortex (Cx) - Occipitotemporal (fusiform) gyrus, temporal part (FuGt) - Temporal area TF",
  "Cerebral cortex (Cx) - Temporal pole (TP) - Temporopolar area - A38",
  "Cerebral cortex (Cx) - Anterior parahippocampal gyrus (AG) - Lateral entorhinal cortex - LEC",
  "Cerebral cortex (Cx) - Middle Temporal Gyrus - MTG",
  "Cerebral cortex (Cx) - Superior Temporal Gyrus - STG",
  "Cerebral cortex (Cx) - Transverse temporal gyrus (TTG) - Primary auditory cortex - A1C",
  "Cerebral cortex (Cx) - Perirhinal gyrus (PRG) - A35-A36"
]

parietal_lobe_regions = [
  "Cerebral cortex (Cx) - Parietal operculum (PaO) - Gustatory cortex - A43",
  "Cerebral cortex (Cx) - Postcentral gyrus (PoCG) - Primary somatosensory cortex - S1C",
  "Cerebral cortex (Cx) - Supramarginal gyrus (SMG) - A40",
  "Cerebral cortex (Cx) - Supraparietal lobule (SPL) - Posterosuperior (dorsal) parietal cortex - A5-A7"
]

regions_of_interest = frontal_lobe_regions + temporal_lobe_regions + parietal_lobe_regions
adata= adata[adata.obs['dissection'].isin(regions_of_interest)].copy()

adata.write('merge_data.h5ad')

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata = adata[adata.obs.pct_counts_mt < 5, :]
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=21000, subset=False, inplace=True)
adata1 = adata[:, adata.var['highly_variable']]
gc.collect()  

import anndata
adata1 = anndata.AnnData(X=adata1.X, obs=adata1.obs, var=adata1.var)
gc.collect()  

print(adata1)
adata1.write('merge_data_select.h5ad')
#R
use_condaenv('sceasy')

time.py2R = system.time({
  sceasy::convertFormat(obj = "./merge_data_select.h5ad", from="anndata",
                        to="seurat",outFile = 'regional_healthy.rds')
})

#process GBM data
#python
brain_regions  = ["left temporal lobe", "right frontal lobe", "right parietal lobe",
                   "right temporal lobe", "left frontal lobe", "temporoparietal junction",
                   "left parietal lobe", "temporal lobe", "frontal lobe", "parietal lobe"]
adata = sc.read_h5ad("D:/GBmap_expend.h5ad")

adata= adata[adata.obs['tissue'].isin(brain_regions)].copy()
adata.write('regional_GBmap_expend.h5ad')

#R
time.py2R = system.time({
  sceasy::convertFormat(obj = "./regional_GBmap_expend.h5ad", from="anndata",
                        to="seurat",outFile = 'regional_GBmap_expend.rds')
})

###Non-Malignant GBM CellChat------
#R
A<-readRDS("D:/COTEX_BRAIN/regional_GBmap_expend.rds")
Idents(object = A) <- "cell_type"
A<-subset(x = A, !idents = c("malignant cell"))

gene_ids <- A@assays$RNA@counts@Dimnames[[1]]
unique(gene_ids)

B<-vroom("D:/humanGTF",col_names = c("external_gene_name","TYPE","ensembl_gene_id","CHR","START","END"))
# B<-vroom("selected_gene_data.txt")
B$ensembl_gene_id <- gsub("\\.\\d+", "", B$ensembl_gene_id)
gene_info<-dplyr::select(B,c(external_gene_name,ensembl_gene_id))

gene_ids<-as.data.frame(gene_ids)
colnames(gene_ids)[colnames(gene_ids) == "gene_ids"] <- "ensembl_gene_id"
gene_info<-left_join(gene_ids,gene_info,by = "ensembl_gene_id")
gene_info[gene_info == ""] <- NA
gene_info$external_gene_name <- ifelse(is.na(gene_info$external_gene_name), gene_info$ensembl_gene_id, gene_info$external_gene_name)
duplicated_rows <- duplicated(gene_info$external_gene_name)

for (i in 1:length(duplicated_rows)) {
  if (duplicated_rows[i]) {
    duplicated_value <- gene_info$external_gene_name[i]
    count <- sum(duplicated_rows[1:i])
    gene_info$external_gene_name[i] <- paste(duplicated_value, count, sep = "_")
  }
}
RenameGenesSeurat <- function(obj, newnames, gene.use = NULL, de.assay = "RNA") {
  print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features, @reductions$pca@feature.loadings")
  lassays <- Assays(obj)
  assay.use <- if (!is.null(obj@reductions$pca)) obj@reductions$pca@assay.used else de.assay
  DefaultAssay(obj) <- de.assay
  
  if (is.null(gene.use)) {
    all_genenames <- rownames(obj)
  } else {
    all_genenames <- gene.use
    obj <- subset(obj, features = gene.use)
  }
  
  order_name <- function(v1, v2, ref) {
    v2 <- make.names(v2, unique = TRUE)
    df1 <- data.frame(v1, v2)
    rownames(df1) <- df1$v1
    df1 <- df1[ref, ]
    return(df1)
  }
  
  df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(obj))
  all_genenames <- df1$v1
  newnames <- df1$v2
  B<-vroom("D:/humanGTF",col_names = c("external_gene_name","TYPE","ensembl_gene_id","CHR","START","END"))
  B$ensembl_gene_id <- gsub("\\.\\d+", "", B$ensembl_gene_id)
  gene_info<-dplyr::select(B,c(external_gene_name,ensembl_gene_id))
  gene_ids <- A@assays$RNA@counts@Dimnames[[1]]
  if ("SCT" %in% lassays) {
    if ("SCTModel.list" %in% slotNames(obj@assays$SCT)) {
      obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames, ]
      rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
    }
  }
  
  change_assay <- function(a1 = de.assay, obj, newnames = NULL, all_genenames = NULL) {
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
      if (length(RNA@var.features)) {
        df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = RNA@var.features)
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        RNA@var.features <- newnames1
      }
      if (length(RNA@scale.data) && !is.null(RNA@scale.data)) {
        df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(RNA@scale.data))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(RNA@scale.data) <- newnames1
      }
    } else {
      message("Unequal gene sets: nrow(RNA) != nrow(newnames)")
    }
    obj@assays[a1][[1]] <- RNA
    return(obj)
  }
  
  for (a in lassays) {
    DefaultAssay(obj) <- a
    df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = rownames(obj))
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    obj <- change_assay(obj = obj, a1 = a, newnames = newnames1, all_genenames = all_genenames1)
  }
  
  if (!is.null(obj@reductions$pca) && length(obj@reductions$pca)) {
    hvg <- VariableFeatures(obj, assay = assay.use)
    df1 <- order_name(v1 = all_genenames, v2 = newnames, ref = hvg)
    df1 <- df1[rownames(obj@reductions$pca@feature.loadings), ]
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    rownames(obj@reductions$pca@feature.loadings) <- newnames1
  }
  
  if (!is.null(obj[[de.assay]]@meta.features)) {
    try(obj[[de.assay]]@meta.features <- data.frame(row.names = rownames(obj[[de.assay]])))
  }
  return(obj)
}
A<- RenameGenesSeurat(A, gene_info$external_gene_name)
data.input = A@assays$RNA@data
meta = A@meta.data
A <- CreateSeuratObject(counts = data.input,project = "GBM",metadata = meta)
A<-NormalizeData(A,normalization.method = "LogNormalize",scale.factor = 10000)
Idents(object = A) <- "cell_type"
unique(meta$cell_type)
meta <- meta %>%
  mutate(new_cell_type = case_when(
    cell_type %in% c("neuron") ~ "Neuron",
    cell_type %in% c("oligodendrocyte") ~ "Oligodendrocyte",
    cell_type %in% c("oligodendrocyte precursor cell") ~ "OPC",
    cell_type %in% c("astrocyte") ~ "Astrocyte",
    cell_type %in% c("endothelial cell") ~ "Endothelial Cell",
    cell_type %in% c("macrophage") ~ "CNS Macrophage",
    cell_type %in% c("microglial cell") ~ "CNS Macrophage",
    cell_type %in% c("neutrophil") ~ "Leukocyte",
    cell_type %in% c("radial glial cell") ~ "Other cell",
    cell_type %in% c("mural cell") ~ "Mural cell",
    cell_type %in% c("monocyte") ~ "Leukocyte",
    cell_type %in% c("mast cell") ~ "Leukocyte",
    cell_type %in% c("B cell") ~ "Leukocyte",
    cell_type %in% c("natural killer cell") ~ "Leukocyte",
    cell_type %in% c("plasma cell") ~ "Leukocyte",
    cell_type %in% c("mature T cell") ~ "Leukocyte",
    cell_type %in% c("dendritic cell") ~ "Leukocyte",
    TRUE ~ "Other cell"
  ))
meta$nCount_RNA<-rownames(meta)
meta1<-select(meta,c(nCount_RNA,new_cell_type))
data.input = A@assays$RNA@data
cellchat <- createCellChat(object = data.input,
                           meta = meta1,
                           group.by = "new_cell_type")

cellchat <- addMeta(cellchat, meta = meta1)
rm(list = setdiff(ls(), c("cellchat")))
gc()
cellchat <- setIdent(cellchat, ident.use = "new_cell_type")
levels(cellchat@idents)
gc()
load('D:/COTEX_BRAIN/CellChatDB_cellphonedb.human_user.rda')
CellChatDB <- db.new
CellChatDB.use<-CellChatDB
cellchat@DB <-CellChatDB.use
cellchat <- subsetData(cellchat)

options(future.globals.maxSize = 250000 * 1024^5)
gc()

future::plan("multisession", workers = 3) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

options(future.globals.maxSize = 3000000 * 1024^2)
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
saveRDS(cellchat,file= "NOmaligant_gbMAP_EXPAND_cellchat.rds.gz")
rm(list = setdiff(ls(), c("A")))
gc()
###Healthy Control CellChat----
A<-readRDS("merge_data_21000.rds")
gc()
gene_ids <- A@assays$RNA@counts@Dimnames[[1]]
# gene_ids <- sub("\\..*", "", gene_ids)
# A@assays$RNA@counts@Dimnames[[1]]<-gene_ids
# A@assays$RNA@data@Dimnames[[1]]<-gene_ids
unique(gene_ids)

B<-vroom("D:/humanGTF",col_names = c("external_gene_name","TYPE","ensembl_gene_id","CHR","START","END"))
# B<-vroom("selected_gene_data.txt")
B$ensembl_gene_id <- gsub("\\.\\d+", "", B$ensembl_gene_id)
gene_info<-dplyr::select(B,c(external_gene_name,ensembl_gene_id))

gene_ids<-as.data.frame(gene_ids)
colnames(gene_ids)[colnames(gene_ids) == "gene_ids"] <- "ensembl_gene_id"
gene_info<-left_join(gene_ids,gene_info,by = "ensembl_gene_id")
gene_info[gene_info == ""] <- NA
gene_info$external_gene_name <- ifelse(is.na(gene_info$external_gene_name), gene_info$ensembl_gene_id, gene_info$external_gene_name)
duplicated_rows <- duplicated(gene_info$external_gene_name)

for (i in 1:length(duplicated_rows)) {
  if (duplicated_rows[i]) {
    duplicated_value <- gene_info$external_gene_name[i]
      count <- sum(duplicated_rows[1:i])
      gene_info$external_gene_name[i] <- paste(duplicated_value, count, sep = "_")
  }
}
A<- RenameGenesSeurat(A, gene_info$external_gene_name)
A<-NormalizeData(A,normalization.method = "LogNormalize",scale.factor = 10000)
gc()
unique(A$supercluster_term)
data.input = A@assays$RNA@data
meta = A@meta.data
unique(meta$cell_type)
meta <- meta %>%
  mutate(new_cell_type = case_when(
    cell_type %in% c("neuron") ~ "Neuron",
    cell_type %in% c("oligodendrocyte") ~ "Oligodendrocyte",
    cell_type %in% c("oligodendrocyte precursor cell") ~ "OPC",
    cell_type %in% c("astrocyte") ~ "Astrocyte",
    cell_type %in% c("endothelial cell") ~ "Endothelial Cell",
    cell_type %in% c("central nervous system macrophage") ~ "CNS Macrophage",
    cell_type %in% c("fibroblast") ~ "Other cell",
    cell_type %in% c("pericyte") ~ "Mural cell",
    cell_type %in% c("leukocyte") ~ "Leukocyte",
    cell_type %in% c("vascular associated smooth muscle cell") ~ "Mural cell",
    cell_type %in% c("choroid plexus epithelial cell") ~ "Other cell",
    TRUE ~ "Other cell"
  ))
meta = A@meta.data
meta$nCount_RNA<-rownames(meta)
meta1<-select(meta,c(nCount_RNA,new_cell_type))
data.input = A@assays$RNA@data
cellchat <- createCellChat(object = data.input,
                           meta = meta1,
                           group.by = "new_cell_type")

cellchat <- addMeta(cellchat, meta = meta1)
rm(list = setdiff(ls(), c("cellchat")))
gc()
cellchat <- setIdent(cellchat, ident.use = "new_cell_type")
levels(cellchat@idents)
gc()
load('D:/COTEX_BRAIN/CellChatDB_cellphonedb.human_user.rda')
CellChatDB <- db.new
CellChatDB.use<-CellChatDB
cellchat@DB <-CellChatDB.use
cellchat <- subsetData(cellchat)

options(future.globals.maxSize = 250000 * 1024^5)
gc()

future::plan("multisession", workers = 3) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

options(future.globals.maxSize = 3000000 * 1024^2)
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
saveRDS(cellchat,file= "merge_data_21000_cellchat.rds.gz")
rm(list = setdiff(ls(), c("A")))

gc()
###Compare Multiple CellChats ----
cellchat1<-readRDS("CellChatDB_cellphonedb_cellchat_frontal_lobe_norm_brain_process.rds.gz")
cellchat2<-readRDS("CellChatDB_cellphonedb_cellchat_frontal_lobe_expand_GBmap.rds.gz")
dim_info <- dim(cellchat1@data)
dim_names <- dimnames(cellchat1@data)

zero_matrix <- Matrix(0, nrow = dim_info[1], ncol = dim_info[2], sparse = TRUE)

dimnames(zero_matrix) <- dim_names

cellchat1@data <- zero_matrix
dim_info <- dim(cellchat2@data)
dim_names <- dimnames(cellchat2@data)
zero_matrix <- Matrix(0, nrow = dim_info[1], ncol = dim_info[2], sparse = TRUE)

dimnames(zero_matrix) <- dim_names
cellchat2@data <- zero_matrix

object.list <- list(Control = cellchat1, GBM = cellchat2)
rm(cellchat1,cellchat2)
gc()
object.list<-load("G:/A/cell_specify/cellchat/cellchat_frontal_lobe_object.list_GBM_NL_LS.RData")
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gc()
cellchat
save(object.list, file = "cellchat_frontal_lobe_object.list_GBM_NL_LS.RData")
save(cellchat, file = "cellchat_frontal_lobe_merged_GBM_NL_LS.RData")
setwd("G:/A/cell_specify/cellchat/")
object.list<-readRDS("cellchat_frontal_lobe_object.list_GBM_NL_LS.rds.gz")
cellchat<-readRDS("cellchat_frontal_lobe_merged_GBM_NL_LS.rds.gz")
ptm = Sys.time()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

###All GBM Map Expand Regional CellChat----
A<-readRDS("D:/COTEX_BRAIN/regional_GBmap_expend.rds")
Idents(object = A) <- "tissue"
rm(list = setdiff(ls(), c("A")))
gc()
gene_ids <- A@assays$RNA@counts@Dimnames[[1]]
unique(gene_ids)
B<-vroom("D:/humanGTF",col_names = c("external_gene_name","TYPE","ensembl_gene_id","CHR","START","END"))
B$ensembl_gene_id <- gsub("\\.\\d+", "", B$ensembl_gene_id)
gene_info<-dplyr::select(B,c(external_gene_name,ensembl_gene_id))

gene_ids<-as.data.frame(gene_ids)
colnames(gene_ids)[colnames(gene_ids) == "gene_ids"] <- "ensembl_gene_id"
gene_info<-left_join(gene_ids,gene_info,by = "ensembl_gene_id")
gene_info[gene_info == ""] <- NA
gene_info$external_gene_name <- ifelse(is.na(gene_info$external_gene_name), gene_info$ensembl_gene_id, gene_info$external_gene_name)
duplicated_rows <- duplicated(gene_info$external_gene_name)

for (i in 1:length(duplicated_rows)) {
  if (duplicated_rows[i]) {
    duplicated_value <- gene_info$external_gene_name[i]
    count <- sum(duplicated_rows[1:i])
    gene_info$external_gene_name[i] <- paste(duplicated_value, count, sep = "_")
  }
}
A<- RenameGenesSeurat(A, gene_info$external_gene_name)
data.input = A@assays$RNA@data
meta = A@meta.data

A <- CreateSeuratObject(counts = data.input,project = "GBM",metadata = meta)

A<-NormalizeData(A,normalization.method = "LogNormalize",scale.factor = 10000)


Idents(object = A) <- "cell_type"
unique(meta$cell_type)
meta <- meta %>%
  mutate(new_cell_type = case_when(
    cell_type %in% c("malignant cell") ~ "Malignant cell",
    cell_type %in% c("neuron") ~ "Neuron",
    cell_type %in% c("oligodendrocyte") ~ "Oligodendrocyte",
    cell_type %in% c("oligodendrocyte precursor cell") ~ "OPC",
    cell_type %in% c("astrocyte") ~ "Astrocyte",
    cell_type %in% c("endothelial cell") ~ "Endothelial Cell",
    cell_type %in% c("macrophage") ~ "CNS Macrophage",
    cell_type %in% c("microglial cell") ~ "CNS Macrophage",
    cell_type %in% c("neutrophil") ~ "Leukocyte",
    cell_type %in% c("radial glial cell") ~ "Other cell",
    cell_type %in% c("mural cell") ~ "Mural cell",
    cell_type %in% c("monocyte") ~ "Leukocyte",
    cell_type %in% c("mast cell") ~ "Leukocyte",
    cell_type %in% c("B cell") ~ "Leukocyte",
    cell_type %in% c("natural killer cell") ~ "Leukocyte",
    cell_type %in% c("plasma cell") ~ "Leukocyte",
    cell_type %in% c("mature T cell") ~ "Leukocyte",
    cell_type %in% c("dendritic cell") ~ "Leukocyte",
    TRUE ~ "Other cell"
  ))
meta$nCount_RNA<-rownames(meta)
meta1<-select(meta,c(nCount_RNA,new_cell_type))
data.input = A@assays$RNA@data
cellchat <- createCellChat(object = data.input,
                           meta = meta1,
                           group.by = "new_cell_type")

cellchat <- addMeta(cellchat, meta = meta1)
rm(list = setdiff(ls(), c("cellchat")))
gc()
cellchat <- setIdent(cellchat, ident.use = "new_cell_type")
levels(cellchat@idents)
gc()
load('D:/COTEX_BRAIN/CellChatDB_cellphonedb.human_user.rda')
CellChatDB <- db.new
showDatabaseCategory(CellChatDB)
CellChatDB.use<-CellChatDB
cellchat@DB <-CellChatDB.use
cellchat <- subsetData(cellchat)

options(future.globals.maxSize = 250000 * 1024^5)
gc()

future::plan("multisession", workers = 3) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

options(future.globals.maxSize = 3000000 * 1024^2)
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
saveRDS(cellchat,file = "maligant_gbmap_regional_expand.rds.gz")

A<-readRDS("maligant_gbmap_regional_expand.rds.gz")
cellchat<-A
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchatV2_sig_interCount <- function(cellchat, 
                                         count=T, 
                                         source_celltype,
                                         showInter=F,
                                         color_set=NULL,
                                         celltype_order=NULL,
                                         celltype.size=F,
                                         flipped=T
){
  
  require("reshape2")
  require("ggplot2")
  require("ggraph")
  require("tidygraph")
  require("dplyr")
  require("igraph")
  require("CellChat")
  
  if(count==T){
    mat <- as.data.frame(cellchat@net$count)
  }else{
    mat <- as.data.frame(cellchat@net$weight)
  }
  
  
  sourecell <- which(colnames(mat)==source_celltype)
  mat <- mat[order(mat[,sourecell], decreasing = TRUE),]
  
  group.size <- as.data.frame(table(cellchat@idents))
  
  
  df <- data.frame(from = rep(source_celltype,nrow(mat)),
                   to = rownames(mat),
                   inter = mat[,source_celltype])
  
  
  nodes <- data.frame(name = df$to)
  nodes$inter <- df$inter
  
  if(celltype.size==F){
    
    size = rep(5, nrow(mat))
    nodes$size <- size
    nodes <- nodes[order(nodes[,"inter"], decreasing = TRUE),]
    
  }else{
    
    colnames(group.size) <- c('name',"size")
    nodes = merge(nodes, group.size, by='name', all=F)
    
    nodes <- nodes[order(nodes[,"size"], decreasing = TRUE),]
    
  }
  
  
  edges <- df[c("from","to","inter")]
  
  if(is.null(color_set)){
    
    color.use <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B", "#FEE500","#8A9FD1","#C06CAB", "#D8A767",
                   "#90D5E4", "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8", "#6E4B9E","#0C727C", "#7E1416")
    
    color.use <- color.use[1:nrow(mat)]
    names(color.use) <- df$to
    
  }else{
    
    color.use <- color_set
    
    if(is.null(celltype_order)){
      
      names(color.use) <- df$to
      
    }else{
      
      names(color.use) <- celltype_order
    }
    
  }
  
  
  net <- tbl_graph(nodes = nodes, edges = edges)
  
  
  #plot
  p=ggraph(net,layout='igraph', algorithm = 'circle') +
    geom_edge_bend(mapping = aes(edge_width = inter),
                   strength = 0.2,alpha = 0.8,
                   flipped =flipped, edge_color = "#A9AAAA",
                   n=50, show.legend = F,
                   check_overlap =T)+
    geom_edge_loop(aes(edge_width = inter,
                       direction = (from - 1)*360 / length(net)),
                   colour = "#A9AAAA",
                   alpha = 0.5, show.legend = F)+
    scale_edge_width_continuous(range = c(0,5))
  
  
  if(showInter==F){
    
    p = p+geom_node_point(aes(size=size,colour = name), show.legend = F) +
      geom_node_point(aes(size=size), show.legend = F,
                      shape=21,colour = 'black',stroke = 1.5)+
      geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                     angle=0,hjust=0, size=3) + # 设置点的注释
      scale_size_continuous(range = c(1, 15))+
      scale_color_manual(values = color.use)+
      theme_graph()+
      theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
    
  }else{
    
    if(celltype.size==T){
      
      p = p+geom_node_point(aes(size=size,colour = inter)) +
        geom_node_point(aes(size=size), show.legend = F,
                        shape=21,colour = 'black',stroke = 1.5)+
        geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                       angle=0,hjust=0, fontface="bold",size=4) + # 设置点的注释
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100))+
        scale_size_continuous(range = c(1, 10))+
        theme_graph()+
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      
    }else{
      
      p = p+geom_node_point(aes(size=inter,colour = inter)) +
        geom_node_point(aes(size=inter), show.legend = F,
                        shape=21,colour = 'black',stroke = 1.5)+
        geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                       angle=0,hjust=0, fontface="bold",size=3) + # 设置点的注释
        scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100))+
        scale_size_continuous(range = c(1, 10))+
        theme_graph()+
        theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
      
    }
    
  }
  
  return(p)
  
}
unique(A@meta$new_cell_type)
A@meta$cell_type<-A@meta$new_cell_type
final_plot <-cellchatV2_sig_interCount(cellchat = A,
                                          source_celltype = "Malignant cell",
                                          count = TRUE,
                                          celltype.size = TRUE,
                                          showInter = TRUE,
                                          flipped = TRUE)
ggsave("output_image.pdf", final_plot, width = 10, height = 8, dpi = 600)
