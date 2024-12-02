pacman::p_load("vroom", "DOSE", "org.Hs.eg.db", "topGO", "tidyr", "dplyr", "clusterProfiler", "pathview",
               "ggplot2", "forcats", "aplot", "ggtree", "reshape2", "ggsci", "scales", "data.table")
gc()
### GO ----
setwd("G:/A/target/enrichment")
genes <- fread("noGBM.txt",header=T,sep='\t')
gene.df <- bitr(genes$gene, fromType = "ENSEMBL", toType = c("ENTREZID","SYMBOL"), OrgDb = org.Hs.eg.db)
go_bp<-enrichGO(gene = genes$gene,
                OrgDb = org.Hs.eg.db,
                keyType = 'ENSEMBL',
                ont  = "BP",
                pAdjustMethod = "bonferroni",
                pvalueCutoff = 0.05)
go_mf<-enrichGO(gene = genes$gene,
                OrgDb = org.Hs.eg.db,
                keyType = 'ENSEMBL',
                ont  = "MF",
                pAdjustMethod = "bonferroni",
                pvalueCutoff = 0.05)
go_cc<-enrichGO(gene = genes$gene,
                OrgDb = org.Hs.eg.db,
                keyType = 'ENSEMBL',
                ont  = "ALL",
                pAdjustMethod = "bonferroni",
                pvalueCutoff = 0.05)
A<-go_bp@result
A$ONTOLOGY<-"go_bp"
B<-go_cc@result
B$ONTOLOGY<-"go_cc"
A<-rbind(A,B)
colnames(B)
B<-go_mf@result
B$ONTOLOGY<-"go_mf"
A<-rbind(A,B)
write_tsv(A,file = "noGBM_GO.txt")
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = './msigdb_v2023.2.Hs.db')
DBI::dbListTables(con)
geneset_db <- dplyr::tbl(con, 'gene_set')                                             
details_db <- dplyr::tbl(con, 'gene_set_details')                                      
geneset_genesymbol_db <- dplyr::tbl(con, 'gene_set_gene_symbol')                    
genesymbol_db <- dplyr::tbl(con, 'gene_symbol')                                       
collection_db <- dplyr::tbl(con, 'collection') %>% dplyr::select(collection_name, full_name)  
msigdb <- geneset_db %>%
  dplyr::left_join(details_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(collection_db, by = 'collection_name') %>%
  dplyr::left_join(geneset_genesymbol_db, by = c('id' = 'gene_set_id')) %>%
  dplyr::left_join(genesymbol_db, by = c('gene_symbol_id' = 'id')) %>%
  dplyr::select(collection = collection_name, subcollection = full_name, geneset = standard_name, description = description_brief, symbol) %>%
  dplyr::as_tibble() 
DBI::dbDisconnect(con)
#### Prepare GBM enrichment gene set -----

msigdb.go.bp <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:BP") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.bp$geneset <- factor(msigdb.go.bp$geneset)
msigdb.go.bp <- msigdb.go.bp %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.bp$geneset))
list_names <- names(msigdb.go.bp)
list_names<-as.data.frame(list_names)
names_to_select <- c("GOBP_REGULATION_OF_TRANSCRIPTION_INVOLVED_IN_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE",
                     "GOBP_REPLICATIVE_SENESCENCE",
                     "GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION",
                     "GOBP_ERBB2_EGFR_SIGNALING_PATHWAY",
                     "GOBP_ERBB_SIGNALING_PATHWAY",
                     "GOBP_ERBB2_SIGNALING_PATHWAY",
                     "GOBP_POSITIVE_REGULATION_OF_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE",
                     "GOBP_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA1_PRODUCTION",
                     "GOBP_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_ACTIVATION",
                     "GOBP_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION",
                     "GOBP_NEGATIVE_REGULATION_OF_PEPTIDYL_TYROSINE_PHOSPHORYLATION",
                     "GOBP_POSITIVE_REGULATION_OF_PEPTIDYL_TYROSINE_PHOSPHORYLATION",
                     "GOBP_PEPTIDYL_TYROSINE_MODIFICATION",
                     "GOBP_NCRNA_TRANSCRIPTION",
                     "GOBP_CELLULAR_SENESCENCE",
                     "GOBP_PROTEIN_LOCALIZATION_TO_NUCLEOLUS",
                     "GOBP_REGULATION_OF_MIRNA_TRANSCRIPTION","GOBP_POSITIVE_REGULATION_OF_MIRNA_TRANSCRIPTION"
)
selected_bp <- msigdb.go.bp[names_to_select]
msigdb.go.mf <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:MF") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.mf$geneset <- factor(msigdb.go.mf$geneset)
msigdb.go.mf <- msigdb.go.mf %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.mf$geneset))
list_names <- names(msigdb.go.mf)
A<-as.data.frame(list_names)
names_to_select <- c("GOMF_CYCLIN_DEPENDENT_PROTEIN_SERINE_THREONINE_KINASE_INHIBITOR_ACTIVITY","GOMF_P53_BINDING",
                     "GOMF_TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_ACTIVATOR_ACTIVITY")
selected_mf <- msigdb.go.mf[names_to_select]
msigdb.go.cc <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:CC") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.cc$geneset <- factor(msigdb.go.cc$geneset)
msigdb.go.cc <- msigdb.go.cc %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.cc$geneset))
merged_list <- c(selected_bp, selected_mf)
merged_list
pbmc3k.final<-readRDS("G:/A/cell_specify/sc-RNA/glioma/glioma_scRNA.rds")
pbmc3k.final@assays$RNA@counts <- new("dgCMatrix")
A<-pbmc3k.final@meta.data
A$ID<-rownames(A)
write_tsv(A,file = "meta.data_glioma.txt")
gc()

#### Prepare Pan-Glioma Enrichment Gene Set -----
msigdb.go.bp <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:BP") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.bp$geneset <- factor(msigdb.go.bp$geneset)
msigdb.go.bp <- msigdb.go.bp %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.bp$geneset))
list_names <- names(msigdb.go.bp)
list_names<-as.data.frame(list_names)
names_to_select <- c("GOBP_REGULATION_OF_TRANSCRIPTION_INVOLVED_IN_G1_S_TRANSITION_OF_MITOTIC_CELL_CYCLE",
                     "GOBP_POSITIVE_REGULATION_OF_MIRNA_TRANSCRIPTION",
                     "GOBP_CELL_CYCLE_G1_S_PHASE_TRANSITION",
                     "GOBP_REGULATION_OF_PROTEIN_SERINE_THREONINE_KINASE_ACTIVITY",
                     "GOBP_REPLICATIVE_SENESCENCE",
                     "GOBP_REGULATION_OF_CELL_CYCLE_G1_S_PHASE_TRANSITION",
                     "GOBP_REGULATION_OF_MIRNA_TRANSCRIPTION",
                     "GOBP_POSITIVE_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_NCRNA_TRANSCRIPTION",
                     "GOBP_NEGATIVE_REGULATION_OF_PHOSPHORUS_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_NITRIC_OXIDE_SYNTHASE_ACTIVITY",
                     "GOBP_POSITIVE_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_NCRNA_TRANSCRIPTION",
                     "GOBP_REGULATION_OF_PEPTIDYL_TYROSINE_PHOSPHORYLATION",
                     "GOBP_MIRNA_METABOLIC_PROCESS",
                     "GOBP_PEPTIDYL_TYROSINE_MODIFICATION",
                     "GOBP_NCRNA_TRANSCRIPTION",
                     "GOBP_REGULATION_OF_MONOOXYGENASE_ACTIVITY",
                     "GOBP_CELLULAR_SENESCENCE")
selected_bp <- msigdb.go.bp[names_to_select]
selected_bp
msigdb.go.mf <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:MF") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.mf$geneset <- factor(msigdb.go.mf$geneset)
msigdb.go.mf <- msigdb.go.mf %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.mf$geneset))

list_names <- names(msigdb.go.mf)
A<-as.data.frame(list_names)
names_to_select <- c("GOMF_CADHERIN_BINDING")
selected_mf <- msigdb.go.mf[names_to_select]
selected_mf
msigdb.go.cc <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:CC") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.cc$geneset <- factor(msigdb.go.cc$geneset)
msigdb.go.cc <- msigdb.go.cc %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.cc$geneset))
list_names <- names(msigdb.go.cc)
A<-as.data.frame(list_names)
names_to_select <- c("GOCC_PML_BODY","GOCC_TRANSCRIPTION_REPRESSOR_COMPLEX")
selected_cc <- msigdb.go.cc[names_to_select]
selected_cc
merged_list <- c(selected_bp, selected_mf,selected_cc)
merged_list
#### Prepare Non-GBM Enrichment Gene Set -----

msigdb.go.bp <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:BP") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.bp$geneset <- factor(msigdb.go.bp$geneset)
msigdb.go.bp <- msigdb.go.bp %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.bp$geneset))
list_names <- names(msigdb.go.bp)
list_names<-as.data.frame(list_names)
names_to_select <- c("GOBP_NEGATIVE_REGULATION_OF_EPIDERMAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY",
                     "GOBP_NEGATIVE_REGULATION_OF_ERBB_SIGNALING_PATHWAY",
                     "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION",
                     "GOBP_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
                     "GOBP_EPITHELIAL_CELL_PROLIFERATION",
                     "GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY",
                     "GOBP_NCRNA_TRANSCRIPTION",
                     "GOBP_MIRNA_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_NCRNA_TRANSCRIPTION",
                     "GOBP_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
                     "GOBP_REGULATION_OF_MIRNA_TRANSCRIPTION",
                     "GOBP_POSITIVE_REGULATION_OF_MIRNA_METABOLIC_PROCESS",
                     "GOBP_POSITIVE_REGULATION_OF_MIRNA_TRANSCRIPTION"
)
selected_bp <- msigdb.go.bp[names_to_select]
selected_bp
msigdb.go.mf <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:MF") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.mf$geneset <- factor(msigdb.go.mf$geneset)
msigdb.go.mf <- msigdb.go.mf %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.mf$geneset))
list_names <- names(msigdb.go.mf)
A<-as.data.frame(list_names)
names_to_select <- c("GOMF_CYCLIN_DEPENDENT_PROTEIN_SERINE_THREONINE_KINASE_INHIBITOR_ACTIVITY","GOMF_KINASE_REGULATOR_ACTIVITY",
                     "GOMF_TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_ACTIVATOR_ACTIVITY")
selected_mf <- msigdb.go.mf[names_to_select]
selected_mf
msigdb.go.cc <- msigdb %>% 
  dplyr::filter(collection=="C5:GO:CC") %>% 
  dplyr::select(c("geneset", "symbol"))
msigdb.go.cc$geneset <- factor(msigdb.go.cc$geneset)
msigdb.go.cc <- msigdb.go.cc %>% 
  dplyr::group_split(geneset, .keep = F) %>%
  purrr::map( ~.x %>% dplyr::pull(symbol) %>% unique(.)) %>%
  purrr::set_names(levels(msigdb.go.cc$geneset))
list_names <- names(msigdb.go.cc)
A<-as.data.frame(list_names)
names_to_select <- c("GOCC_BASOLATERAL_PLASMA_MEMBRANE","GOCC_BASAL_PART_OF_CELL")
selected_cc <- msigdb.go.cc[names_to_select]
selected_cc
merged_list <- c(selected_bp, selected_mf,selected_cc)
merged_list
rm(list = setdiff(ls(), c("merged_list")))

### scGSEA Analysis ------
pbmc3k.final<-readRDS("G:/process_MARKER_nonGBM.rds.gz")
pbmc3k.final@assays$RNA@counts <- new("dgCMatrix")
gc()
options(future.globals.maxSize = 100000 * 1024^4)

pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", slot = "data",
                             custom = T, geneset = merged_list, ncores = 10,
                             method = c("JASMINE"),
                             kcdf = 'Gaussian')
merged_list
# saveRDS(pbmc3k.final,file ="nonGBM_irGSEA.rds.gz")
result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "cell2",
                               metadata = NULL, col.name = NULL,
                               method = c("JASMINE"))
geneset.show <- result.dge$JASMINE %>% 
  dplyr::filter(p_val_adj <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)

pdf("GO-PATHWAY_GLIOMA.pdf", width = unit(8.5, "cm"), height = unit(11, "cm"), res = 600, compression = "lzw")
irGSEA.heatmap(object = result.dge, 
               method = "JASMINE",heatmap.width = 22,heatmap.heigh = 8,
               show.geneset = geneset.show)

dev.off() 
### scGSEA Analysis ------
inputGenes<-c("FAM181B", "RORB", "RAVER2", "SLC12A7", "TP53", "PLA2G6", "SOX18", "AXIN1", "DLG4", "CDKN2B", "CYLD", "NLGN2", "TERT", "EGFR", "CDKN2A", "JAK1", "STMN3", "TGFA", "PICK1", "RASA4", "PLD5", "MAP2", "KBTBD8", "VSTM2A", "ASAP1", "TCF7L2", "JRKL", "NCAM1", "ATXN7L3B", "EIF5A", "GNG7", "ZGPAT", "NPAS3", "IDH1", "MYC", "ZBTB16", "CBL", "AKT3", "DDX6", "ISL2", "SOX8", "STN1", "GADD45B", "MDM4", "PTPN18", "HEATR3", "ADRM1", "LEPR", "ATP2B4", "PDGFB", "LAMA5", "SPATA24", "PHLDB1", "TPPP", "POLR2F", "SLC16A8", "RHBDF1", "HEXD")
###ALL
gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes = gs,rnaSeqDataset =1)
seEnrichmentOutput<-output[[1]]
enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),row.names = rowData(seEnrichmentOutput)[,1]), colData(seEnrichmentOutput)[,1])
enrichmentOutput$Tissue<-row.names(enrichmentOutput)
head(enrichmentOutput)
ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x='', y = '-LOG10(P-Adjusted)')+
  theme_bw()+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
