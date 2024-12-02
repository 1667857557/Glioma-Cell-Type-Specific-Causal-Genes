pacman::p_load("mapgen", "susieR", "tidyverse", "ggplot2", "GWAS.utils", "RadialMR", 
               "vroom", "TwoSampleMR", "data.table", "readr", "tidyr", "dplyr", 
               "devtools", "MungeSumstats", "bigsnpr", "IRanges", "ChIPseeker", 
               "clusterProfiler", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
               "rtracklayer", "BRGenomics", "GenomicFeatures")
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
setwd("G:/A/target/GWAS")
gwas.sumstats<-vroom("GBM_GLIOMA.txt.gz")
LD_blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
n <- as.numeric(names(sort(table(gwas.sumstats$N),decreasing=TRUE)[1]))
gc()
### Prepare GWAS summary data----
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

# write_tsv(gwas.sumstats,file = "GBM_fine_mapping.txt.gz")
if(max(gwas.sumstats$pval) <= 1){
  gwas.sumstats <- gwas.sumstats %>% dplyr::mutate(pval = -log10(pval))
}
sig.loci <- gwas.sumstats %>% dplyr::filter(pval > -log10(5e-8)) %>% dplyr::pull(locus) %>% unique()
sig.loci
### Diagnosed and excluded unmatched SNPs; running SUSIE fine-mapping----
results_df <- data.frame()
for (value in sig.loci) {
  locus <- value
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
  # ggplot(condz$conditional_dist, aes(x = condmean, y = z, col = factor(detected))) +
  #   geom_point() +
  #   scale_colour_manual(values = c("0" = "black", "1" = "red")) + 
  #   labs(x = "Expected value", y = "Observed z scores", color = "Detected allele flipping") + 
  #   theme_bw()
  sumstats.locus.filtered <- sumstats.locus[-detected_index, ]
  if (length(detected_index) == 0) {
    susie.locus.sumstats <- merge_susie_sumstats(susie.locus.res, sumstats.locus)
  } else {
    sumstats.locus.filtered <- sumstats.locus[-detected_index, ]
    matched.sumstat.LD <- match_gwas_LDREF(sumstats.locus.filtered, LDREF$R, LDREF$var_info)
    sumstats.locus <- matched.sumstat.LD$sumstats
    R.locus <- matched.sumstat.LD$R
    LD_matrices <- list(R.locus)
    names(LD_matrices) <- locus
    
    susie.locus.res <- run_finemapping(sumstats.locus.filtered, LD_matrices = LD_matrices, priortype = 'uniform', n = n, L = 10)
    susie_plot(susie.locus.res[[1]], y='PIP')
    susie.locus.sumstats <- merge_susie_sumstats(susie.locus.res, sumstats.locus.filtered)
    # condz <- LD_diagnosis_susie_rss(sumstats.locus.filtered$zscore, R = R.locus, n = n)
    # condz$plot
  }
  results_df<-rbind(results_df,susie.locus.sumstats)
  rm(LDREF,susie.locus.sumstats,susie.locus.res,R.locus,LD_matrices,sumstats.locus.filtered,condz,matched.sumstat.LD,detected_index)
  gc()
}
write_tsv(results_df,file = "seperated_GBM_UK10K_SNP_fine_mapping_result.txt.gz")
results_df <- readRDS("G:/A/target/fine-mapping/UKBB/seperated_GBM_UK10K_SNP_fine_mapping_result.rds.gz")
finemapstats <- process_finemapping_sumstats(results_df, 
                                             snp = 'snp', 
                                             chr = 'chr', 
                                             pos = 'pos', 
                                             pip = 'susie_pip', 
                                             pval = 'pval', 
                                             zscore = 'zscore', 
                                             cs = 'cs', 
                                             locus = 'locus',  
                                             pip.thresh = 1e-5)

cols.to.keep <- c('snp','chr','pos', 'pip', 'pval', 'zscore', 'cs', 'locus')
finemapstats <- finemapstats[, cols.to.keep]
head(finemapstats, 3)
# gtf_file <- 'G:/Article/glioma/torus_input/gencode.v45lift37.annotation.gtf.gz'
# genomic.annots <- make_genomic_annots(gtf_file)
saveRDS(genomic.annots,file = "gencode.v45lift37.rds.gz")
genomic.annots<-readRDS("gencode.v45lift37.rds.gz")
gene.annots <- genomic.annots$genes

### Generate ATAC annotated----
###PMID:31727856
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
granges_1 <- import.bed("PU1_optimal_peak_IDR_ENCODE.ATAC.bed")
granges_2 <- import.bed("Olig2_optimal_peak_IDR_ENCODE.ATAC.bed")
granges_3 <- import.bed("NeuN_optimal_peak_IDR_ENCODE.ATAC.bed")
granges_4 <- import.bed("LHX2_optimal_peak_IDR_ENCODE.ATAC.bed")
# rm(granges_1,granges_2,granges_3,granges_4)
peaks <- list(Microglia =granges_1, Oligodendrocyte=granges_2,Neuron = granges_3, Astrocyte = granges_4)
head(OCRs, 3)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
Microglia<-peakAnnoList$Microglia@anno
Microglia@ranges@NAMES <- rep("Microglia", length(Microglia@ranges@start))
Microglia@elementMetadata$annotation[Microglia@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
Microglia@elementMetadata$peakType <- Microglia@elementMetadata$annotation
Microglia@elementMetadata$nearestGene <- Microglia@elementMetadata$SYMBOL
Microglia@elementMetadata$SYMBOL <- NULL
Microglia@elementMetadata$annotation <- NULL
head(Microglia, 3)
Oligodendrocyte<-peakAnnoList$Oligodendrocyte@anno
Oligodendrocyte@ranges@NAMES <- rep("Oligodendrocyte", length(Oligodendrocyte@ranges@start))
Oligodendrocyte@elementMetadata$annotation[Oligodendrocyte@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
Oligodendrocyte@elementMetadata$peakType <- Oligodendrocyte@elementMetadata$annotation
Oligodendrocyte@elementMetadata$nearestGene <- Oligodendrocyte@elementMetadata$SYMBOL
Oligodendrocyte@elementMetadata$SYMBOL <- NULL
Oligodendrocyte@elementMetadata$annotation <- NULL
head(Oligodendrocyte, 3)
Neuron<-peakAnnoList$Neuron@anno
Neuron@ranges@NAMES <- rep("Neuron", length(Neuron@ranges@start))
Neuron@elementMetadata$annotation[Neuron@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
Neuron@elementMetadata$peakType <- Neuron@elementMetadata$annotation
Neuron@elementMetadata$nearestGene <- Neuron@elementMetadata$SYMBOL
Neuron@elementMetadata$SYMBOL <- NULL
Neuron@elementMetadata$annotation <- NULL
head(Neuron, 3)
Astrocyte<-peakAnnoList$Astrocyte@anno
Astrocyte@ranges@NAMES <- rep("Astrocyte", length(Astrocyte@ranges@start))
Astrocyte@elementMetadata$annotation[Astrocyte@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
Astrocyte@elementMetadata$peakType <- Astrocyte@elementMetadata$annotation
Astrocyte@elementMetadata$nearestGene <- Astrocyte@elementMetadata$SYMBOL
Astrocyte@elementMetadata$SYMBOL <- NULL
Astrocyte@elementMetadata$annotation <- NULL
head(Astrocyte, 3)
OCRs <- list(Astrocyte,Oligodendrocyte,Microglia,Neuron)
OCRs<-mergeGRangesData(Astrocyte,Oligodendrocyte,Microglia,Neuron, ncores = 1,field = 4)
OCRs<-unlist(as(OCRs, "GRangesList"))
head(OCRs, 3)
setwd("G:/glioma/torus_input")
saveRDS(object = OCRs, file = "brian_ATAC_OCRs.rds")
OCRs <- readRDS("G:/A/target/fine-mapping/brian_ATAC_OCRs.rds")
head(OCRs, 3)
class(genomic.annots$promoters)
active_promoters <- IRanges::subsetByOverlaps(genomic.annots$promoters, 
                                              OCRs, 
                                              minoverlap = 100)

### Generate pcHi-C annotated----
###PMID:31367015
txdb <- loadDb("G:/A/target/fine-mapping/gencode.v45lift37.annotation.gtf.sqlite")
lincRNA_txdb=TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
setwd("G:/A/cell_specify/pcHiC/")
granges_1 <- import.bed("GSM3106832_cortical.cutoff.5.washU.txt")
granges_2 <- import.bed("GSM3598051_astrocyte.cutoff.5.washU.txt")
granges_3 <- import.bed("GSM3598046_hippocampal.cutoff.5.washU.txt")
granges_4 <- import.bed("GSM3598048_motor.cutoff.5.washU.txt")
# rm(granges_1,granges_2,granges_3,granges_4)
peaks <- list(cortical =granges_1, astrocyte=granges_2,hippocampal = granges_3, motor = granges_4)
head(peaks, 3)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,tssRegion=c(-3000, 3000), verbose=FALSE,addFlankGeneInfo=TRUE, flankDistance=5000,annoDb="org.Hs.eg.db")
cortical<-peakAnnoList$cortical@anno
cortical@ranges@NAMES <- rep("cortical", length(cortical@ranges@start))
cortical@elementMetadata$annotation[cortical@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
cortical@elementMetadata$peakType <- cortical@elementMetadata$annotation
cortical@elementMetadata$nearestGene <- cortical@elementMetadata$SYMBOL
cortical@elementMetadata$SYMBOL <- NULL
cortical@elementMetadata$annotation <- NULL
head(cortical, 3)
motor<-peakAnnoList$motor@anno
motor@ranges@NAMES <- rep("motor", length(motor@ranges@start))
motor@elementMetadata$annotation[motor@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
motor@elementMetadata$peakType <- motor@elementMetadata$annotation
motor@elementMetadata$nearestGene <- motor@elementMetadata$SYMBOL
motor@elementMetadata$SYMBOL <- NULL
motor@elementMetadata$annotation <- NULL
head(motor, 3)
hippocampal<-peakAnnoList$hippocampal@anno
hippocampal@ranges@NAMES <- rep("hippocampal", length(hippocampal@ranges@start))
hippocampal@elementMetadata$annotation[hippocampal@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
hippocampal@elementMetadata$peakType <- hippocampal@elementMetadata$annotation
hippocampal@elementMetadata$nearestGene <- hippocampal@elementMetadata$SYMBOL
hippocampal@elementMetadata$SYMBOL <- NULL
hippocampal@elementMetadata$annotation <- NULL
head(hippocampal, 3)
astrocyte<-peakAnnoList$astrocyte@anno
astrocyte@ranges@NAMES <- rep("astrocyte", length(astrocyte@ranges@start))
astrocyte@elementMetadata$annotation[astrocyte@elementMetadata$annotation == "Promoter (<=1kb)"] <- "Promoter"
astrocyte@elementMetadata$peakType <- astrocyte@elementMetadata$annotation
astrocyte@elementMetadata$nearestGene <- astrocyte@elementMetadata$SYMBOL
astrocyte@elementMetadata$SYMBOL <- NULL
astrocyte@elementMetadata$annotation <- NULL
head(astrocyte, 3)
OCRs <- list(astrocyte,motor,cortical,hippocampal)
OCRs<-unlist(as(OCRs, "GRangesList"))
head(OCRs, 3)
setwd("G:/glioma/torus_input")
saveRDS(object = OCRs, file = "GSE113481_brian_pcHiC_enhancer_loops.rds")
gc()
enhancer<-readRDS("GSE113481_brian_pcHiC_enhancer_loops.rds")
head(enhancer, 3)
chr <- seqnames(enhancer)
start <- start(enhancer)
end <- end(enhancer)
promoter_start <- mcols(enhancer)$geneStart
promoter_end <- mcols(enhancer)$geneEnd
gene_name <- mcols(enhancer)$nearestGene
new_object <- data.frame(
  chr = chr,
  start = start,
  end = end,
  promoter_start = promoter_start,
  promoter_end = promoter_end,
  gene_name = gene_name,
  score = NA  
)
new_granges <- GRanges(
  seqnames = new_object$chr,
  ranges = IRanges(start = new_object$start, end = new_object$end),
  promoter_start = new_object$promoter_start,
  promoter_end = new_object$promoter_end,
  gene_name = new_object$gene_name,
  score = new_object$score
)
saveRDS(object = new_granges, file = "GSE113481_brian_pcHiC_OCRs.rds")
gc()
pcHiC<-readRDS("G:/A/cell_specify/ATAC/GSE113481_brian_pcHiC_OCRs.rds")
pcHiC <- pcHiC[pcHiC$gene_name %in% gene.annots$gene_name, ]
head(pcHiC, 3)
genomic.annots$exons_promoters <- list(exons = genomic.annots$exons, 
                                       active_promoters = active_promoters)

### Generate ABC annotated----
###PMID:31784727
setwd("G:/Article/glioma/torus_input")
ABC<-fread("AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
ABC <- ABC[ABC$CellType %in% c("H1_Derived_Neuronal_Progenitor_Cultured_Cells-Roadmap", "bipolar_neuron_from_iPSC-ENCODE", "astrocyte-ENCODE"), ]
ABC <- process_ABC(ABC, full.element = TRUE)
ABC <- ABC[ABC$gene_name %in% gene.annots$gene_name, ]
gc()
setwd("G:/A/target/fine-mapping/")
ABC<-vroom("ABC_for_glioma.txt.gz")
ABC <- process_ABC(ABC, full.element = TRUE)
ABC <- ABC[ABC$gene_name %in% gene.annots$gene_name, ] # restrict to protein coding genes
head(ABC, 3)
gc()

### Combined all annotated----
setwd("G:/A/target/fine-mapping/")
genomic.annots$enhancer_regions <- OCRs[OCRs$peakType!="Promoter",]
nearby20kb <- nearby_interactions(genomic.annots$enhancer_regions,
                                  active_promoters, 
                                  max.dist = 20000)
head(nearby20kb, 3)
genomic.annots$enhancer_loops <- list(ABC = ABC, nearby20kb = nearby20kb, pcHiC = pcHiC)
head(genomic.annots, 3)
summary(genomic.annots)
# saveRDS(genomic.annots,file = "genomic.annots.rds")
# genomic.annots<-readRDS("genomic.annots.rds")
gene.annots <- genomic.annots$genes
rm(ABC,gene.annots,OCRs,promoters,nearby20kb)
rm(active_promoters)
gc()

###Run Gene mapping----
gene.mapping.res <- compute_gene_pip(finemapstats, 
                                     genomic.annots, 
                                     intron.mode = FALSE, 
                                     d0 = 50000,
                                     exon.weight = 1, 
                                     loop.weight = 1)

# write_tsv(gene.mapping.res,file = "fine_mapping_result_GBM.txt")
saveRDS(object = results_df, file = "seperated_GBM_UK10K_SNP_fine_mapping_result_gene.cs.rds")
# gene.mapping.res<-vroom("fine_mapping_result_GBM.txt")
# genomic.annots<-readRDS("genomic.annots.rds")
gene.annots <- genomic.annots$genes
gc()
gene.pip.res <- extract_gene_level_result(gene.mapping.res, gene.annots)
head(gene.pip.res)
cat(sprintf("%d genes with PIP >= 0.8", 
            length(gene.pip.res$gene_name[gene.pip.res$gene_pip >= 0.8])))
gene.cs <- gene_cs(gene.mapping.res, by.locus = TRUE, gene.cs.coverage = 0.8)
write_tsv(gene.cs,file = "seperated_GBM_UK10K_SNP_fine_mapping_result_gene.cs.txt")
# gene.cs<-vroom(gene.cs_GBM.txt)

tiff('Figure 1.tiff',height = 15000,width = 12000,res= 600,compression = "lzw")
gene_manhattan_plot(gene.pip.res, gene.pip.thresh = 0.8)
dev.off ()
gc()
rm(gene.cs,genomic.annots,gene.mapping.res)