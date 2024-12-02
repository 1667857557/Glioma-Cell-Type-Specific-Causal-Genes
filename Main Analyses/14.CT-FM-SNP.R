pacman::p_load(dplyr, mapgen, susieR, GWAS.utils, CARMA, data.table, vroom, readr)
### CARMA----
# export LD_PRELOAD=/mnt/wslg/distro/root/miniconda3/lib/libmkl_rt.so
# conda install -c miniconda3 mkl
cd /mnt/g/CARMA
R
tempdir()
tempfile()
tempdir <- function() "/mnt/g/rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
rm (list=ls ())
A<-vroom("GBM_gwas.sumstats.tsv.gz")
colnames(A)
gwas.sumstats<-A

if(max(gwas.sumstats$pval) <= 1){
  gwas.sumstats <- gwas.sumstats %>% dplyr::mutate(pval = -log10(pval))
}
sig.loci <- gwas.sumstats %>% dplyr::filter(pval > -log10(1e-5)) %>% dplyr::pull(locus) %>% unique()
files <- list.files(path = "/mnt/g/CARMA/BC_genome_wide", pattern = "\\.txt$", full.names = TRUE)
root_names <- basename(files)
root_names <- sub("^CARMA_", "", root_names)
root_names <- sub(".txt", "", root_names)
unique_root_names <- unique(root_names)
print(unique_root_names)
sig.loci <- sig.loci[!sig.loci %in% unique_root_names]
sig.loci
for (value in sig.loci) {
  locus <- value
  gwas.sumstats.locus <- gwas.sumstats[gwas.sumstats$locus == locus, ]
  LD_blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
  LDREF <- load_UKBB_LDREF(LD_blocks, 
                           locus = locus, 
                           LDREF.dir = "/mnt/g/Database/fine_mapping",
                           prefix = "ukb_b37_0.1")
  matched.sumstat.LD <- match_gwas_LDREF(gwas.sumstats.locus, LDREF$R, LDREF$var_info)
  n=24381
  condz <- LD_diagnosis_susie_rss(matched.sumstat.LD$sumstats$zscore, R = matched.sumstat.LD$R, n = n)
  
  detected_index <- which(condz$conditional_dist$logLR > 2 & abs(condz$conditional_dist$z) > 2)
  cat(sprintf("Detected %d variants with possible allele flipping", length(detected_index)), "\n")
  condz$conditional_dist$detected <- 0
  condz$conditional_dist$detected[detected_index] <- 1
  if (length(detected_index) > 0) {
    gwas.sumstats.locus <- gwas.sumstats.locus[-detected_index, ]
    matched.sumstat.LD <- match_gwas_LDREF(gwas.sumstats.locus, LDREF$R, LDREF$var_info)
  }
  
  rm(LDREF,gwas.sumstats.locus,condz)
  gc()
  sumstats <- matched.sumstat.LD$sumstats
  ld.list <- list(matched.sumstat.LD$R)
  sumstats<-select(sumstats,c(chr,pos,a0,a1,snp,zscore,pval))
  colnames(sumstats)<-c("CHR","POS","Ref","Alt","SNP","Z","Pval")
  z.list<-list()
  lambda.list<-list()
  z.list[[1]]<-sumstats$Z
  lambda.list[[1]]<-1
  CARMA.results_no_annot<-CARMA(z.list,ld.list,lambda.list = lambda.list,
                                outlier.switch=T)
  
  sumstats = sumstats %>% mutate(PIP = CARMA.results_no_annot[[1]]$PIPs, CS = 0)
  
  tryCatch({
    sumstats$CS[CARMA.results_no_annot[[1]]$`Credible set`[[2]][[1]]] = 1
    
  }, error = function(e) {
    cat("在迭代中出现错误", ":", conditionMessage(e), "\n")
  })
  output_file <- paste0("/mnt/g/CARMA/",sub("", "CARMA_GBM_EGFR_", value), ".txt")
  write_tsv(x = sumstats, file = output_file)
  output_file <- paste0("/mnt/g/CARMA/",sub("", "CARMA_GBM_EGFR_", value), ".rds.gz")
  saveRDS(CARMA.results_no_annot, file = output_file)
  print(locus)
  rm(list = setdiff(ls(), c("gwas.sumstats","sig.loci")))
  gc()
} 

setwd("G:/CARMA/GBM_genome_wide")

data_files_rds <- list.files(pattern = "CARMA_\\d+\\.rds\\.gz")
data_files_txt <- list.files(pattern = "CARMA_\\d+\\.txt")
results_df <- data.frame()
get_number <- function(file_name) {
  return(sub(".*_(\\d+)\\..*", "\\1", file_name))
}
for (file_rds in data_files_rds) {
  file_number <- get_number(file_rds)
  file_txt <- paste0("CARMA_", file_number, ".txt")
  
  if (!(file_txt %in% data_files_txt)) {
    next
  }
  CARMA.results_no_annot <- readRDS(file_rds)
  sumstats <- vroom(file_txt)
  sumstats <- sumstats %>%
    mutate(PIP = CARMA.results_no_annot[[1]]$PIPs, CS = 0)
  tryCatch({
    sumstats$CS[CARMA.results_no_annot[[1]]$`Credible model`[[3]]] <- 1
  }, error = function(e) {
    return(NULL)  
  })  
  print(head(sumstats))
  results_df <-rbind(results_df,sumstats)
  rm(sumstats)
  gc()
  # write.csv(sumstats, file=paste0("Processed_", file_txt), row.names = FALSE)
}
sig.loci <- gwas.sumstats %>% dplyr::filter(pval > -log10(5e-8)) %>% dplyr::pull(locus) %>% unique()
results_df$Pval>-log10(1e-5)
results_df<- subset(results_df, Pval>-log10(1e-5))
results_df<- subset(results_df, CS==1)
results_df<- subset(results_df, PIP<0.1)
colnames(results_df)
results_df$POS1<-results_df$POS+1
results_df<-select(results_df,c(CHR,POS,POS1,SNP))
results_df$CHR=paste0('chr',results_df$CHR)
write_tsv(results_df,file = "CARMA_GBM_low_PIP.bed")
### CT-FM-SNP----
source activate ldsc
cd /mnt/g/linux/CTFM
###step1
perl scripts/1_create_sumstats.pl --filein ./sumstats_in/GBM_GLIOMA.txt \
--fileout ./sumstats_ready/GBM_GLIOMA.sumstats.gz \
--chr 2 \
--pos 3 \
--beta 6 \
--se 7 \
--A1 4 \
--A2 5 \
--hg 19 \
--pval 8 \
--N 9

###step2
bash scripts/2_launch_SLDSC_default.sh GBM_TCSC
###step3
Rscript scripts/3_launch_susie_default.R /mnt/g/linux/CTFM GBM_TCSC
conda deactivate
###step4
script_path="scripts/4_CTFMSNP_default_annots.sh"
input_file="/mnt/g/linux/CTFM/CARMA_GBM_low_PIP.bed"
output_dir="/mnt/g/linux/CTFM/out/"
max_jobs=15
counter=0
run_parallel() {
  for ((i = 1; i <= $max_jobs; i++)); do
  $script_path $input_file $output_dir &
    counter=$((counter + 1))
    if [ $counter -ge $max_jobs ]; then
    wait
    counter=0
    fi
    done
    wait
}

run_parallel
###step5
output_folder="/mnt/g/linux/CTFM/out/"
rscript_path="scripts/5_CTFMSNP_default_susie.R"
fixed_param1="/mnt/g/linux/CTFM/"
fixed_param2="GBM_TCSC"

for file in "$output_folder"*.out; do
if [[ -e $file ]]; then
echo "Processing $file"
Rscript "$rscript_path" "$fixed_param1" "$fixed_param2" "$file"
else
  echo "No .out files found in $output_folder"
fi
done

### Process Data----
setwd("G:/Linux/CTFM")
sumstats<-vroom("CARMA_782.txt")
CARMA.results_no_annot<-readRDS("CARMA_782.rds.gz")
sumstats$CS[CARMA.results_no_annot[[1]]$`Credible model`[[3]]]= 1
B<-sumstats[sumstats$CS = "1" ]
sumstats <- subset(sumstats, CS == 1)
sumstats$POS1<-sumstats$POS+1
sumstats<-select(sumstats,c(CHR,POS,POS1,SNP))
sumstats$CHR <- paste("chr", sumstats$CHR, sep = "")
write_tsv(sumstats, file = "CARMA_782.bed", col_names = FALSE)
