pacman::p_load("ieugwasr", "TwoSampleMR", "vroom", "data.table", "readr", "tidyr", "dplyr", "devtools", "ggplot2", "MungeSumstats", "tidyverse", "ldscr", "plinkbinr", "MRcML", "doParallel")

tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
### LDSC ----
set.seed(1000)
setwd("G:/Database/phenotyes/buildGRCh37")
setwd("D:/phenotyes/process/buildGRCh37")
file_list <- list.files(pattern = "\\.tsv.gz$")
results_df <- data.frame()
for (i in seq_along(file_list)) {
  setwd("D:/phenotyes/process/buildGRCh37")
  B <- vroom(file_list[i],col_select = c("SNP","A1","A2","N","BETA","SE"))
  colnames(B)[colnames(B) == "A2"] <- "effect_allele"
  colnames(B)[colnames(B) == "A1"] <- "A2"
  colnames(B)[colnames(B) == "effect_allele"] <- "A1"
  B$Z=B$BETA/B$SE
  B <- B %>% drop_na(Z)
  h2_res <- ldsc_h2(munged_sumstats = B, ancestry = "EUR")
  h2_res$id.exposure <- file_list[i]
  set.seed(1000)
  results_df <- rbind(results_df, h2_res)
  rm("B")
  gc()
}
write_tsv(results_df,file = "noGBM_ldsc_h2.txt")

### Extract significant SNPs ----
setwd("G:/phenotyes")
data_files <- list.files(pattern = "\\.gz$")
for (file in data_files) {
  setwd("G:/phenotyes")
  exp_dat<-vroom(file)
  exp_dat<- subset(exp_dat, P<5e-8)
  if (nrow(exp_dat)>0){
    colnames(exp_dat)[colnames(exp_dat) == "SNP"] <-"rsid"
    colnames(exp_dat)[colnames(exp_dat) == "P"] <-"pval"
    tryCatch({
      exp_dat <- ieugwasr::ld_clump(
        dplyr::tibble(exp_dat),clump_r2 = 0.001,
        plink_bin = "G:/plink_Windows.exe",
        bfile = "G:/1000G/1kg.v3/EUR"
      )},
      error = function(e)
      {cat("Error:",conditionMessage(e),"\n")
        cat("Retrying...\n")})
  }
  if (nrow(exp_dat)>0){  
    exp_dat <- TwoSampleMR::format_data(exp_dat,type = "exposure", header = TRUE,snp_col = "rsid",beta_col = "BETA",se_col = "SE",effect_allele_col = "A2",other_allele_col = "A1", pval_col = "pval", chr_col = "CHR", pos_col = "BP",samplesize_col = "N")
    
    exp_dat$exposure <-sub(".tsv.gz", "", file)
    setwd("G:/phenotyes")
    output_file <- paste0(sub(".tsv.gz", "",file), ".tsv")
    write_tsv(exp_dat, output_file)
  }
}
### Forward MR-PheWAS ----
setwd("G:/Database/phenotyes/buildGRCh37/significant")
gc()
file_list <- list.files(pattern = "*.tsv")
A<-fread("G:/Database/phenotyes/buildGRCh37/A2.txt")
file_list<-A$ID
B <- vroom("G:/A/target/GWAS/GBM_GLIOMA.txt.gz")

results_df <- data.frame()
for (i in seq_along(file_list)) {
  setwd("G:/Database/phenotyes/buildGRCh37/significant")
  exp_dat <- fread(file_list[i], head = TRUE, data.table = FALSE)
  exp_dat$id.exposure <- file_list[i]
  if (nrow(exp_dat) > 3) {
    
    outcome_dat <- format_data(B, type = "outcome", snps = exp_dat$SNP,  snp_col = "SNP",
                               beta_col = "BETA",  se_col = "SE",chr_col = "CHR",pos_col = "POS",samplesize_col = "N",
                               effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "P") %>%
      mutate(outcome = "GBM")
    if (nrow(outcome_dat) > 3) {
      set.seed(1000)
      dat3 <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat, action = 3)
      dat3<-subset(dat3, mr_keep != "FALSE")
      if (nrow(dat3) > 3) {
        dat4 <- tsmr_to_rmr_format(dat3)
        if (nrow(dat4) > 3) {
          
          eggrad <- egger_radial(r_input = dat4, alpha = 0.05,weights = 1, summary = TRUE)
          if(!is.null(eggrad[["outliers"]]) && !identical(eggrad[["outliers"]], "No significant outliers")) {dat3 <- dat3[!dat3$SNP %in% eggrad$outliers$SNP, ]
          }
          
        }
        dat4 <- tsmr_to_rmr_format(dat3)
        if (nrow(dat4) > 3) {
          eggrad <- egger_radial(r_input = dat4, alpha = 0.05,weights = 1, summary = TRUE)
          if(!is.null(eggrad[["outliers"]]) && !identical(eggrad[["outliers"]], "No significant outliers")) {dat3 <- dat3[!dat3$SNP %in% eggrad$outliers$SNP, ]
          }
        }
        set.seed(1000)
        results <- mr(dat3, method_list = c("mr_ivw"))
        select=results[1,]
        results
        presso_res = run_mr_presso(dat3,NbDistribution=10000)
        C<- presso_res[[1]][1]$`Main MR results` %>%
          tibble::as_tibble() 
        C<-C[1,]
        colnames(C)
        C<-select(C,c(`Causal Estimate`,Sd,`P-value`))
        C$Global_Pvalue<-presso_res[[1]][["MR-PRESSO results"]] $`Global Test`$Pvalue
        cML_result_DP = mr_cML_DP(dat3$beta.exposure,
                                  dat3$beta.outcome,
                                  dat3$se.outcome,
                                  dat3$se.exposure,
                                  n = dat3$samplesize.exposure,
                                  random_start = 100,
                                  random_start_pert = 100,
                                  random_seed = 100,
                                  num_pert = 200)
        select1 <- data.frame(cML_result_DP$MA_BIC_theta,cML_result_DP$MA_BIC_se,cML_result_DP$MA_BIC_p)
        select1 
        dat3$r2=2*(dat3$beta.exposure)^2/(2*(dat3$beta.exposure)^2+2*(dat3$se.exposure)^2*dat3$samplesize.exposure)
        select$R2xz=sum(dat3$r2)
        dat3$F=(dat3$r2/(1-dat3$r2))*(dat3$samplesize.exposure-2)
        select$F <- paste(
          sprintf("%.1f", min(dat3$F, na.rm = TRUE)), 
          sprintf("%.1f", max(dat3$F, na.rm = TRUE)), 
          sep = "-"
        )        
        N = 24381 
        K = 6191 /18190
        alpha = 0.05
        epower = NA
        A<-results_binary(N, alpha, select$R2xz, K, exp(select$b), epower)
        select$power <- A$Value[A$Parameter == "Power"]
        
        or = exp(select$b)
        lower_ci=exp(select$b-1.96*select$se)
        upper_ci=exp(select$b+1.96*select$se)
        select2$OR <-paste0("OR =",round(or,3)," 95% CI (",round(lower_ci,3)," to ",round(upper_ci,3),")")
        df <- data.frame(select,select1,C,select2)
        results_df <- rbind(results_df, df)

      } else {
        next
      }
    }
    else {
      next
    }
  }
  else {
    next
  }
  rm(list = setdiff(ls(), c("B","results_df","file_list")))
  gc()
  gc()
}
write_tsv(results_df,file = "MR-PheWAS.txt")
### Reverse MR-PheWAS -----
setwd("G:/glioma") 
exp_dat<-vroom("GBM_GLIOMA.tsv")
setwd("G:/phenotyes/buildGRCh37")
data_files <- list.files(pattern = "*.gz")
results_df <- data.frame()
for (file in data_files) {
  outcome_dat <- vroom(file) %>%
    format_data(., type = "outcome", snps = exp_dat$SNP,  snp_col = "SNP",
                beta_col = "BETA",  se_col = "SE",chr_col = "CHR",pos_col = "BP",
                effect_allele_col = "A2", other_allele_col = "A1", pval_col = "P") %>%
    mutate(outcome = "GLIOMA")
  outcome_dat$outcome<-file
  Multiple_dat <- harmonise_data(exp_dat,outcome_dat, action = 3)
  results_df<-rbind(results_df,Multiple_dat)
  rm(Multiple_dat)
  gc()
}
Multiple_dat<-results_df
Multiple_dat<-subset(Multiple_dat, mr_keep != "FALSE")

list_dat <- split(Multiple_dat, Multiple_dat$outcome)

mr_IVW <- function(dat) {
  if (nrow(dat) > 1) {
    set.seed(1000)
    res=TwoSampleMR::mr(
      dat,
      method_list =c("mr_ivw"))
    return(res)
  }
}
mr_wald <- function(dat) {
  if (nrow(dat) < 2) {
    set.seed(1000)
    res=TwoSampleMR::mr(
      dat,
      method_list =c("mr_wald_ratio"))
    return(res)
  }
} 
detectCores(logical = FALSE) 
cl <- makeCluster(7) 1
registerDoParallel(cl) 
clusterExport(cl = cl ,varlist = ls())
start1=Sys.time()
res=parLapply(cl = cl,
              X = list_dat,
              fun= mr_IVW)
end1=Sys.time();print(end1-start1)
res_df <- do.call(rbind, res)
cl <- makeCluster(2)
start1=Sys.time()
registerDoParallel(cl)
res=parLapply(cl = cl,
              X = list_dat,
              fun= mr_wald)

res_df1 <- do.call(rbind, res)
res_df2<-rbind(res_df,res_df1)


setwd("G:/")
res_df2$bonferroni<-p.adjust(res_df2$pval, method = "bonferroni")
res_df2$fdr<-p.adjust(res_df2$pval, method = "fdr")
res_df2 <- arrange(res_df2, fdr)
write_tsv(res_df2, file = "10000,0.001,BID_2SMR_phenotypes_GBM.txt")
gc()
