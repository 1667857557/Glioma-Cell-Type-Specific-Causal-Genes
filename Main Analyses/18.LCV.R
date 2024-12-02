### LCV----
WeightedRegression <- function(X,y,w=array(1,length(y),1) ){
  beta<-solve(t(X*w)%*%X,t(X*w)%*%y)
  return(beta)
}

WeightedMean <- function(y,w=array(1,c(length(y),1)) ){
  mu<-sum(y%*%w)/sum(w)
  return(mu)
}

EstimateK4 <-function(ell,z.1,z.2,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0,nargout=3){
  
  if(ldsc.intercept==0){
    intercept.1<-1/n.1
    intercept.2<-1/n.2
    temp<-WeightedRegression(ell,z.1^2-intercept.1,weights)
    h2g.1<-temp[1]
    temp<-WeightedRegression(ell,z.2^2-intercept.2,weights)
    h2g.2<-temp[1]
  }
  else{
    keep.snps.1<-z.1^2 <= sig.threshold*mean(z.1^2)
    temp<-WeightedRegression(cbind(ell[keep.snps.1],matrix(1,sum(keep.snps.1))),z.1[keep.snps.1]^2,weights[keep.snps.1]);
    intercept.1<-temp[2];
    
    temp<-WeightedRegression(ell,z.1^2-intercept.1,weights);
    h2g.1<-temp[1];
    
    keep.snps.2<-z.2^2 <= sig.threshold*mean(z.2^2)
    temp<-WeightedRegression(cbind(ell[keep.snps.2],matrix(1,sum(keep.snps.2))),z.2[keep.snps.2]^2,weights[keep.snps.2]);
    intercept.2<-temp[2];
    
    temp<-WeightedRegression(ell,z.2^2-intercept.2,weights);
    h2g.2<-temp[1];  
  }  
  
  if(crosstrait.intercept==0) {
    temp<-WeightedRegression(ell,z.1*z.2-intercept.12,weights)
    rho<-temp/sqrt(h2g.1*h2g.2)
  }
  else{
    keep.snps.12=(z.1^2<sig.threshold*mean(z.1^2))*(z.2^2<sig.threshold*mean(z.2^2))==1
    temp<-WeightedRegression(cbind(ell[keep.snps.12],matrix(1,sum(keep.snps.12))),z.1[keep.snps.12]*z.2[keep.snps.12],weights[keep.snps.12])
    intercept.12<-temp[2]
    temp<-WeightedRegression(ell,z.1*z.2-intercept.12,weights)
    rho<-temp[1]/sqrt(h2g.1*h2g.2)
  }
  
  s.1<-sqrt(WeightedMean(z.1^2,weights)-intercept.1)
  s.2<-sqrt(WeightedMean(z.2^2,weights)-intercept.2)
  nz.1<-z.1/s.1
  nz.2<-z.2/s.2
  
  k41<-WeightedMean(nz.2*nz.1^3-3*nz.1*nz.2*(intercept.1/s.1^2)-3*(nz.1^2-intercept.1/s.1^2)*intercept.12/s.1/s.2,weights)
  k42<-WeightedMean(nz.1*nz.2^3-3*nz.1*nz.2*(intercept.2/s.2^2)-3*(nz.2^2-intercept.2/s.2^2)*intercept.12/s.1/s.2,weights)
  
  argout<-c(rho,k41,k42,intercept.12,s.1,s.2,intercept.1,intercept.2)
  return(argout[1:nargout])
  
}
RunLCV <- function(ell,z.1,z.2,no.blocks=100,crosstrait.intercept=1,ldsc.intercept=1,weights=1/pmax(1,ell),sig.threshold=.Machine$integer.max,n.1=1,n.2=1,intercept.12=0){
  mm=length(ell)
  if(length(z.1)!=mm||length(z.2)!=mm){
    stop('LD scores and summary statistics should have the same length')
  }
  
  source("MomentFunctions.R")
  grid<- (-100:100)/100;
  
  size.blocks=floor(mm/no.blocks)
  jackknife=matrix(0,no.blocks,8)
  for(jk in 1:no.blocks){
    if(jk==1)
    {ind<-(size.blocks+1):mm}
    else if(jk==no.blocks)
    {ind <- 1:(size.blocks*(jk-1))}
    else
    {ind<-c(1:((jk-1)*size.blocks), (jk*size.blocks+1):mm)}
    jackknife[jk,] <- EstimateK4(ell[ind],z.1[ind],z.2[ind],crosstrait.intercept,ldsc.intercept,weights[ind],sig.threshold,n.1,n.2,intercept.12,8)
    if(any(is.nan(jackknife))){
      stop('NaNs produced, probably due to negative heritability estimates. Check that summary statistics and LD scores are ordered correctly.')
    }
  }
  rho.est<-mean(jackknife[,1])
  rho.err=sd(jackknife[,1])*sqrt(no.blocks+1)
  flip=sign(rho.est)
  
  jackknife[,2:3]<-jackknife[,2:3]-3*jackknife[,1]
  
  gcp.likelihood=grid;gcp.likelihood[]=0
  for(kk in 1:length(grid)){
    xx<-grid[kk]
    fx<-abs(jackknife[,1])^(-xx)
    numer<-jackknife[,2]/fx-fx*jackknife[,3]
    denom=pmax(1/abs(jackknife[,1]),sqrt(jackknife[,2]^2/fx^2 + jackknife[,3]^2*fx^2 ))
    pct.diff<-numer/denom
    
    gcp.likelihood[kk]<-dt(mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    
    if(xx==-1){
      pval.fullycausal.2<-pt(-flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==1){
      pval.fullycausal.1<-pt(flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1),no.blocks-2)
    }
    if(xx==0){
      zscore<- flip*mean(pct.diff)/sd(pct.diff)/sqrt(no.blocks+1)
    }
    
  }
  
  pval.gcpzero.2tailed=pt(-abs(zscore),no.blocks-1)*2
  
  gcp.pm<-WeightedMean(grid,gcp.likelihood)
  gcp.pse<-sqrt(WeightedMean(grid^2,gcp.likelihood)-gcp.pm^2)
  
  h2.zscore.1<-mean(jackknife[,5])/sd(jackknife[,5])/sqrt(no.blocks+1)
  h2.zscore.2<-mean(jackknife[,6])/sd(jackknife[,6])/sqrt(no.blocks+1)
  
  if(h2.zscore.1<4 || h2.zscore.2<4){
    warning('Very noisy heritability estimates potentially leading to false positives')
  }
  else{
    if(h2.zscore.1<7 || h2.zscore.2<7){
      warning('Borderline noisy heritability estimates potentially leading to false positives')
    }
  }
  if(abs(rho.est/rho.err)<2){
    warning('No significantly nonzero genetic correlation, potentially leading to conservative p-values')
  }
  
  lcv.output<-list(zscore=zscore,pval.gcpzero.2tailed=pval.gcpzero.2tailed,gcp.pm=gcp.pm,gcp.pse=gcp.pse,rho.est=rho.est,rho.err=rho.err,
                   pval.fullycausal=c(pval.fullycausal.1,pval.fullycausal.2),h2.zscore=c(h2.zscore.1,h2.zscore.2))
  
  return(lcv.output)
}
setwd("G:/")
pacman::p_load("vroom", "data.table", "readr", "tidyr", "dplyr", "devtools", "ggplot2", "corrr", "MungeSumstats", "tidyverse", "here", "ldscr")
tempdir()
tempfile()
tempdir <- function() "G:\\rtemp"
unlockBinding("tempdir", baseenv())
utils::assignInNamespace("tempdir", tempdir, ns="base", envir=baseenv())
assign("tempdir", tempdir, baseenv())
lockBinding("tempdir", baseenv())
set.seed(1000)
setwd("G:/A/target/GWAS")
df1 <- fread("GBM_GLIOMA.txt.gz")
colnames(df1)[colnames(df1) == "effect_allele"] <- "A1"
colnames(df1)[colnames(df1) == "other_allele"] <- "A2"
df1<-dplyr::select(df1,c(SNP,A1,A2,Z,N))
gc()

path_ld <- "G:/Linux/ldsc/weights_hm3_no_MHC/" 
setwd("G:/Database/phenotyes/buildGRCh37")
file_list <- list.files(pattern = "\\.tsv.gz$")
results_df <- data.frame()
results_df<-fread("G:/A/phenotypes/LCV_GBM.txt")
file_list <- file_list[!file_list %in% results_df$phenotypes]
for (i in seq_along(file_list)) {
  results <- data.frame()
  
  setwd("G:/Database/phenotyes/buildGRCh37")
  df2 <- vroom(file_list[i],col_select = c("SNP","A1","A2","N","BETA","SE"))
  colnames(df2)[colnames(df2) == "A2"] <- "effect_allele"
  colnames(df2)[colnames(df2) == "A1"] <- "A2"
  colnames(df2)[colnames(df2) == "effect_allele"] <- "A1"
  df2$Z=df2$BETA/df2$SE
  df2<-dplyr::select(df2,c(SNP,A1,A2,Z,N))
  df2 <- df2 %>% drop_na(Z)
  set.seed(1000)
  setwd("G:/")
  source("RunLCV.R")
  source("MomentFunctions.R")
  df3 <- data.frame()
  for (chr in 1:22) {
    file <- paste0(path_ld, chr, ".l2.ldscore.gz")
    sub <- read.table(gzfile(file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df3 <- rbind(df3, sub)
  }
  m <- merge(df3, df1, by = "SNP")
  m2 <- merge(m, df2, by = "SNP")
  data <- m2[order(m2[, "CHR"], m2[, "BP"]), ]
  rm(df2,df3,m,m2)
  gc()
  mismatch <- which(data$A1.x != data$A1.y, arr.ind = TRUE)
  data[mismatch, ]$Z.y <- data[mismatch, ]$Z.y * -1
  data[mismatch, ]$A1.y <- data[mismatch, ]$A1.x
  data[mismatch, ]$A2.y <- data[mismatch, ]$A2.x
  
  tryCatch({
    LCV <- RunLCV(data$L2, data$Z.x, data$Z.y)
    sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)", LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed) / log(10), LCV$rho.est, LCV$rho.err)
    gc()
    
    select <- c(LCV$gcp.pm, LCV$gcp.pse, LCV$pval.gcpzero.2tailed, LCV$rho.est, LCV$rho.err)
    results <- rbind(results, select)
    colnames(results) <- c("gcp.pm","gcp.pse","pval.gcpzero.2tailed","rho.est","rho.err")
    results$phenotypes <- file_list[i]
    
    results_df <- rbind(results_df, results)
    rm(LCV,results)
  }, error = function(e) {
    cat("Error occurred with file:", file_list[i], "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
  rm(data,sub)
  gc()
}
write_tsv(results_df,file = "G:/A/phenotypes/LCV_GBM.txt")

