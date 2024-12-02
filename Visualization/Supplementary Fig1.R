pacman::p_load("vroom","data.table","readr","tidyr","dplyr","devtools","ggplot2","tidyverse","GWAS.utils","CMplot")
setwd("G:/A/target/GWAS")
A<-vroom("ALL_GLIOMA.txt.gz",col_select = c("SNP","CHR","POS","P"))
colnames(A)[colnames(A) == "P"] <- "Pan-Glioma"

B<-vroom("GBM_GLIOMA.txt.gz",col_select = c("SNP","P"))
colnames(B)[colnames(B) == "P"] <- "GBM"

C<-vroom("nonGBM_GLIOMA.txt.gz",col_select = c("SNP","P"))
colnames(C)[colnames(C) == "P"] <- "non-GBM"
A<-left_join(A,B,by = "SNP")
gc()
A<-left_join(A,C,by = "SNP")
colnames(A)<-c("SNP","Chromosome","Position","Pan-Glioma","GBM","non-GBM")
rm(list = setdiff(ls(), c("A")))
gc()
setwd("G:/")
MVP.Report(A, plot.type=c("m","q"), multracks=TRUE, threshold=c(1e-5,5e-8),threshold.lty=c(1,2), 
           threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e3,
           chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),signal.cex=c(1,1),
           file.type="tiff",memo="",dpi=600,file.output = TRUE)
