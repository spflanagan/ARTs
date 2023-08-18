library(vcfR)
library(dplyr)
# Run from R directory
setwd("../results/")
source("../R/genomic_summaries.R")

QTLS<-TRUE
SUPS<-TRUE
if(isTRUE(QTLS)){
  path <- "/mnt/BigData/fixedART-results/qtls"
  qtl_outliers<-dplyr::bind_rows(lapply(
    list.files(path=path, 
               pattern = "*.vcf",
               full.names = TRUE),
    get_summary,
    path = path))
  write.csv(qtl_outliers,"qtl_outliers.csv",quote=FALSE,row.names = FALSE)
  
}

if(isTRUE(SUPS)){
  path <- "/mnt/BigData/fixedART-results/supergene/"
  supergene_outliers<-dplyr::bind_rows(lapply(
    list.files(pattern=".*q1\\_.*vcf",
               path= path, 
               full.names = TRUE),
    get_summary,
    path = path))
  write.csv(supergene_outliers,"supergene_outliers_q1.csv",quote=FALSE,row.names = FALSE)
  
}
