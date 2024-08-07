# This file extracts phenotype and genotype information from vcf files
# and then conducts GWAS using the GWASpoly library separately for courter and parent traits

# run from ARTs/R/

library(GWASpoly)
source("gwas_preparation.R")

vcf_files<-c(
  #list.files(path="./single_locus/",pattern="vcf",full.names = TRUE),
  list.files(path="../fixedART-results/qtls",pattern="vcf",full.names = TRUE),
  list.files(path="../fixedART-results/supergene",pattern="vcf",full.names = TRUE)
)

vcf_files<-vcf_files[212:length(vcf_files)]


sum_out<-lapply(vcf_files, function(vcfname) try({
  
  # process files
  genofile <- gsub("\\.vcf","_gt.csv",vcfname)
  trait_name<-gsub("_pop_\\d+\\.vcf","_traits\\.txt",vcfname)
  pop<-as.numeric(gsub("^.*pop_(\\d+).*$","\\1",vcfname))

  qtlfile<-gsub("_traits.txt","_qtlinfo.txt",trait_name)
  qtls<-read.delim(qtlfile)
  VCF2dosage(vcfname, 
             dosage.file = gsub("\\.vcf","_gt.csv",vcfname), 
             geno.code = "GT", 
             ploidy = 2, 
             samples=NULL,
             min.DP=1, 
             max.missing=0.9, 
             min.minor=1)
  phenofile<-get_phenos(vcfname)
  
  # data analysis thing
  data <- read.GWASpoly(ploidy=2, 
                        pheno.file=phenofile, 
                        geno.file=genofile,
                        format="numeric", 
                        n.traits=2, 
                        delim=",")
  
  # specify some parameter settings
  N<-nrow(data@pheno)
  params <- set.params(geno.freq = 1 - 5/N, fixed = "Sex", fixed.type = "factor")
  
  # check for population structure
  data.loco <- set.K(data,LOCO=TRUE,n.core=2)
  data.original <- set.K(data,LOCO=FALSE,n.core=2)
  
  fit_courter<-tryCatch({
    scan_courter(data.loco, qtls, params)
  },  error = function(e) {
    scan_courter(data.original, qtls, params)
  }, finally = {
    data.frame(Marker=NA,
               Chrom=NA,
               Position=NA,
               Model=NA,
               R2=NA,
               pval=NA,
               NearQTL=NA)
  })
  fit_courter$Trait<-"Courter"
  
  fit_parent<-tryCatch({
    scan_parent(data.loco, qtls, params)
  },  error = function(e) {
    scan_parent(data.original, qtls, params)
  }, finally = {
    data.frame(Marker=NA,
               Chrom=NA,
               Position=NA,
               Model=NA,
               R2=NA,
               pval=NA,
               NearQTL=NA)
  })
  fit_parent$Trait<-"Parent"
  
  fit_results<-dplyr::bind_rows(fit_courter,fit_parent)
  fit_results$vcf<-vcfname
  
  write.table(fit_results, 
              "../results/GWAS_summary.csv",
              append = TRUE,
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              sep=",")
  
  return(fit_results)
  
  
}))
