# This file summarises Tajima's D output from vcftools in the context of 
# the simulation data, specifically comparing averages from regions
# with true QTLs and regions without

if(length(system.file(package='spfTools'))<=1){
  devtools::install_github("https://github.com/spflanagan/spfTools/")
}
library(spfTools)

source("genomic_summaries.R")

td_files<-c(list.files(path="../fixedART-results/supergene/", pattern=".Tajima.D", full.names = TRUE),
            list.files(path="../fixedART-results/qtls/", pattern=".Tajima.D", full.names = TRUE))
td_summary<-dplyr::bind_rows(lapply(td_files, summarize_tajimaD_vcftools))
write.csv(td_summary, "../results/tajimaD_summary.csv",row.names = td_files)
