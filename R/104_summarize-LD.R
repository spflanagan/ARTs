# This file compiles information on the linkage disequilibrium analysis
# and calculates summary statistics for each locus in each analysis

# run from ARTs/results/

ld_files<-c(list.files(path="../fixedART-results/supergene/", pattern=".LD.geno.ld", full.names = TRUE),
            list.files(path="../fixedART-results/qtls/", pattern=".LD.geno.ld", full.names = TRUE))

ld_summary<-dplyr::bind_rows(lapply(ld_files, function(ld_file){
  qtl_file<-gsub("_pop_\\d.LD.geno.ld","_qtlinfo.txt", ld_file)
  pop<-gsub("^.*pop_(\\d).*$","Pop\\1",ld_file)
  
  ld_info<-read.delim(ld_file)
  qtl_info<-read.delim(qtl_file)
  qtl_info<-qtl_info[qtl_info$Pop==pop,-1]
  qtls<-data.frame(
    Chrom=gsub("(\\d)\\.(\\d+)","\\1",qtl_info[!is.na(qtl_info)]),
    Locus=gsub("(\\d)\\.(\\d+)","\\2",qtl_info[!is.na(qtl_info)])
  )
  
  ld_summary<-data.frame(CHR=NA,POS=NA,AvgR2=NA, AvgR2Q=NA, AvgR2NQ=NA)
  counter<-1
  for(chrom in unique(ld_info$CHR)){
    chld<-ld_info[ld_info$CHR==chrom,]
    locs<-qtls[qtls$Chrom ==chrom,"Locus"]
    for(pos in unique(chld$POS1)){
      ld_summary[counter,"CHR"]<-chrom
      ld_summary[counter,"POS"]<-pos
      ld_summary[counter,"AvgR2"]<-mean(chld[chld$POS1==pos,"R.2"], na.rm=TRUE)
      ld_summary[counter,"AvgR2Q"]<-mean(chld[chld$POS1==pos & chld$POS2 %in% locs,"R.2"], na.rm=TRUE)
      ld_summary[counter,"AvgR2NQ"]<-mean(chld[chld$POS1==pos & !(chld$POS2 %in% locs),"R.2"], na.rm=TRUE)
      counter<-counter + 1
    }
  }
  ld_summary$file<-ld_file
  write.table(ld_summary, "ld_summary_all.csv", append=TRUE, quote=FALSE, sep=",", row.names=FALSE, col.names = FALSE)
  return(ld_summary)
}))
