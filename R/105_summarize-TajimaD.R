
sem<-function(x){
  if(all(is.na(x)) || is.null(x) || length(x)==0) return(x)
  x<-as.numeric(x)
  x<-x[!is.na(x)]
  n<-length(x)
  se<-sd(x)/sqrt(n)
  return(se)
}

summarize_tajimaD_vcftools<-function(td_name){
  td<-read.delim(td_name)
  qtl_file<-gsub("_pop_\\d.Tajima.D","_qtlinfo.txt",td_name)
  pop<-gsub("^.*pop_(\\d).*$","Pop\\1",td_name)
  
  qtl_info<-read.delim(qtl_file)
  qtl_info<-qtl_info[qtl_info$Pop==pop,-1]
  qtls<-data.frame(
    Chrom=gsub("(\\d)\\.(\\d+)","\\1",qtl_info[!is.na(qtl_info)]),
    Locus=gsub("(\\d)\\.(\\d+)","\\2",qtl_info[!is.na(qtl_info)])
  )
  
  td$QTL<-FALSE
  for(chrom in unique(td$CHROM)){
    chtd<-td[td$CHROM==chrom,]
    
    locs<-qtls[qtls$Chrom ==chrom,"Locus"]
    for(i in 1:nrow(chtd)){
      for(loc in locs){
        if(loc >= chtd[i,"BIN_START"] & loc <= chtd[i,"BIN_START"]+chtd[i,"N_SNPS"]){
          # then there's a qtl there
          chtd$QTL[i] <- TRUE
        }
      }
    }
    # Save it back to the overall df
    td[td$CHROM==chrom,"QTL"]<-chtd$QTL
  }
  
  tres<-t.test(td$TajimaD[td$QTL==TRUE], td$TajimaD[td$QTL==FALSE])
  
  return(c(
    avgTD=mean(td$TajimaD),
    semTD=sem(td$TajimaD),
    avgTD_qtls=mean(td$TajimaD[td$QTL==TRUE]),
    semTD_qtls=sem(td$TajimaD[td$QTL==TRUE]),
    avgTD_neutral=mean(td$TajimaD[td$QTL==FALSE]),
    semTD_neutral=sem(td$TajimaD[td$QTL==FALSE]),
    t = tres$statistic,
    df = tres$parameter,
    p = tres$p.value
  ))
}

td_files<-c(list.files(path="../fixedART-results/supergene/", pattern=".Tajima.D", full.names = TRUE),
            list.files(path="../fixedART-results/qtls/", pattern=".Tajima.D", full.names = TRUE))
td_summary<-dplyr::bind_rows(lapply(td_files, summarize_tajimaD_vcftools))
write.csv(td_summary, "../results/tajimaD_summary.csv",row.names = td_files)
