#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)

#' Find peaks in a series of data
#' Code modified from https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
#' answer/code by stas g, posted 5 Aug 2015, accessed 6 Dec 2021
#' @param x Some series of continuous data
#' @param m A local maximum, defaults to 3
#' @return A vector of locations in the dataset that are within a peak in a series.
#' @export
find_peaks <- function (x, m = 3){
  
  # remove NAs (my modification)
  x<-x[!is.na(x)]
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#' Summarize peaks that have been identified in a series of data
#' @param path The file path to a vcf file
#' @param pattern A naming pattern to use to identify the vcf file
#' @param pop The population ID number that has been appended to the end of the filename before ".vcf"
#' @param m The local maximum to use in the find_peaks() function. Default is 50.
#' @return A data.frame of summary statistics from population genetics comparisons, separated out by genome-wide loci versus QTLs
#' @export
summarize_peaks<-function(path, pattern, pop, m=50){
  
  # get vcf data
  vcfname<-paste0(path,pattern, "_",pop,".vcf")
  vcf1<-vcfR::read.vcfR(vcfname)
  
  # get qtl info
  qtl_info<-read.delim(paste0(path,pattern, "_qtlinfo.txt"))
  qtl_locs<-which(vcf1@fix[,"ID"] %in% qtl_info[1,])
  
  # define male types
  mf_pop<-as.factor(gsub("^(\\w{3}).*$","\\1",colnames(vcf1@gt)[-1]))
  mt_pop<-as.factor(gsub("^(\\w{3}).*_(\\w+)$","\\2",colnames(vcf1@gt)[-1]))
  mt_mal<-factor(mt_pop[mf_pop=="MAL"])
  
  # get differentiation info
  mal_vcf<-vcf1
  mal_vcf@gt<-vcf1@gt[,c(1,which(mf_pop=="MAL")+1)]
  mtm_diff<-genetic_diff(mal_vcf, pops = mt_mal, method = 'nei')
  mtm_diff$Gprimest[is.na(mtm_diff$Gprimest)]<-0
  
  # find the peaks
  peaks<-find_peaks(mtm_diff$Gprimest,m = m)
  Gcutoff<-qnorm(p = 0.99,
                 mean=mean(mtm_diff$Gprimest),
                 sd = sd(mtm_diff$Gprimest))
  sigPeaks<-peaks[which(mtm_diff$Gprimest[peaks] >= Gcutoff)]
  if(is.null(sigPeaks) | length(sigPeaks)==0) { 
    sigPeaks<-0
    QTLsInPeaks <- 0
  }  else {

    QTLsInPeaks<-sapply(qtl_locs,function(q,peaks, m){
      n<-sapply(peaks, function (p){
        if(p - m <= q & p + m >= q){
          return(TRUE)
        }  else{
          return(FALSE)
        }
      })
      return(sum(n))
    }, peaks=sigPeaks, m = m)
  }
  
  
  # summarize data
  simSum<-data.frame(
    sim = vcfname,
    # Heterozygosity
    gHtMean = mean(mtm_diff$Ht),
    gHtMedian = median(mtm_diff$Ht),
    gHtVar = var(mtm_diff$Ht),
    qHtMean = mean(mtm_diff$Ht[qtl_locs]),
    qHtMedian = median(mtm_diff$Ht[qtl_locs]),
    qHtVar = var(mtm_diff$Ht[qtl_locs]),
    # Gprime
    gGpMean = mean(mtm_diff$Gprimest),
    gGpMedian = median(mtm_diff$Gprimest),
    gGpVar = var(mtm_diff$Gprimest),
    qGpMean = mean(mtm_diff$Gprimest[qtl_locs]),
    qGpMedian = median(mtm_diff$Gprimest[qtl_locs]),
    qGpVar = var(mtm_diff$Gprimest[qtl_locs]),
    # peaks
    nPeaks = length(peaks),
    nSigPeaks = length(sigPeaks),
    nQTLS = length(qtl_locs),
    propQTLsInPeaks = length(QTLsInPeaks[QTLsInPeaks > 0])/length(QTLsInPeaks)
  )
  
  return(simSum)
  
}

#' A function to get the summary data from summarize_peaks() given a filename and a path
#' @param filename The base filename to use (without the population or vcf extension)
#' @param path The path where the file is located
#' @return A data.frame that is output from summarize_peaks()
#' @export
get_summary<-function(filename, path){
  print(filename)
  pattern <- gsub(paste0(path,"(.*)_(pop_\\d).vcf"),"\\1",filename)
  pop <-gsub(paste0(path,"(.*)_(pop_\\d).vcf"),"\\2",filename)
  summarize_peaks(path, pattern=pattern, pop=pop)
}


#' A function to summarise the vcftools Tajima's D output, 
#' specifically to identify regions with and without QTLs
#' @param td_name A character name of a Tajima's D output file from vcftools
#' @return A list of summary statistics for all markers, QTLs, and neutral markers
#' @import spfTools
#' @export
summarize_tajimaD_vcftools<-function(td_name){
  library(spfTools)
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
