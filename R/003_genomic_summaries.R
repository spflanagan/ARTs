# From https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
# answer/code by stas g, posted 5 Aug 2015, accessed 6 Dec 2021
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
  QTLsInPeaks<-sapply(qtl_locs,function(q,peaks){
    n<-sapply(peaks, function (p){
      if(p - m <= q & p + m >= q){
        return(TRUE)
      }  else{
        return(FALSE)
      }
    })
    return(sum(n))
  }, peaks=sigPeaks)
  
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

get_summary<-function(filename, path){
  pattern <- gsub(paste0(path,"(.*)_(pop_\\d).vcf"),"\\1",filename)
  pop <-gsub(paste0(path,"(.*)_(pop_\\d).vcf"),"\\2",filename)
  summarize_peaks(path, pattern=pattern, pop=pop)
}