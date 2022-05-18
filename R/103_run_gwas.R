
# run from ARTs/R/

library(GWASpoly)


vcf_files<-c(
  #list.files(path="./single_locus/",pattern="vcf",full.names = TRUE),
  list.files(path="../fixedART-results/qtls",pattern="vcf",full.names = TRUE),
  list.files(path="../fixedART-results/supergene",pattern="vcf",full.names = TRUE)
)

vcf_files<-vcf_files[212:length(vcf_files)]

get_phenos<-function(vcfname){
  
  pop<-as.numeric(gsub("^.*pop_(\\d+).*$","\\1",vcfname))
  
  trait_name<-gsub("_pop_\\d+\\.vcf","_traits\\.txt",vcfname)
  traits<-read.delim(trait_name)
  
  phenofile <- gsub("traits\\.txt",paste0("pop",pop,"_pheno.csv"),trait_name)
  
  traits <- traits[traits$Gen==12000 & traits$Pop == pop,]
  
  phenos<-traits[,c("Individual","CourtTrait","ParentTrait","PrefTrait","Sex")]
  
  genofile <- gsub("\\.vcf","_gt.csv",vcfname)
  gt<-read.csv(genofile)
  ids<-colnames(gt)[6:length(colnames(gt))]
  phenos<-data.frame(ID=ids,
                     CourtTrait=NA,
                     ParentTrait=NA,
                     Sex=gsub("^(\\w{3}).*$","\\1",ids),
                     Morph=gsub("^(\\w{3})\\d+_(.*)$","\\2",ids))
  
  
  phenos$CourtTrait[phenos$Sex=="MAL"]<-gsub("^.*_(\\w+)$","\\1",
                                             phenos$ID[phenos$Sex=="MAL"])
  phenos$ParentTrait[phenos$CourtTrait %in% c("CP","P")] <-"parent"
  phenos$ParentTrait[phenos$CourtTrait %in% c("C","NON")] <-"parent"
  phenos$ParentTrait <- as.numeric(as.factor(phenos$ParentTrait))
  
  phenos$CourtTrait[phenos$CourtTrait %in% c("CP","C")] <-"courter"
  phenos$CourtTrait[phenos$CourtTrait %in% c("P","NON")] <-"noncourter"
  phenos$CourtTrait <- as.numeric(as.factor(phenos$CourtTrait))
  
  write.csv(phenos,
            phenofile,
            row.names = FALSE)
  return(phenofile)
  
}

scan_courter<-function(data, qtls, params){
   
  # test markers for significance
  court_scan <- GWASpoly(data=data,
                         models=c("additive","1-dom"),
                         traits=c("CourtTrait"),
                         params=params
  )
  
  courtQTLs<-qtls[1,grep("CourterQTL",colnames(qtls))]
  qtllocs<-as.numeric(gsub("(\\d)\\.(\\d+)","\\1",courtQTLs))*1000 + 
    as.numeric(gsub("(\\d)\\.(\\d+)","\\2",courtQTLs))
  
  data2 <- set.threshold(court_scan,method="M.eff",level=0.05)
  qtl <- get.QTL(data=data2,traits="CourtTrait",models="additive",bp.window=10)
  if(nrow(qtl)>0){
    fit.ans <- fit.QTL(data=data2,trait="CourtTrait",
                       qtl=qtl[,c("Marker","Model")])
    
    fit.ans$NearQTL<-0
    for(i in 1:nrow(fit.ans)){
      chromqtl<-courtQTLs[as.numeric(gsub("(\\d)\\.(\\d+)","\\1",courtQTLs)) == fit.ans$Chrom[i]]
      chromlocs<- as.numeric(gsub("(\\d)\\.(\\d+)","\\2",chromqtl))
      for(loc in chromlocs){
        
        if(fit.ans$Position[i] +10 >= loc & fit.ans$Position[i] -10 <= loc){
          fit.ans$NearQTL[i] <- fit.ans$NearQTL[i] +1
        }
      }
    }
    
  } else{
    fit.ans<-data.frame(Marker=NA,
               Chrom=NA,
               Position=NA,
               Model=NA,
               R2=NA,
               pval=NA,
               NearQTL=NA)
  }
 
  return(fit.ans)
  
}

scan_parent<-function(data, qtls, params){
  # test markers for significance
  parent_scan <- GWASpoly(data=data,
                          models=c("additive","1-dom"),
                          traits=c("ParentTrait"),
                          params=params
  )
  
  
  data2 <- set.threshold(parent_scan,method="M.eff",level=0.05)
  
  qtl <- get.QTL(data=data2,traits="ParentTrait",models="additive",bp.window=10)
  if(nrow(qtl) > 0){
    fit.ans <- fit.QTL(data=data2,trait="ParentTrait",
                       qtl=qtl[,c("Marker","Model")])
    
    parentQTLs<-qtls[1,grep("ParentQTL",colnames(qtls))]
    qtllocs<-as.numeric(gsub("(\\d)\\.(\\d+)","\\1",parentQTLs))*1000 + 
      as.numeric(gsub("(\\d)\\.(\\d+)","\\2",parentQTLs))
    
    
    fit.ans$NearQTL<-0
    for(i in 1:nrow(fit.ans)){
      chromqtl<-parentQTLs[as.numeric(gsub("(\\d)\\.(\\d+)","\\1",parentQTLs)) == fit.ans$Chrom[i]]
      chromlocs<- as.numeric(gsub("(\\d)\\.(\\d+)","\\2",chromqtl))
      for(loc in chromlocs){
        
        if(fit.ans$Position[i] +10 >= loc & fit.ans$Position[i] -10 <= loc){
          fit.ans$NearQTL[i] <- fit.ans$NearQTL[i] +1
        }
      }
    }
    
  } else{
    fit.ans<- data.frame(Marker=NA,
                         Chrom=NA,
                         Position=NA,
                         Model=NA,
                         R2=NA,
                         pval=NA,
                         NearQTL=NA)
  }
  
  return(fit.ans)
}

sum_out<-lapply(vcf_files, function(vcfname){
  
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
    NA
  })
  fit_courter$Trait<-"Courter"
  
  fit_parent<-tryCatch({
    scan_parent(data.loco, qtls, params)
  },  error = function(e) {
    scan_parent(data.original, qtls, params)
  }, finally = {
    NA
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
  
  
})
