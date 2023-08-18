# Author: Sarah P. Flanagan

#' A function to extract phenotype information from simulation model output and relate it to vcf files
#' @param vcfname This is the name of the vcf file to use. A corresponding _traits.txt file must exist in the same path
#' @return The name of the phenotype file that has been saved as a csv
#' @export
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

#' A function to run a genome-wide association scan with the courter trait
#' @param data A data output from GWASpoly's set.K function
#' @param qtls The list of true QTLs from the simulation output (read in from qtlinfo.txt)
#' @param params The output of GWASpoly's set.params function
#' @return A data.frame with statistics about the GAWS outputs for the marker
scan_courter<-function(data, qtls, params){
  
  # test markers for significance
  court_scan <- GWASpoly::GWASpoly(data=data,
                         models=c("additive","1-dom"),
                         traits=c("CourtTrait"),
                         params=params
  )
  
  courtQTLs<-qtls[1,grep("CourterQTL",colnames(qtls))]
  qtllocs<-as.numeric(gsub("(\\d)\\.(\\d+)","\\1",courtQTLs))*1000 + 
    as.numeric(gsub("(\\d)\\.(\\d+)","\\2",courtQTLs))
  
  data2 <- GWASpoly::set.threshold(court_scan,method="M.eff",level=0.05)
  qtl <- GWASpoly::get.QTL(data=data2,traits="CourtTrait",models="additive",bp.window=10)
  if(nrow(qtl)>0){
    fit.ans <- GWASpoly::fit.QTL(data=data2,trait="CourtTrait",
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

#' A function to run a genome-wide association scan with the parent trait
#' @param data A data output from GWASpoly's set.K function
#' @param qtls The list of true QTLs from the simulation output (read in from qtlinfo.txt)
#' @param params The output of GWASpoly's set.params function
#' @return A data.frame with statistics about the GAWS outputs for the marker
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
