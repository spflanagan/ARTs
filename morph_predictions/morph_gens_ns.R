# object-oriented version of morph_gens
# also now in terms of actual numbers!
source("check_freqs.R")
nest_fertilize<-function(freqs, Nm, max_off,ws=c(1,1,0,0), rs=c(4,4,8,8)){
  # number of nests 
  n<-freqs*ws*Nm
  # number of sperm
  sperm<-freqs*Nm*rs
  # max offspring and the difference
  NnM<-sum(n)*max_off
  diffs<-NnM-sperm
  # any negatives become zero
  diffs[diffs<0] <-0
  # the courters get to use all of theirs
  ferts<-sperm
  ferts[which(ws==1)]<-sperm[which(ws==1)]
  # non-courters get to keep either all of their sperm or as much as can be used
  ferts[which(ws==0)]<-apply(rbind(sperm[which(ws==0)],
                                   diffs[which(ws==0)]),2,min)
  return(ferts)
}



nest_survival<-function(ferts,morph_sneak=c(0,0,1,1),morph_wn=c(1,0,1,0)){
  # get the surviving nests
  nests<-ferts*morph_wn
  nests[which(morph_sneak==1)]<-0
  # modify sneaker survivals based on nest reductions
  nest_loss<-sum(nests)/sum(ferts[which(morph_sneak==0)])
  surv<-ferts*morph_sneak*nest_loss
  # modify non-sneakers
  surv[morph_sneak==0]<-nests[which(morph_sneak==0)]

  return(surv)
}

viability<-function(surv,morph_wvs){
  # reduce based on viability
  via_surv<-surv*morph_wvs
  
  # convert to frequencies
  via_surv<-via_surv/sum(via_surv)
  return(via_surv)
}


one_gen<-function(rs=c(4,4,8,8),
                  freqs=c(CP=0.25,CN=0.25,NP=0.25,NN=0.25),
                  Nm=500,
                  cv=0.5,
                  wn=1,
                  wv=exp(-0.5/(2*50))){
  max_off<-cv*3*rs[4]
  freqs<-check_freqs(freqs)
  eggs<-nest_fertilize(freqs, Nm=Nm,max_off,ws=c(1,1,0,0), rs=rs)
  nests<-nest_survival(eggs,morph_sneak=c(0,0,1,1),morph_wn=c(wn,0,wn,0))
  adults<-viability(nests,morph_wvs=c(wv*wv,wv,wv,1))
  return(adults)
}




#' Function to loop over some number of generations to get the equilibiral morph frequencies
#' @param gens The number of generations
#' @param freqs A vector containing the frequencies of the other morphs (CP,CN,NP,NN). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @return Returns the vector of allele frequencies in the final generation
morph_gens_ns<-function(gens,freqs, ...){
  freqs<-check_freqs(freqs)
  output<-data.frame(matrix(nrow=gens+1,ncol=4))
  colnames(output)<-c("CP","CN","NP","NN")
  output[1,]<-freqs
  for (i in 1:(gens+1)){
    output[i+1,]<-one_gen(freqs=output[i,], ...)
    if(all(is.na(output[i+1,]))){
      #then the population has crashed
      output[(i+1):(gens+1),]<-c(0,0,0,0)
      break
      
    }
  }
  return(output[gens+1,])
}


