#' Function to generate predictions based on my life cycle diagram
#' @param freqs A vector containing the frequencies of the other morphs (CP,CS,NP,NS). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @param Nm number of males in the population. Default is 500.
#' @param Nf number of females in the population. Default is 500.
#' @param r Relative reproductive effort/contribution within each clutch for courters. Must be between 0 and 1. Default is 1/3.
#' @param c Sperm competition coefficient. Default is 0.5
#' @param ws Sexual selection strength, aka female preference for male type. Default is 1 (unidirectional preference for courters).
#' @param wn Selection strength on nesting trait in males, aka nest survival. Default is 1 (parental male nests survive and non-parental nests all die).
#' @param wv Viability selection against courtship and nesting traits. Default is exp(-0.5/(2*50)).
#' @return Returns a new vector of allele frequencies
morph_predictions<-function(
  freqs=c(CP=0.25,CS=0.25,NP=0.25,NS=0.25),
  Nm=500,
  Nf=500,
  r=2/3,
  c=0.5,
  ws=1,
  wn=1,
  wv=exp(-0.5/(2*50))
){
  
  # check the frequency inputs with the function
  freqs<-check_freqs(freqs)
  
  # run the calculations for each of the morphs given the inputs
  out<-unlist(lapply(c("CP","CS","NP","NS"),function(morph){
    
    # assign various weights to add to things
    if(tolower(morph) %in% c("courter-parent","cp")){
      # courters are preferred
      morph_ws<-ws
      # parents can't sneak
      morph_sneak<-0
      # parents have nests that survive
      morph_wn<-1
      # courter-parents get doubly hit by viability selection
      morph_wv<-wv*wv
    } else if(tolower(morph) %in% c("courter-sneaker","cs")){
      morph_ws<-ws
      morph_sneak<-1
      # sneakers have nests that die
      morph_wn<-0
      morph_wv<-wv
    } else if(tolower(morph) %in% c("noncourter-parent","np")){
      morph_ws<-1-ws
      # parents can't sneak
      morph_sneak<-0
      # parents have nests that survive
      morph_wn<-1
      morph_wv<-wv
    } else{
      morph_ws<-1-ws
      morph_sneak<-1
      # sneakers have nests that die
      morph_wn<-0
      # no viability selection against noncourter-sneakers
      morph_wv<-1
    }
    
    # proportion of nests
    ## nests per chosen male / total number of nests (which is Nf)
    n<-Nf/((Nm*freqs[["CP"]]+Nm*freqs[["CS"]])*ws + (Nm*freqs[["NP"]]+Nm*freqs[["NS"]])*(1-ws))
    pn<-(Nm*freqs[[morph]]*n*morph_ws)/Nf
    
    # proportion of eggs that are fertilized in your nest
    pfn <- pn*(r*(freqs[["CP"]] + freqs[["CS"]]) + (1-r)*c*(freqs[["NP"]] + freqs[["NP"]])  )
    
    # proportion of eggs that are fertilized in other nests
    pfs <- morph_sneak*((1-r)*( (freqs[["CP"]]+freqs[["CS"]]) + c*(freqs[["NP"]]+freqs[["NS"]]) ))
    
    # proportion of offspring that survive nest abandonment in your nest
    psn <- pfn*morph_wn
    # proportion of offspring that survive nest abandonment in sneaker nests
    psf <- pfs*(wn*(freqs[["CP"]] + freqs[["NP"]]) + (1-wn)*(freqs[["CS"]] + freqs[["NS"]]))
    # overall proportion that survive to juveniles
    ps <- psn + psf
    
    # probability of juvenile survival
    pv<-ps*morph_wv
    
    #summary<-data.frame(cbind(pn,pfn,pfs,psn,psf,ps,pv))
    #rownames(summary)<-morph
    
    return(pv)
  }))
  
  out<-out/sum(out)
  names(out)<-c("CP","CS","NP","NS")
  
  return(out)
}
