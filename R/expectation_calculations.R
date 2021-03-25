#' Function to generate predictions based on my life cycle diagram
#' @param morph Specify which morph you want to know about. Options are:
#'              "courter-parent" or "CP"
#'              "courter-sneaker" or "CS"
#'              "noncourter-parent" or "NP"
#'              "noncourter-sneaker" or "NS"
#' @param freqs A vector containing the frequencies of the other morphs (fCP,fCS,fNP,fNS). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @param Nm number of males in the population. Default is 500.
#' @param Nf number of females in the population. Default is 500.
#' @param pf Proportion of males sampled by each female. Default is 0.1.
#' @param r Relative reproductive effort/contribution within each clutch for parents. Default is 2/3.
#' @param c Sperm competition coefficient. Default is 0.5
#' @param ws Sexual selection strength, aka female preference for male type. Default is 1 (unidirectional preference for courters).
#' @param wn Selection strength on nesting trait in males, aka nest survival. Default is 1 (parental male nests survive and non-parental nests all die).
#' @param wv Viability selection against courtship and nesting traits. Default is exp(-0.5/(2*50)).
morph_predictions<-function(
  morph,
  freqs,
  Nm=500,
  Nf=500,
  pf=0.1,
  r=2/3,
  c=0.5,
  ws=1,
  wn=1,
  wv=exp(-0.5/(2*50))
){
  ###### work in progress ######
  # these are not properly integrated
  # possibly need to run each step for each morph
  
  # sanity checks of morph and freqs
  
  
  # probability of finding a nest partner for courters
  # really this is the number of times each male can be sampled
  # should somehow make this the number of nests?
  if(tolower(morph) %in% c("courter-parent","cp")){
    pm <- ws*( (pf*Nm*Nf) / (Nm*freqs[["fCP"]]+Nm*freqs[["fCS"]]) )
  } else if(tolower(morph) %in% c("courter-parent","cp")){
    pm <- ws*( (pf*Nm*Nf) / (Nm*freqs[["fCP"]]+Nm*freqs[["fCS"]]) )
  } else if(tolower(morph) %in% c("noncourter-parent","np")){
    pm <- (1-ws)*( (pf*Nm*Nf) / (Nm*freqs[["fNP"]]+Nm*freqs[["fNS"]]) )
  } else{
    pm <- (1-ws)*( (pf*Nm*Nf) / (Nm*freqs[["fNP"]]+Nm*freqs[["fNS"]]) )
  }
  # this is hacky and mimicking the simulation rather than being properly mathy
  if(pm > 1){
    pm <- 1
  }

  # proportion of eggs that are fertilized
  # this should include the proportion of nests
  if(tolower(morph) %in% c("courter-parent","cp","courter-sneaker","cs")){
    pf <- r*( (freqs[["fCP"]]+freqs[["fCS"]]) + c*(freqs[["fNP"]]+freqs[["fNP"]]) ) 
  } else{
    pf <-  (1-r)*( (freqs[["fCP"]]+freqs[["fCS"]]) + c*(freqs[["fNP"]]+freqs[["fNP"]]) )
  }

  
  # proportion of nests that survive
  pn <- wn*(freqs[["fCP"]] + freqs[["fNP"]]) + 
    (1-wn)*(freqs[["fCS"]] + freqs[["fNS"]])
  
  
  
}

