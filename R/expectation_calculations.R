#' Function to generate predictions based on my life cycle diagram
#' @param freqs A vector containing the frequencies of the other morphs (CP,CS,NP,NS). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @param Nm number of males in the population. Default is 500.
#' @param Nf number of females in the population. Default is 500.
#' @param r Relative reproductive effort/contribution within each clutch for parents. Default is 2/3.
#' @param c Sperm competition coefficient. Default is 0.5
#' @param ws Sexual selection strength, aka female preference for male type. Default is 1 (unidirectional preference for courters).
#' @param wn Selection strength on nesting trait in males, aka nest survival. Default is 1 (parental male nests survive and non-parental nests all die).
#' @param wv Viability selection against courtship and nesting traits. Default is exp(-0.5/(2*50)).
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
  
  # sanity checks of freqs
  if(length(freqs)<4){
    if(length(freqs)==3){
      if(length(which((names(freqs) %in% c("CP","CS","NP"))==TRUE))==3){
        freqs["NS"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CP","CS","NS"))==TRUE))==3){
        freqs["NP"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CP","NP","NS"))==TRUE))==3){
        freqs["CS"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CS","NP","NS"))==TRUE))==3){
        freqs["CP"]<-1-sum(freqs)
      }
    }else{
      stop("Too few frequencies were provided to the function.")
    }
  } else{
    if(length(which((names(freqs) %in% c("CP","CS","NP","NS"))==TRUE))<4){
      names(freqs)<- c("CP","CS","NP","NS")
      warning("The given frequencies have been renamed in the order CP, CS, NP, and NS.")
    }  
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

