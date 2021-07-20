#' Function to check frequency inputs
#' @param freqs A vector containing the frequencies of the other morphs (CP,CS,NP,NS). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @return Returns the checked/corrected vector of allele frequencies
check_freqs<-function(freqs){
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
  if(round(sum(freqs,na.rm=TRUE)) != 1){
    stop("The sum of the frequencies is not 1.")
  }
  return(freqs)
}
