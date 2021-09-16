#' Function to check frequency inputs
#' @param freqs A vector containing the frequencies of the other morphs (CP,CN,NP,NN). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @return Returns the checked/corrected vector of allele frequencies
check_freqs<-function(freqs){
  # sanity checks of freqs
  if(length(freqs)<4){
    if(length(freqs)==3){
      if(length(which((names(freqs) %in% c("CP","CN","NP"))==TRUE))==3){
        freqs["NN"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CP","CN","NN"))==TRUE))==3){
        freqs["NP"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CP","NP","NN"))==TRUE))==3){
        freqs["CN"]<-1-sum(freqs)
      } else if (length(which((names(freqs) %in% c("CN","NP","NN"))==TRUE))==3){
        freqs["CP"]<-1-sum(freqs)
      }
    }else{
      stop("Too few frequencies were provided to the function.")
    }
  } else{
    if(length(which((names(freqs) %in% c("CP","CN","NP","NN"))==TRUE))<4){
      names(freqs)<- c("CP","CN","NP","NN")
      warning("The given frequencies have been renamed in the order CP, CN, NP, and NN.")
    }  
  }
  if(round(sum(freqs,na.rm=TRUE)) != 1){
    stop("The sum of the frequencies is not 1.")
  }
  return(freqs)
}
