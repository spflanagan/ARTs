#' Function to loop over some number of generations to get the equilibiral morph frequencies
#' @param gens The number of generations
#' @param freqs A vector containing the frequencies of the other morphs (CP,CS,NP,NS). 
#'              The frequencies must be labelled or all four provided in the order above.
#'              If it is a vector with three values, the fourth is calculated as the frequency of the morph of interest.
#' @return Returns the vector of allele frequencies in the final generation
morph_gens<-function(gens,freqs,...){
  freqs<-check_freqs(freqs)
  output<-data.frame(matrix(nrow=gens+1,ncol=4))
  colnames(output)<-c("CP","CN","NP","NN")
  output[1,]<-freqs
  for (i in 1:(gens+1)){
    output[i+1,]<-morph_predictions(freqs=output[i,], ...)
    if(all(is.na(output[i+1,]))){
      #then the population has crashed
      output[(i+1):(gens+1),]<-c(0,0,0,0)
      break
      
    }
  }
  return(output[gens+1,])
}
