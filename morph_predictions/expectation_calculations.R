# code to calculate expectations for all sorts of morph frequencies

# set to directory of this script
source("morph_predictions.R")
source("check_freqs.R")
source("morph_gens.R")

init<-FALSE
testing<-FALSE
create_outputs<-TRUE

if(isTRUE(init)){
  # create all of the intersections
  freqs_list<-expand.grid(CP=seq(0,1,0.05),
                          CS=seq(0,1,0.05),
                          NP=seq(0,1,0.05), 
                          NS=seq(0,1,0.05))
  freqs_list<-freqs_list[rowSums(freqs_list)==1,]
  
  write.table(freqs_list,"freqs_list.txt",col.names=TRUE,row.names=TRUE)
  
  # create some expectations
  expectations_list<-do.call(rbind,apply(freqs_list,1,morph_gens,gens=100)) # this is slow
  
  # saved on 25 May 2021 to save time for defaults.
  saveRDS(expectations_list,"expectations_100gens.RDS")
  
  
  
} else{
  freqs_list<-read.delim("freqs_list.txt",header=TRUE)
}


