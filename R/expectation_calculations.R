# code to calculate expectations for all sorts of morph frequencies

# set to directory of this script
source("morph_predictions.R")

# create all of the intersections
freqs_list<-expand.grid(CP=seq(0,1,0.05),
                        CS=seq(0,1,0.05),
                        NP=seq(0,1,0.05), 
                        NS=seq(0,1,0.05))
freqs_list<-freqs_list[rowSums(freqs_list)==1,]

# calculate the expected values for each combination
expectations_list<-do.call(rbind,apply(freqs_list,1,morph_predictions))

saveRDS(expectations_list,"../results/expectations_list.RDS")

# plot the expectated RS values
library(plotly)
