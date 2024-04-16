# This script creates a list of all combinatorial frequencies
# It is called by run_morph_gens.R and app.R

# create all of the intersections
freqs_list<-expand.grid(CP=seq(0,1,0.05),
                        CN=seq(0,1,0.05),
                        NP=seq(0,1,0.05), 
                        NN=seq(0,1,0.05))
freqs_list<-freqs_list[rowSums(freqs_list)==1,]

write.table(freqs_list,"freqs_list.txt",col.names=TRUE,row.names=FALSE,sep='\t')