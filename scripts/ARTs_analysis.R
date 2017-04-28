# Analyzing results from ARTs model
## Written by Sarah P. Flangan (spflanagan.phd@gmail.com)
### Date Started: 27 April 2017

setwd("~/Projects/ARTs/results/")

rm.sum<-read.table("random_mating_summary.txt",header=T)
plot(c(1,nrow(rm.sum)),c(0,1),type='n')
for(i in 7:ncol(rm.sum)){ points(rownames(rm.sum),rm.sum[,i],type="l")}

plot(rm.sum$CourterFreq, type="l")

rm.traits<-read.delim("random_mating_traits.txt")
