# code to calculate expectations for all sorts of morph frequencies

# set to directory of this script
source("morph_predictions.R")
source("check_freqs.R")
source("morph_gens.R")

init<-FALSE
testing<-FALSE
create_outputs<-TRUE

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

# for NP = 0 and NS =0
sub_calcs<-expectations_list[which(expectations_list$NP_freq==0 & expectations_list$NS_freq==0),]
fig <- plot_ly(
  x = sub_calcs$CP_freq, 
  y = sub_calcs$CS_freq, 
  z = as.matrix(sub_calcs[,c("CP_rs","CS_rs")]), 
  type = "contour"
)
# add axis labels
x<-list(title="Courter-Parent frequency")
y<-list(title="Courter-Sneaker frequency")
fig <- fig %>% layout(xaxis=x,yaxis=y)
# add label to contour names
fig <- fig %>% colorbar(title = "Relative RS")
# plot
fig


# looping

output<-data.frame(matrix(nrow=11,ncol=4))
output[1,]<-c(CP=0.25,CS=0.25,NP=0.25,NS=0.25)
for (i in 1:11){
  output[i+1,]<-morph_predictions(output[i,])
}


