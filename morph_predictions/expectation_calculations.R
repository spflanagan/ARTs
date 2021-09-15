# code to calculate expectations for all sorts of morph frequencies

library(tidyr)
library(vegan)
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
  
  write.table(freqs_list,"freqs_list.txt",col.names=TRUE,row.names=FALSE,sep='\t')
  
} else{
  freqs_list<-read.delim("freqs_list.txt",header=TRUE, sep='\t')
}


if(isTRUE(testing)){
  # create some expectations
  expectations_list<-do.call(rbind,apply(freqs_list,1,morph_gens,gens=100)) # this is slow
  
  # saved on 25 May 2021 to save time for defaults.
  saveRDS(expectations_list,"expectations_100gens.RDS")
  # calculate the expected values for each combination - one generation
  expectations_list<-do.call(rbind,apply(freqs_list,1,morph_predictions))
  
  
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
  
  # testing the case where the pop goes extinct.
  morph_gens(gens=100,freqs=freqs_list[21,])
  morph_gens(gens=100,freqs=c(0.25,0.25,0.25,0.25))
}


if(isTRUE(create_outputs)){
  # create a HUGE data.frame with the different outputs, 
  # so that they can be loaded into the shiny app.
  rs<-seq(0.11,2.11,0.1)
  cs<-seq(0,1,0.25)
  
  morph_results<-as.data.frame(matrix(ncol=10,nrow=0))
  colnames(morph_results)<-c("initial_CP","initial_CS","initial_NP","initial_NS",
                             "CP","CS","NP","NS","r","c")
  
  #expectations_list<-dplyr::bind_rows(apply(freqs_list,1,morph_gens,gens=100, r=0.11,c=1))
  
  for(r in rs){
    for(cv in cs){
      outputs<-dplyr::bind_rows(apply(freqs_list,1,morph_gens,gens=100, r=r,c=cv))
      print(paste("r=",r,"c=",cv))
      to_save<-dplyr::bind_cols(freqs_list,outputs,.name_repair = "minimal")
      colnames(to_save)[1:4]<-paste0("initial_",colnames(to_save)[1:4])
      to_save$r<-r
      to_save$c<-cv
      morph_results<-dplyr::bind_rows(morph_results,to_save)
    }
  }
  
  saveRDS(morph_results,"morph_results.RDS")
}


morph_results<-readRDS("morph_results.RDS")
morph_results$diversity<-vegan::diversity(round(morph_results[,c("CP","CS","NP","NS")],4))

vars<-morph_results[morph_results$CS>0 & morph_results$NS>0,]
sub_calcs<-morph_results[morph_results$r == 0.11 & morph_results$c == 0 &
                  morph_results$initial_CP==0.05,]
data_wide <- tidyr::spread(
  sub_calcs[,c("initial_CS","initial_NS","diversity")],
  initial_CS,
  diversity
)
rownames(data_wide)<-data_wide[,1]
data_wide<-data_wide[,-1]

# fig 1: diversity with NS vs CS
fig1<-plot_ly(
  x = as.numeric(colnames(data_wide)), 
  y = as.numeric(rownames(data_wide)), 
  z = as.matrix(data_wide),
  colorscale=list(seq(0,1,length.out = 9),
                  c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
  type = "contour"
)
# add axis labels
x<-list(title="Initial Noncourter-Sneaker frequency")
y<-list(title="Initial Courter-Sneaker frequency")
fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
# add label to contour names
fig1 <- fig1 %>% colorbar(title = "Diversity of the population")

fig1



