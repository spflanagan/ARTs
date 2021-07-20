# code to calculate expectations for all sorts of morph frequencies

# set to directory of this script
source("morph_predictions.R")
source("check_freqs.R")
source("morph_gens.R")

init<-TRUE
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


if(isTRUE(testing)){
  
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
  
  expectations_list<-dplyr::bind_rows(apply(freqs_list,1,morph_gens,gens=100, r=0.11,c=1))
  
  for(r in rs){
    for(c in cs){
      outputs<-dplyr::bind_rows(apply(freqs_list,1,morph_gens,gens=100, r=r,c=c))
      to_save<-dplyr::bind_cols(freqs_list,outputs,.name_repair = "minimal")
      colnames(to_save)[1:4]<-paste0("initial_",colnames(to_save)[1:4])
      to_save$r<-r
      to_save$c<-c
      morph_results<-dplyr::bind_rows(morph_results,to_save)
    }
  }
  
  saveRDS(morph_results,"morph_results.RDS")
}



