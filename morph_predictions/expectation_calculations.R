# code to calculate expectations for all sorts of morph frequencies

library(tidyr)
library(vegan)
# set to directory of this script
source("morph_predictions.R")
source("check_freqs.R")
source("morph_gens_ns.R")

init<-FALSE
testing<-FALSE
create_outputs<-FALSE
create_outputs_Ns<-TRUE
plot_test<-FALSE
step_by_step<-FALSE
plot_generations<-FALSE

if(isTRUE(init)){
  # create all of the intersections
  freqs_list<-expand.grid(CP=seq(0,1,0.05),
                          CN=seq(0,1,0.05),
                          NP=seq(0,1,0.05), 
                          NN=seq(0,1,0.05))
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
    z = as.matrix(sub_calcs[,c("CP_rs","CN_rs")]), 
    type = "contour"
  )
  # add axis labels
  x<-list(title="Courter-Parent frequency")
  y<-list(title="Courter-Nonparent frequency")
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
  rs<-seq(0,1,0.1)
  cs<-seq(0,1,0.25)
  
  morph_results<-as.data.frame(matrix(ncol=10,nrow=0))
  colnames(morph_results)<-c("initial_CP","initial_CN","initial_NP","initial_NN",
                             "CP","CN","NP","NN","r","c")
  
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
  
  saveRDS(morph_results,"morph_results_20210915.RDS")
}

if(isTRUE(create_outputs_Ns)){
  # create a HUGE data.frame with the different outputs, 
  # so that they can be loaded into the shiny app.
  rs<-seq(0,2,0.1)
  cs<-seq(0,1,0.25)
  
  morph_results<-as.data.frame(matrix(ncol=10,nrow=0))
  colnames(morph_results)<-c("initial_CP","initial_CN","initial_NP","initial_NN",
                             "CP","CN","NP","NN","r","c")
  
  
  for(r in rs){
    for(cv in cs){
      outputs<-dplyr::bind_rows(apply(freqs_list,1,morph_gens_ns,gens=100, 
                                      rs=c(r*8,r*8,8,8),cv=cv))
      print(paste("r=",r,"c=",cv))
      to_save<-dplyr::bind_cols(freqs_list,outputs,.name_repair = "minimal")
      colnames(to_save)[1:4]<-paste0("initial_",colnames(to_save)[1:4])
      to_save$r<-r
      to_save$c<-cv
      morph_results<-dplyr::bind_rows(morph_results,to_save)
    }
  }
  
  saveRDS(morph_results,"morph_results_Ns.RDS")
}


if(isTRUE(plot_test)){
  morph_results<-readRDS("morph_results.RDS")
  morph_results$diversity<-vegan::diversity(round(morph_results[,c("CP","CN","NP","NN")],4))
  
  vars<-morph_results[morph_results$CS>0 & morph_results$NS>0,]
  sub_calcs<-morph_results[morph_results$r == 0.11 & morph_results$c == 0 &
                             morph_results$initial_CP==0.05,]
  data_wide <- tidyr::spread(
    sub_calcs[,c("initial_CN","initial_NN","diversity")],
    initial_CS,
    diversity
  )
  rownames(data_wide)<-data_wide[,1]
  data_wide<-data_wide[,-1]
  
  # fig 1: diversity with NN vs CN
  fig1<-plot_ly(
    x = as.numeric(colnames(data_wide)), 
    y = as.numeric(rownames(data_wide)), 
    z = as.matrix(data_wide),
    colorscale=list(seq(0,1,length.out = 9),
                    c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
    type = "contour"
  )
  # add axis labels
  x<-list(title="Initial Noncourter-Nonparent frequency")
  y<-list(title="Initial Courter-Nonparent frequency")
  fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
  # add label to contour names
  fig1 <- fig1 %>% colorbar(title = "Diversity of the population")
  
  fig1
  
}



# walking through an example
if(isTRUE(step_by_step)){
  freqs=c(CP=0.25,CN=0.25,NP=0.25,NN=0.25)
  Nm=500
  Nf=500
  r=1/4
  c=0.5
  ws=1
  wn=1
  wv=exp(-0.5/(2*50))
  
  freqs<-check_freqs(freqs)
  
  out<-t(do.call(rbind,lapply(c("CP","CN","NP","NN"),function(morph){
    
    # assign various weights to add to things
    if(tolower(morph) %in% c("courter-parent","cp")){
      # courters are preferred
      morph_ws<-ws
      # courters can't sneak
      morph_sneak<-0
      # parents have nests that survive
      morph_wn<-1
      # courter-parents get doubly hit by viability selection
      morph_wv<-wv*wv
    } else if(tolower(morph) %in% c("courter-nonparent","cn")){
      morph_ws<-ws
      # courters can't sneak
      morph_sneak<-0
      # non-parents have nests that die
      morph_wn<-0
      morph_wv<-wv
    } else if(tolower(morph) %in% c("noncourter-parent","np")){
      morph_ws<-1-ws
      # noncourters sneak
      morph_sneak<-1
      # parents have nests that survive
      morph_wn<-1
      morph_wv<-wv
    } else{
      morph_ws<-1-ws
      # noncourters sneak
      morph_sneak<-1
      # sneakers have nests that die
      morph_wn<-0
      # no viability selection against noncourter-nonparents
      morph_wv<-1
    }
    
    # proportion of nests
    ## nests per chosen male / total number of nests (which is Nf)
    n<-Nf/((Nm*freqs[["CP"]]+Nm*freqs[["CN"]])*ws + (Nm*freqs[["NP"]]+Nm*freqs[["NN"]])*(1-ws))
    pn<-(Nm*freqs[[morph]]*n*morph_ws)/Nf
    
    # update to 'actual' r, scaled by frequencies
    actR<-(r*(freqs[["CP"]] + freqs[["CN"]]))/(r*(freqs[["CP"]] + freqs[["CN"]]) + (1-r)*(freqs[["NP"]] + freqs[["NP"]])  )
    # proportion of eggs that are fertilized in your nest
    pfn <- pn*actR
    
    # proportion of eggs that are fertilized in other nests
    pfs <- morph_sneak*(1-pfn)*actR
    
    # proportion of offspring that survive nest abandonment in your nest
    psn <- pfn*morph_wn
    # proportion of offspring that survive nest abandonment in sneaker nests
    psf <- pfs*(wn*(freqs[["CP"]] + freqs[["NP"]]) + (1-wn)*(freqs[["CN"]] + freqs[["NN"]]))
    # overall proportion that survive to juveniles
    ps <- psn + psf
    
    # probability of juvenile survival
    pv<-ps*morph_wv
    
    summary<-data.frame(cbind(n,pn,pfn,pfs,psn,psf,ps,pv))
    rownames(summary)<-morph
    
    return(summary)
  })))
}

if(isTRUE(plot_generations)){
  
  gens=100
  freqs=c(0.25,0.25,0.25,0.25)
  
  freqs<-check_freqs(freqs)
  output<-data.frame(matrix(nrow=gens+1,ncol=4))
  colnames(output)<-c("CP","CN","NP","NN")
  output[1,]<-freqs
  for (i in 1:(gens+1)){
    output[i+1,]<-one_gen(freqs=output[i,],rs=c(0.5*8,0.5*8,8,8),c=0.5)
    if(all(is.na(output[i+1,]))){
      #then the population has crashed
      output[(i+1):(gens+1),]<-c(0,0,0,0)
      break
      
    }
  }
  plot(output$NN,type='l',ylim=c(0,1))
  lines(output$NP,type='l',col=2)
  lines(output$CP,type='l',col=3)
  
}
