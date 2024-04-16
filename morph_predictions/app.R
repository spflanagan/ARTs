library(shiny) 
library(shinydashboard)
library(rsconnect)
library(plotly)
library(dplyr)
library(vegan)
library(tidyr)
library(rmarkdown)
source("morph_gens_ns.R")
source("check_freqs.R")


get_freqs<-function(){
  # get all of the possible frequency combos
  freqs_list<-read.table("freqs_list.txt",sep='\t',header = TRUE)
  return(freqs_list)
}

get_results<-function(survival){
  if(survival=="0%"){
    morph_results<-readRDS("morph_results_Ns_10000.RDS")  
  } else if(survival=="10%"){
    morph_results<-readRDS("morph_results_Ns_10000_ten-survive.RDS")  
  } else {
    morph_results<-readRDS("morph_results_Ns_10000-fifty-survive.RDS")  
  }
  
  morph_results$diversity<-vegan::diversity(round(morph_results[,c("CP","CN","NP","NN")],4))
  return(morph_results)
}

create_predictions<-function(gens,
                             Nm,
                             Nf,
                             r,
                             c,
                             ws,
                             wn,
                             wv,
                             freqs_list){

  
  data<-dplyr::bind_rows(apply(freqs_list,1,
                      morph_gens,
                      gens=gens,
                      Nm=Nm,
                      Nf=Nf,
                      r=r,
                      c=c,
                      ws=ws,
                      wn=wn,
                      wv=wv))
  # remove ones where population will crash
  data<-data[complete.cases(data),]
  return(data)
}


# Function to check if the current combination has any relevant rows
no_rows<-function(data, whichSlider, sliderFreq, chosenR, chosenC, chosenNs){
  if(whichSlider=="sliderCP"){
    if(nrow(data[round(data$initial_CP,2)==round(sliderFreq,2) & 
                 round(data$r,1) == round(chosenR,1) &
                 round(data$c,2) == round(chosenC,2) &
                 data$num_sneak == chosenNs,]) == 0){
      "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
    } else{
      NULL
    }
  }else {
    if(nrow(data[round(data$initial_NP,2)==round(sliderFreq,2) & 
                 round(data$r,1) == round(chosenR,1) &
                 round(data$c,2) == round(chosenC,2) &
                 data$num_sneak == chosenNs,]) == 0){
      "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
    } else{
      NULL
    }
  }

}

# Function to check if the combination results in zero NN, which means there is no plot.
no_plots<-function(data, whichSlider, sliderFreq, chosenR, chosenC, chosenNs){
  if(whichSlider=="sliderCP"){
    subdat<-data[round(data$initial_CP,2)==round(sliderFreq,2) & 
             round(data$r,1) == round(chosenR,1) &
             round(data$c,2) == round(chosenC,2) &
               data$num_sneak == chosenNs,]
    if(sum(subdat$NN)==0){
      "The chosen combination results in 0 noncourter-nonparents, so there are no graphs to show."
    } else{
      NULL
    }
  } else {
    subdat<-data[round(data$initial_NP,2)==round(sliderFreq,2) & 
                   round(data$r,1) == round(chosenR,1) &
                   round(data$c,2) == round(chosenC,2) &
                   data$num_sneak == chosenNs,]
    if(sum(subdat$NN)==0){
      "The chosen combination results in 0 noncourter-nonparents, so there are no graphs to show."
    } else{
      NULL
    }
  }
}

# Function to generate our subset of data
create_subset<-function(data,whichSlider, sliderFreq,chosenR, chosenC, chosenNs){
  if(whichSlider=="sliderCP"){
    sub<-data[which(
      round(data$initial_CP,2)==round(sliderFreq,2) &
        round(data$r,1) == round(chosenR,1) &
        round(data$c,2) == round(chosenC,2) &
        data$num_sneak == chosenNs),] 
    return(sub)
  }else {
    sub<-data[which(
        round(data$initial_NP,2)==round(sliderFreq,2) &
        round(data$r,1) == round(chosenR,1) &
        round(data$c,2) == round(chosenC,2) &
          data$num_sneak == chosenNs),]  
    return(sub)
  }
  
}

# create UI
ui <- dashboardPage(
  dashboardHeader(
    title="Frequency-dependence",
    titleWidth=250
  ),
  dashboardSidebar(
    
    hr(),
    sidebarMenu(id="tabs",
               # menuItem("Parameter Settings", tabName = "params", icon=icon("table")),
                menuItem("Plots", tabName="plot", icon=icon("line-chart"), selected=TRUE),
                menuItem("ReadMe", tabName = "readme", icon=icon("book"))
    ),
    hr(),
    conditionalPanel("input.tabs == 'plot'",
                     selectInput(
                       "survival",
                       "How many of the disfavoured morph survive?",
                       c("0%","10%","50%")
                     )
    ),
    hr(),
    conditionalPanel("input.tabs == 'plot'",
                     sliderInput("sliderCP","CP freq", min=0, max=1, step=0.05, value=0.25,round=2),
                     sliderInput("sliderNP","NP freq", min=0, max=1, step=0.05, value=0.25,round=2),
                     sliderInput("r","Relative reproductive input (courters/non-courters)",min=0, max=2,step=0.1,value=0.5,round=2),
                     sliderInput("c","Sperm competition coefficient",min=0, max=1,step=0.25,value=0.5, round=2),
                     sliderInput("Ns","Maximum number of sneakers",min=1, max=5,step=1,value=3, round=0),
                     width=250
    )
   
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "readme",
              fluidPage(
                includeMarkdown( 'README.md')
              )
      ),
      tabItem(tabName="plot",
              fluidRow(
                column(6,
                       plotlyOutput('contoursCP')),
                column(6,
                       plotlyOutput('contoursNP'))
              ),
              fluidRow( column(12,verbatimTextOutput('info'))),
              fluidRow(
                column(6, tableOutput('tableCP')),
              
                column(6, tableOutput('tableNP'))
              )
      )#,
      # tabItem(tabName = "params",
      #         fluidRow(
      #           column(3,
      #                  numericInput("gens",
      #                               h3("Number of generations"),
      #                               value=100))
      #         ),
      #         fluidRow(
      #           column(3, 
      #                  numericInput("Nm", 
      #                               h3("Number of males"), 
      #                               value = 500)
      #                  ),
      #           column(3, 
      #                  numericInput("Nf", 
      #                               h3("Number of females"), 
      #                               value = 500)
      #           )
      #         ),
      #         fluidRow(
      #           column(3, 
      #                  numericInput("ws", 
      #                               h3("Sexual selection strength"), 
      #                               value = 1)
      #           ),
      #           column(3, 
      #                  numericInput("wn", 
      #                               h3("Nest survival selection"), 
      #                               value = 1)
      #           ),
      #           column(3, 
      #                  numericInput("wv", 
      #                               h3("Viability selection strength"), 
      #                               value = exp(-0.5/(2*50)))
      #           )
      #         )
      # )
    )
  )
 
)

# server to show the plot
server <- function(input, output, session) { 
  
 
  get_data<-reactive({
    
    withProgress(
      message = "Loading... Please wait", {
        freqs_list<-get_freqs()
        morph_results<-get_results(input$survival)
      }
    )
    
  })
  
  # check the inputs and create the subset
  subdatCP<-reactive({
    data<-get_data()
    validate(no_plots(data, "sliderCP",input$sliderCP, input$r, input$c, input$Ns),
             no_rows(data, "sliderCP",input$sliderCP,  input$r, input$c, input$Ns)
    )
    create_subset(data, "sliderCP", input$sliderCP,  input$r, input$c, input$Ns)

  }) 
  subdatNP<-reactive({
    data<-get_data()
    validate(no_plots(data,"sliderNP", input$sliderNP, input$r, input$c, input$Ns),
             no_rows(data, "sliderNP",input$sliderNP,  input$r, input$c, input$Ns)
    )
    create_subset(data, "sliderNP", input$sliderNP,  input$r, input$c, input$Ns)
    
  }) 
  
  
  # create the plot
  output$contoursCP <- renderPlotly({
    
    # fig 1: adjusting the CP bar
   sub_calcsCP<-subdatCP()
   
   
   data_wide <- tidyr::spread(
     sub_calcsCP[,c("initial_CN","initial_NN","diversity")],
     initial_CN,
     diversity
    )
   rownames(data_wide)<-data_wide[,1]
   data_wide<-data_wide[,-1]
   
   # fig 1: diversity with NN vs CN
   fig1<-plot_ly(
     x = as.numeric(colnames(data_wide)), 
     y = as.numeric(rownames(data_wide)), 
     z = as.matrix(data_wide),
     customdata = sub_calcsCP[,c("CP","CN","NP","NN")],
     colorscale=list(seq(0,1,length.out = 9),
                     c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
     type = "contour",
     source="contoursCP"
   )
   # add axis labels
   x<-list(title="Initial Noncourter-Nonparent frequency")
   y<-list(title="Initial Courter-Nonparent frequency")
   fig1 <- fig1 %>% layout(xaxis=x,yaxis=y, annotations=list(text="initial CP from slider, initial NP=0",
                                                             x=0.5,y=1,
                                                             showarrow=FALSE)
   )
   # add label to contour names
   fig1 <- fig1 %>% colorbar(title = "Diversity of the population") %>%
     event_register('plotly_brushed')
    
  })
  
  # fig 2: adjusting the NP bar
  output$contoursNP <- renderPlotly({
   sub_calcsNP<-subdatNP()
   
   
   data_wide <- tidyr::spread(
     sub_calcsNP[,c("initial_CN","initial_NN","diversity")],
     initial_CN,
     diversity
   )
   rownames(data_wide)<-data_wide[,1]
   data_wide<-data_wide[,-1]
   
   fig2<-plot_ly(
     x = as.numeric(colnames(data_wide)), 
     y = as.numeric(rownames(data_wide)), 
     z = as.matrix(data_wide),
     colorscale=list(seq(0,1,length.out = 9),
                     c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
     type = "contour", 
     hoverinfo="none"
   )
   # add axis labels
   x<-list(title="Initial Noncourter-Nonparent frequency")
   y<-list(title="Initial Courter-Nonparent frequency")
   fig2 <- fig2 %>% layout(xaxis=x,yaxis=y, annotations=list(
     text="initial NP from slider, initial CP=0",
     x=0.5,y=1,
     showarrow=FALSE)) 
   # add label to contour names
   fig2 <- fig2 %>% colorbar(title = "Diversity of the population")
   
    
  }) 
  

  test<-observe({
    
    brush<-event_data('plotly_brushed',source="contoursCP",priority="event")
  
    if(is.null(brush)){
      renderPrint("Waiting for brush")
    }
    else{
      output$info<-renderPrint(brush)
    }
  })  
 
  
  
  output$tableCP <- renderTable({
    d<-subdatCP()
    d<-d[round(d[,"initial_CP"],2)==0.25 & 
           round(d[,"initial_CN"],2) == 0.25 & 
           round(d[,"initial_NP"],2)==0.25 & 
           round(d[,"initial_NN"],2)==0.25,]
  })
  output$tableNP<-renderTable({
    d<-subdatNP()
    d<-d[round(d[,"initial_CP"],2)==0.25 & 
           round(d[,"initial_CN"],2) == 0.25 & 
           round(d[,"initial_NP"],2)==0.25 & 
           round(d[,"initial_NN"],2)==0.25,]
  })
}

shinyApp(ui, server)
