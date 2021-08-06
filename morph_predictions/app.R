library(shiny) 
library(shinydashboard)
library(rsconnect)
library(plotly)
library(dplyr)
library(vegan)
library(tidyr)
library(rmarkdown)
source("morph_predictions.R")
source("morph_gens.R")
source("check_freqs.R")


get_freqs<-function(){
  # get all of the possible frequency combos
  freqs_list<-read.table("freqs_list.txt",sep='\t',header = TRUE)
  return(freqs_list)
}

get_results<-function(){
  morph_results<-readRDS("morph_results.RDS")
  morph_results$diversity<-vegan::diversity(round(morph_results[,c("CP","CS","NP","NS")],4))
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
no_rows<-function(data, whichSlider, sliderFreq, chosenR, chosenC){
  if(whichSlider=="sliderCP"){
    if(nrow(data[as.character(data$initial_CP)==as.character(sliderFreq) & 
                 as.character(data$r) == as.character(chosenR) &
                 as.character(data$c) == as.character(chosenC),]) == 0){
      "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
    } else{
      NULL
    }
  }else {
    if(nrow(data[as.character(data$initial_NP)==as.character(sliderFreq) & 
                 as.character(data$r) == as.character(chosenR) &
                 as.character(data$c) == as.character(chosenC),]) == 0){
      "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
    } else{
      NULL
    }
  }

}

# Function to check if the combination results in zero NS, which means there is no plot.
no_plots<-function(data, whichSlider, sliderFreq, chosenR, chosenC){
  if(whichSlider=="sliderCP"){
    subdat<-data[as.character(data$initial_CP)==as.character(sliderFreq) & 
             as.character(data$r) == as.character(chosenR) &
             as.character(data$c) == as.character(chosenC),]
    if(sum(subdat$NS)==0){
      "The chosen combination results in 0 noncourter-sneakers, so there are no graphs to show."
    } else{
      NULL
    }
  } else {
    subdat<-data[as.character(data$initial_NP)==as.character(sliderFreq) & 
                   as.character(data$r) == as.character(chosenR) &
                   as.character(data$c) == as.character(chosenC),]
    if(sum(subdat$NS)==0){
      "The chosen combination results in 0 noncourter-sneakers, so there are no graphs to show."
    } else{
      NULL
    }
  }
}

# Function to generate our subset of data
create_subset<-function(data,whichSlider, sliderFreq,chosenR, chosenC){
  if(whichSlider=="sliderCP"){
    sub<-data[which(
      as.character(data$initial_CP)==as.character(sliderFreq) & 
        as.character(data$r) == as.character(chosenR) &
        as.character(data$c) == as.character(chosenC)),]  
    return(sub)
  }else {
    sub<-data[which(
        as.character(data$initial_NP)==as.character(sliderFreq) &
        as.character(data$r) == as.character(chosenR) &
        as.character(data$c) == as.character(chosenC)),]  
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
                     sliderInput("sliderCP","CP freq", min=0, max=1, step=0.05, value=0.9),
                     sliderInput("sliderNP","NP freq", min=0, max=1, step=0.05, value=0),
                     sliderInput("r","Relative reproductive input (parents)",min=0.11, max=0.211,step=0.01,value=0.11),
                     sliderInput("c","Sperm competition coefficient",min=0, max=1,step=0.25,value=0),
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
                column(12,
                       plotlyOutput('contours'))
              ),
              
              fluidRow(
                column(12, tableOutput('table'))
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
        morph_results<-get_results()
      }
    )
    
  })
  
  # check the inputs and create the subset
  subdatCP<-reactive({
    data<-get_data()
    validate(no_plots(data, "sliderCP",input$sliderCP, input$r, input$c),
             no_rows(data, "sliderCP",input$sliderCP,  input$r, input$c)
    )
    create_subset(data, "sliderCP", input$sliderCP,  input$r, input$c)

  }) 
  subdatNP<-reactive({
    data<-get_data()
    validate(no_plots(data,"sliderNP", input$sliderNP, input$r, input$c),
             no_rows(data, "sliderNP",input$sliderNP,  input$r, input$c)
    )
    create_subset(data, "sliderNP", input$sliderNP,  input$r, input$c)
    
  }) 
  
  
  # create the plot
  output$contours <- renderPlotly({
    
    # fig 1: adjusting the CP bar
   sub_calcs<-subdatCP()
   
   
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
   fig1 <- fig1 %>% layout(xaxis=x,yaxis=y, annotations=list(text="initial CP from slider, initial NP=0",
                                                             x=0.5,y=1,
                                                             showarrow=FALSE)
   )
   # add label to contour names
   fig1 <- fig1 %>% colorbar(title = "Diversity of the population")
   
   
   
    
    # fig 2: adjusting the NP bar
   sub_calcs<-subdatNP()
   
   
   data_wide <- tidyr::spread(
     sub_calcs[,c("initial_CS","initial_NS","diversity")],
     initial_CS,
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
     type = "contour"
   )
   # add axis labels
   x<-list(title="Initial Noncourter-Sneaker frequency")
   y<-list(title="Initial Courter-Sneaker frequency")
   fig2 <- fig2 %>% layout(xaxis=x,yaxis=y, annotations=list(
     text="initial NP from slider, initial CP=0",
     x=0.5,y=1,
     showarrow=FALSE))
   # add label to contour names
   fig2 <- fig2 %>% colorbar(title = "Diversity of the population")
   
   
    fig<-subplot(fig1,fig2,titleX=TRUE,titleY=TRUE,margin=0.1)
    
   fig 
    
  }) 
  
  #output$table<-renderTable(subdat())
}

shinyApp(ui, server)
