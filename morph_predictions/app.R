library(shiny) 
library(shinydashboard)
library(rsconnect)
library(plotly)
library(dplyr)
library(vegan)
library(tidyr)
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
no_plots<-function(data, sliderCP, sliderNP, chosenR, chosenC){
  subdat<-data[as.character(data$initial_CP)==as.character(sliderCP) & 
         as.character(data$initial_NP)==as.character(sliderNP) & 
           as.character(data$r) == as.character(chosenR) &
           as.character(data$c) == as.character(chosenC),]
  if(sum(subdat$NS)==0){
    "The chosen combination results in 0 noncourter-sneakers, so there are no graphs to show."
  } else{
    NULL
  }
}

# Function to generate our subset of data
create_subset<-function(data,sliderCP,sliderNP, chosenR, chosenC){
  sub<-data[which(
    as.character(data$initial_CP)==as.character(sliderCP) & 
      as.character(data$initial_NP)==as.character(sliderNP)) &
      as.character(data$r) == as.character(chosenR) &
      as.character(data$c) == as.character(chosenC),]  
  return(sub)
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
                menuItem("Parameter Settings", tabName = "params", icon=icon("table"), selected=TRUE),
                menuItem("Plots", tabName="plot", icon=icon("line-chart")),
                menuItem("ReadMe", tabName = "readme", icon=icon("book"))
    ),
    hr(),
    conditionalPanel("input.tabs == 'plot'",
                     sliderInput("sliderCP","CP freq", min=0, max=1, step=0.05, value=0.9),
                     sliderInput("sliderNP","NP freq", min=0, max=1, step=0.05, value=0),
                     sliderInput("r","Relative reproductive input (parents)",min=0.11, max=0.211,step=0.1,value=0.11),
                     sliderInput("c","Sperm competition coefficient",min=0, max=1,step=0.25,value=0),
                     width=250
    )
   
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "readme",
              fluidPage(
                tags$iframe(src = './README.md', 
                            width = '100%', height = '800px', 
                            frameborder = 0, scrolling = 'auto'
                )
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
      ),
      tabItem(tabName = "params",
              fluidRow(
                column(3,
                       numericInput("gens",
                                    h3("Number of generations"),
                                    value=100))
              ),
              fluidRow(
                column(3, 
                       numericInput("Nm", 
                                    h3("Number of males"), 
                                    value = 500)
                       ),
                column(3, 
                       numericInput("Nf", 
                                    h3("Number of females"), 
                                    value = 500)
                )
              ),
              fluidRow(
                column(3, 
                       numericInput("ws", 
                                    h3("Sexual selection strength"), 
                                    value = 1)
                ),
                column(3, 
                       numericInput("wn", 
                                    h3("Nest survival selection"), 
                                    value = 1)
                ),
                column(3, 
                       numericInput("wv", 
                                    h3("Viability selection strength"), 
                                    value = exp(-0.5/(2*50)))
                )
              )
      )
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
  subdat<-reactive({
    data<-get_data()
    validate(no_plots(data, input$sliderCP, input$sliderNP, input$r, input$c),
             no_rows(data, input$sliderCP, input$sliderNP, input$r, input$c)
    )
    create_subset(data, input$sliderCP, input$sliderNP, input$r, input$c)

  }) %>%
    bindCache(data)
  
  
  
  # create the plot
  output$contours <- renderPlotly({
    
   sub_calcs<-subdat()
    
    # fig 1: NS vs CS
    fig1 <- plot_ly(
      x = sub_calcs$initial_NS, 
      y = sub_calcs$initial_CS, 
      z = as.matrix(sub_calcs[,c("CS","NS")]), 
      colorscale=list(seq(0,1,length.out = 9),
                      c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
      type = "contour"
    )
    # add axis labels
    x<-list(title="Initial Noncourter-Sneaker frequency")
    y<-list(title="Initial Courter-Sneaker frequency")
    fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
    # add label to contour names
    fig1 <- fig1 %>% colorbar(title = "Frequency at equilibrium")
    
    fig1
    # # fig 1: CP and NP RS
    # fig2 <- plot_ly(
    #   x = sub_calcs$NS_freq, 
    #   y = sub_calcs$CS_freq, 
    #   z = as.matrix(sub_calcs[,c("CP_rs","NP_rs")]), 
    #   type = "contour",
    #   colorscale=list(seq(0,1,length.out = 9),
    #                   c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
    # )
    # # add axis labels
    # x<-list(title="Noncourter-Sneaker frequency (NP RS)")
    # y<-list(title="Courter-Sneaker frequency (CP RS)")
    # fig2 <- fig2 %>% layout(xaxis=x,yaxis=y)
    # # add label to contour names
    # fig2 <- fig2 %>% colorbar(title = "Relative RS")
    # 
    # fig<-subplot(fig1,fig2,titleX=TRUE,titleY=TRUE,margin=0.1)
    
    
    
  }) 
  
  output$table<-renderTable(subdat())
}

shinyApp(ui, server)
