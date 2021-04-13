library(shiny) 
library(shinydashboard)
library(plotly)
source("morph_predictions.R")

data<-readRDS("../results/expectations_list.RDS")
# remove ones where population will crash
data<-data[complete.cases(data),]


no_rows<-function(data, sliderCP, sliderNP){
  if(nrow(data[as.character(data$CP_freq)==as.character(sliderCP) & 
               as.character(data$NP_freq)==as.character(sliderNP),]) == 0){
    "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
  } else{
    NULL
  }
}

no_plots<-function(data, sliderCP, sliderNP){
  subdat<-data[as.character(data$CP_freq)==as.character(sliderCP) & 
         as.character(data$NP_freq)==as.character(sliderNP),]
  if(sum(subdat$NS_rs)==0){
    "The chosen combination results in 0 noncourter-sneakers, so there are no graphs to show."
  } else{
    NULL
  }
}


# create UI
ui <- dashboardPage(
  dashboardHeader(
    title="Frequency-dependent predictions of alternative reproductive tactics"
  ),
  dashboardSidebar(
    sliderInput("sliderCP","CP freq", min=0, max=1, step=0.05, value=1),
    sliderInput("sliderNP","NP freq", min=0, max=1, step=0.05, value=1)
  ),
  dashboardBody(
    plotlyOutput('contours'),
    
    fluidRow(
      box("Every combination results in noncourter-parents (NP) having a relative reproductive success of 0.")
    )
  ))

# create the plots
server <- function(input, output, session) { 
  

  
  sliders<-reactive({
    validate(no_plots(data, input$sliderCP, input$sliderNP),
             no_rows(data, input$sliderCP, input$sliderNP)
    )
    CP<-get(input$sliderCP)
    NP<-get(input$sliderNP)
  })
  
  output$contours <- renderPlotly({
    sub_calcs<-data[which(
      as.character(data$CP_freq)==as.character(sliders()$CP) & 
        as.character(data$NP_freq)==as.character(sliders()$NP)),]  
    
    # fig 1: NS vs CS
    fig1 <- plot_ly(
      x = sub_calcs$NS_freq, 
      y = sub_calcs$CS_freq, 
      z = as.matrix(sub_calcs[,c("NS_rs","CS_rs")]), 
      type = "contour"
    )
    # add axis labels
    x<-list(title="Noncourter-Sneaker frequency")
    y<-list(title="Courter-Sneaker frequency")
    fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
    # add label to contour names
    fig1 <- fig1 %>% colorbar(title = "Relative RS")
    
    fig1
  })
}

shinyApp(ui, server)
