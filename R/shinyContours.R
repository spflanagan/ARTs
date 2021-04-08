library(shiny) 
library(shinydashboard)
library(plotly)
source("morph_predictions.R")

data<-readRDS("../results/expectations_list.RDS")

ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(sliderInput("sliderCP","CP freq", min=0, max=1, step=0.05, value=1),
                   sliderInput("sliderNS","NS freq", min=0, max=1, step=0.05, value=1)),
  dashboardBody(
    plotlyOutput('contours')
  ))

server <- function(input, output, session) { 
  output$contours <- renderPlotly({
    sub_calcs<-expectations_list[which(expectations_list$CP_freq==input$sliderCP & expectations_list$NS_freq==input$sliderNS),]
    
    # fig 1: NP vs CS
    fig1 <- plot_ly(
      x = sub_calcs$NP_freq, 
      y = sub_calcs$CS_freq, 
      z = as.matrix(sub_calcs[,c("NP_rs","CS_rs")]), 
      type = "contour"
    )
    # add axis labels
    x<-list(title="Noncourter-Parent frequency")
    y<-list(title="Courter-Sneaker frequency")
    fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
    # add label to contour names
    fig1 <- fig1 %>% colorbar(title = "Relative RS")
    
    # fig 2: NP vs NS
    fig2 <- plot_ly(
      x = sub_calcs$NP_freq, 
      y = sub_calcs$NS_freq, 
      z = as.matrix(sub_calcs[,c("NP_rs","NS_rs")]), 
      type = "contour"
    )
    # add axis labels
    x<-list(title="Noncourter-Parent frequency")
    y<-list(title="Noncourter-Sneaker frequency")
    fig2 <- fig2 %>% layout(xaxis=x,yaxis=y)
    # add label to contour names
    fig2 <- fig2 %>% colorbar(title = "Relative RS")
    
    # overall figure
    fig <- subplot(fig1,fig2)
    
    fig
  })
}

shinyApp(ui, server)
