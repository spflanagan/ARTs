library(shiny) 
library(shinydashboard)
library(rsconnect)
library(plotly)
source("morph_predictions.R")

create_predictions<-function(Nm,
                             Nf,
                             r,
                             c,
                             ws,
                             wn,
                             wv){
  # create all of the intersections
  freqs_list<-expand.grid(CP=seq(0,1,0.05),
                          CS=seq(0,1,0.05),
                          NP=seq(0,1,0.05), 
                          NS=seq(0,1,0.05))
  freqs_list<-freqs_list[rowSums(freqs_list)==1,]
  data<-do.call(rbind,apply(freqs_list,1,
                      morph_predictions,
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

read_predictions<-function(){
  
  data<-readRDS("expectations_list.RDS")
  # remove ones where population will crash
  data<-data[complete.cases(data),]
  return(data)
}


# Function to check if the current combination has any relevant rows
no_rows<-function(data, sliderCP, sliderNP){
  if(nrow(data[as.character(data$CP_freq)==as.character(sliderCP) & 
               as.character(data$NP_freq)==as.character(sliderNP),]) == 0){
    "The chosen combination does not have any results. Ensure the sum of frequencies is <= 1."
  } else{
    NULL
  }
}

# Function to check if the combination results in zero NS, which means there is no plot.
no_plots<-function(data, sliderCP, sliderNP){
  subdat<-data[as.character(data$CP_freq)==as.character(sliderCP) & 
         as.character(data$NP_freq)==as.character(sliderNP),]
  if(sum(subdat$NS_rs)==0){
    "The chosen combination results in 0 noncourter-sneakers, so there are no graphs to show."
  } else{
    NULL
  }
}

# Function to generate our subset of data
create_subset<-function(data,sliderCP,sliderNP){
  sub<-data[which(
    as.character(data$CP_freq)==as.character(sliderCP) & 
      as.character(data$NP_freq)==as.character(sliderNP)),]  
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
                       numericInput("r", 
                                    h3("Relative reproductive input (parents)"), 
                                    value = 2/3)
                ),
                column(3, 
                       numericInput("c", 
                                    h3("Sperm competition coefficient"), 
                                    value = 0.5)
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
    if(round(input$Nm)==500 &
       round(input$Nf)==500 &
       round(input$r,2)==round(2/3,2) &
       round(input$c,2)==0.5 &
       input$ws==1 &
       input$wn==1 &
       round(input$wv,3)==round(exp(-0.5/(2*50)),3)){
      
      read_predictions()
    }else{
      create_predictions(Nm=input$Nm,
                         Nf=input$Nf,
                         r=input$r,
                         c=input$c,
                         ws=input$ws,
                         wn=input$wn,
                         wv=input$wv)
    }
  })
  
  # check the inputs and create the subset
  subdat<-reactive({
    data<-get_data()
    validate(no_plots(data, input$sliderCP, input$sliderNP),
             no_rows(data, input$sliderCP, input$sliderNP)
    )
    create_subset(data, input$sliderCP, input$sliderNP)

  })
  
  
  # create the plot
  output$contours <- renderPlotly({
    
   sub_calcs<-subdat()
    
    # fig 1: NS vs CS
    fig1 <- plot_ly(
      x = sub_calcs$NS_freq, 
      y = sub_calcs$CS_freq, 
      z = as.matrix(sub_calcs[,c("CS_rs","NS_rs")]), 
      colorscale=list(seq(0,1,length.out = 9),
                      c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')),
      type = "contour"
    )
    # add axis labels
    x<-list(title="Noncourter-Sneaker frequency (NS RS)")
    y<-list(title="Courter-Sneaker frequency (CS RS)")
    fig1 <- fig1 %>% layout(xaxis=x,yaxis=y)
    # add label to contour names
    fig1 <- fig1 %>% colorbar(title = "Relative RS")
    
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
