
# This is a Shiny web application. You can run the application at web:
# https://sherriying.shinyapps.io/BulkRNAseq_Shiny_demo/


library(shiny)
library(bslib)
library(ggplot2)
library(org.Hs.eg.db)
library(GOstats)
library(openxlsx)
library(plyr)
library(BiocManager)
library(dplyr)
options(repos = BiocManager::repositories())

load("./exprData_metaData.RData")

# Define UI for application that draws barplot and boxplot
ui <- 
  navbarPage("Bulk-RNAseq Shiny app demo",
             theme = bs_theme(version = 4,bootswatch = "flatly"),         
  tabPanel("Individual gene expression",
        # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectizeInput('genelist',"Select/enter a gene symbol:", choices = NULL),
          downloadButton("download", "Download")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(type = "tabs",id = "tabselected",
                      tabPanel("Barplot",value=1,
                               plotOutput("barplot")),
                      tabPanel("Boxplot",value=2,
                               plotOutput("boxplot"))
          ) # end with tabsetPanel
        ) # end with mainPanel
    ) # end with sidebarLayout
  ) # end with tabPanel1 in navbar
) # finish navbarPage


# Define server logic required to draw plots
server <- function(input, output,session) {
  # update "selectsize" to the server version,The server-side selectize input uses R to process searching, and R will return the filtered data to selectize. 
  updateSelectizeInput(session, 'genelist', choices = row.names(counts.allsamples.normalized),server = TRUE)
  # call function from other file "utils.R" which including barplot and boxplot
  source("./utils.R")
  # run funciton 
  ind.gene<-reactive({
    individual.gene.barplot(input$genelist)})
  
  output$barplot<-renderPlot({
    validate(
      need(input$genelist != "", "Please select a gene")
    )
    print(ind.gene())
  })
  
 
  # run funciton 
  ind.gene.boxplot<-reactive({individual.gene.boxplot(input$genelist)})
  output$boxplot<-renderPlot({
    print(ind.gene.boxplot())
  })
  
  ####################
  ### download files ##
  ###################
  #  barplot
  ###################
  observe({
    if (input$tabselected==1){
      output$download<- downloadHandler(
        filename = function() { paste(input$genelist, '_barplot.pdf', sep='') },
        content = function(file) {
          #ggsave(file,plot = ind.gene(), device = "pdf")
          pdf(file)
          print(ind.gene())
          dev.off()
        }  
      ) # finish downloadHandler
    } # finish barplot download script
  ###############
  #  boxplot
  ###############
  if (input$tabselected==2){
    output$download<- downloadHandler(  
      filename = function() { paste(input$genelist, '_boxplot.pdf', sep='') },
      content = function(file) {
        ggsave(file,plot = ind.gene.boxplot(), device = "pdf")
      }
    ) # finish downloadHandler
   } # finish boxplot 
 })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
