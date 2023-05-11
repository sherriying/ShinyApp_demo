
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
                          #selectizeInput("genelist", "Select/enter a gene symbol:",
                          #               choices=row.names(counts.allsamples.normalized)
                          #),
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
  ## barplot function
  individual.gene.barplot<-function(input.gene){
    tem<-counts.allsamples.normalized[row.names(counts.allsamples.normalized)==input.gene,]
    tem.df<-data.frame(gene=t(tem),group=meta.data[,"samples"],
                       geno_treat=meta.data$geno.treat.combine)
    tem.df$group<-factor(tem.df$group,levels=unique(tem.df$group))
    tem.df$geno_treat<-factor(tem.df$geno_treat,levels=c("A_ctrl","B_ctrl","A_drug","B_drug"))
    names(tem.df)[1]<-"gene"
    plots<-ggplot(tem.df, aes(x=group, y=gene,fill=geno_treat) )+ 
      geom_bar(stat="identity",color="black",width=.5)+
      theme_minimal()+
      ggtitle(input.gene)+
      ylab("Normalized expression")+
      #labs(title="Apln", y = "Counts")+
      theme(
        plot.title = element_text(face="bold",hjust=0.5),
        legend.position = "none",
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=7,angle = 50,hjust=1),
        text = element_text(size=9,face="bold"))
  }
  # run funciton 
  ind.gene<-reactive({
    individual.gene.barplot(input$genelist)})
  
  output$barplot<-renderPlot({
    validate(
      need(input$genelist != "", "Please select a gene")
    )
    print(ind.gene())
  })
  
  ## boxplot function
  individual.gene.boxplot<-function(input.gene){
    tem<-counts.allsamples.normalized[row.names(counts.allsamples.normalized)==input.gene,]
    tem.df<-data.frame(gene=t(tem),group=meta.data[,"samples"],
                       geno_treat=meta.data$geno.treat.combine)
    tem.df$group<-factor(tem.df$group,levels=unique(tem.df$group))
    tem.df$geno_treat<-factor(tem.df$geno_treat,levels=c("A_ctrl","B_ctrl","A_drug","B_drug"))
    names(tem.df)[1]<-"gene"
    ggplot(tem.df, aes(geno_treat, gene,fill=geno_treat) )+ 
      geom_boxplot(varwidth=T)+
      #geom_smooth(aes(group = 1), colour = "blue")+
      theme_minimal()+
      ggtitle(input.gene)+
      ylab("Normalized expression")+
      #labs(title="Apln", y = "Counts")+
      theme(
        plot.title = element_text(face="bold",hjust=0.5),
        legend.position = "none",
        #panel.grid.major = element_blank(), 
        #panel.grid.minor = element_blank(),
        #panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_text(size=9,angle = 50,hjust=1),
        text = element_text(size=9,face="bold")
      )
  }
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