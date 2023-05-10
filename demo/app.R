# This is a Shiny web application. You can run the application at web:
# thttps://sherriying.shinyapps.io/BulkRNAseq_Shiny_demo/


library(shiny)
library(bslib)
library(ggplot2)
library(org.Hs.eg.db)
library(GOstats)
library(openxlsx)
library(plyr)
library(dplyr)
library(BiocManager)
options(repos = BiocManager::repositories())

load("./exprData_metaData.RData")
load("./GO_input.RData")

# Define UI for application that draws a histogram
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
             ), # end with tabPanel1 in navbar
             #tabPanel("Component 2"),
             tabPanel("GO enrichment analysis",
                      sidebarLayout(
                        sidebarPanel(
                          selectizeInput("datasets", "Select a dataset:",
                                         choices=names(input.genelist)[1:4]
                          ), # end with selectizeInput
                          br(),
                          radioButtons("GOtype", "GO category",
                                       c("Biological Process"="BP","Molecular Function"="MF","Cellular Component"="CC")
                          ), # end with radioButtons
                          br(),
                          sliderInput("category_small","size of smallest GO term",
                                      value=10,
                                      min=5,
                                      max=20
                          ),# endwith sliderInput for smallest categore
                          sliderInput("category_largest","size of largest GO term",
                                      value=300,
                                      min=50,
                                      max=1000
                          ),# endwith sliderInput for largest categore
                          br(),
                          sliderInput("p","p value",
                                      value=0.001,
                                      min=0,
                                      max=0.05
                          ),# endwith sliderInput
                          downloadButton("downloadGO", "Download")
                        ), # end with sidebarpanel
                        mainPanel(
                          tabsetPanel(type = "tabs",id = "tabselected2",
                                      tabPanel("Table",value=3,
                                               tableOutput("GOtable")),
                                      tabPanel("Plot",value=4,
                                               plotOutput("GOplot"))
                          ) # end with tabsetPanel
                        )# end mainPanel in second tab
                      ) # end with sidebarLayout
             ) # end with tabPanel2 in navbar
  ) # finish navbarPage


# Define server logic required to draw a histogram
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
  
  #######################
  #####GO analysis ######
  ######################
  
  ## GO enrichment funciton
  GOanalysis=function(x,y,z,t1,q1,pCutoff,smallSize,bigSize){
    pCutoff=pCutoff
    params=new("GOHyperGParams",geneIds=x,universeGeneIds=y,annotation=z,ontology=t1,pvalueCutoff=pCutoff,conditional=TRUE,testDirection="over")
    print("Please waiting a bit...")
    GO_analysis=hyperGTest(params)
    GO_analysis_summary=summary(GO_analysis)
    GO_analysis_summary=GO_analysis_summary[GO_analysis_summary[,"Size"]>smallSize & GO_analysis_summary[,"Size"]<bigSize,]
    #write.xlsx(GO_analysis_summary, file=paste(otherName,t1,"_",q1,".xlsx",sep=""))
    ## following for getting gene symbols
    selectedGO_symbol_name_GOterm=NULL
    x_selectGeneId=geneIdsByCategory(GO_analysis)
    x_selectGeneId_selectedGo=x_selectGeneId[names(x_selectGeneId) %in% GO_analysis_summary[,1]]
    GOID_slectedGO=data.frame(GOID=names(x_selectGeneId_selectedGo))
    selectedGO_GOID_term=merge(GOID_slectedGO,GO_analysis_summary[,c(1,7)],by.x="GOID",by.y=names(GO_analysis_summary)[1])
    thesymbs<-toTable(org.Hs.egSYMBOL)
    thegeneName<-toTable(org.Hs.egGENENAME)
    temp<-lapply(x_selectGeneId_selectedGo,function(x){
      selectedGo_symbol=thesymbs[thesymbs[,1] %in% x,]
      selectedGo_genename=thegeneName[thegeneName[,1] %in% x,]
      selectedGo_symbol_name=merge(selectedGo_symbol,selectedGo_genename,by.x="gene_id",by.y="gene_id")
    })
    temp<-ldply(temp,data.frame)
    selectedGO_symbol_name_GOterm=merge(temp,selectedGO_GOID_term,by.x=".id",by.y="GOID")
    names(selectedGO_symbol_name_GOterm)[1]<-"GOID"
    #write.xlsx(selectedGO_symbol_name_GOterm, file=paste(otherName,w,"GO_genelist",".xlsx",sep=""))
    ## collapse gene symbol and output with GO term
    GOID.symbol<-lapply(GO_analysis_summary[,1],function(x){
      temp.symbol=paste(selectedGO_symbol_name_GOterm[selectedGO_symbol_name_GOterm[,"GOID"]==x,"symbol"],collapse=" ")
    })
    names(GOID.symbol)<-GO_analysis_summary[,1]
    GOID.symbol<-ldply(GOID.symbol,data.frame)
    names(GOID.symbol)<-c("GOID","symbol")
    selectedGO_symbol_name_GOterm_symbolCollapse=merge(GO_analysis_summary,GOID.symbol,by.x="GOBPID",by.y="GOID")
    selectedGO_symbol_name_GOterm_symbolCollapse=selectedGO_symbol_name_GOterm_symbolCollapse[order(selectedGO_symbol_name_GOterm_symbolCollapse$Pvalue),]
    #write.xlsx(selectedGO_symbol_name_GOterm_symbolCollapse, file=paste(otherName,w,"GOterm_genesymbol",".xlsx",sep=""))
    return(selectedGO_symbol_name_GOterm_symbolCollapse)
  }
  # run funciton 
  GO.result<-reactive({GOanalysis(x=input.genelist[[input$datasets]],y=genes.universal,z="org.Hs.eg.db",t1=input$GOtype,q1="over",pCutoff=input$p,
                                  smallSize =input$category_small,bigSize = input$category_largest)
  }
  )
  # output bable
  output$GOtable <- renderTable({
    GO.result()
  })
  
  # GO barplot
  GObarplot<-function(genelist){
    if(dim(GO.result())[1]>=50){
      genelist<-GO.result()[1:50,]
    }
    else(genelist<-genelist)
    genelist<-genelist %>% mutate (-log10(Pvalue))
    genelist$Term<-factor(genelist$Term,levels=genelist$Term)
    p<-ggplot(genelist,aes(x=-log10(Pvalue),y=Term))+
      geom_bar(stat="identity",width=0.5,fill="#78afc8")+
      ylab("")+
      theme(
        plot.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line=element_line(color="black"),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(size=8,hjust=1),
        text = element_text(size=9)
      )
  }
  # run funciton 
  GO.barplot<-reactive({GObarplot(GO.result())})
  output$GOplot<-renderPlot({
    print(GO.barplot())
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
    
    ###############
    #  GO table
    ###############
    if (input$tabselected2==3){
      output$downloadGO<- downloadHandler(
        filename = function() {
          paste(input$datasets,"_",input$GOtype, '.xlsx', sep = "")
        },
        content = function(file) {
          write.xlsx(GO.result(), file)
        }
      )
    }# finish GO table download
    
    ###############
    #  GOplot
    ###############
    if (input$tabselected2==4){
      output$downloadGO<- downloadHandler(  
        filename = function() { paste(input$datasets, '_GO_barplot.pdf', sep='') },
        content = function(file) {
          ggsave(file,plot = GO.barplot(), width=7,height=5,device = "pdf")
        }
      ) # finish downloadHandler
    } # finish boxplot
  }) # finish observe
  
}

# Run the application 
shinyApp(ui = ui, server = server)