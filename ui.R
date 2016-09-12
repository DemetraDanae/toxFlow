library("shiny")
library(shinythemes)
library(igraph)
library("DT")
library(extrafont)

library(GSEABase)
library(GSVAdata)
library("GSEAlm")
library(org.Hs.eg.db)
library(GSVA)
library(limma)
library(GO.db)
library(GOstats)
library("Rgraphviz")

ui <- (fluidPage(theme = shinytheme("Flatly"), navbarPage(p(strong("toxFlow"), style = "font-family: 'Agency FB'; Agency FB Bold"), #footer =p("National Technical University of Athens", img(src="emp.png", align="center")),
                                                          
      #Tab 1
      tabPanel("GSVA",
      sidebarLayout(
      #INPUTS
      sidebarPanel(
      #Files
      tags$div(title="Requiering .csv file. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("rawData","Biological data:", accept = "text/csv")),
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' phenotype or other categorical variable of interest based on which GSVA values will be calculated.",fileInput("classifData","Data classification:", accept = "text/csv")),
      #Scaling?
      h4("Parameters of analysis"),
      checkboxInput("scaling","Scaling of raw data"),
      #Duplicated values?
      checkboxInput("double","Average duplicate values"),
      selectInput("ids","Accession ID:", c("UNIPROT"="UNIPROT", "REFSEQ"="REFSEQ","ENTREZID"="ENTREZID","SYMBOL"="SYMBOL")),
      #Gene Set Collection
      selectInput("gsc","Gene set collection:", c("C5: GO MF gene sets (from MSigDB)"="c5",
                        "Chemical-GO enriched associations (from Comparative Toxicogenomis Database)"="Tox","Other..."="Other")),
      uiOutput("Other"),
      #Bootstraps
      numericInput("boot","Number of bootstraps:",value=1,min=1,max=NA),
      #Gene limits
      sliderInput("geneLimits","Gene set size:",min=1,max=10000,value=c(1,1000),step=1),
      #pvalue cut off
      sliderInput("adjPvalueCutoff","Adjusted p-value cut off:",min=0.001, max=0.05, value=0.001,step=0.001),
      selectInput("map","GO mapping for acyclic graph:",c("GOMFPARENTS","GOMFCHILDREN")),
      actionButton("run_gsva","Run analysis")            
      ),
                                                                     
     #OUTPUTS
     mainPanel(uiOutput("image1"),
     h4(textOutput("resultsText")),dataTableOutput("resultsTable_gsva"), uiOutput("DownGSVA"), br(),
     h4(textOutput("heatmapText")),plotOutput("heatmap"),br(),
     #textOutput("prot"),
     h4(textOutput("GOgraphText")),plotOutput("GOGraph"))#, uiOutput("DownGraph"))
     )
     ),
                                                          
     #Tab 2
     tabPanel("Read across training", 
     sidebarLayout(
     #INPUTS
     sidebarPanel(h4("Read across training of model, using leave-one-out cross-validation method"),
     #Files
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes). The 1st row must contain the toxicity endpoint for each sample.",fileInput("descrRaw","Physicochemical data:", accept = "text/csv")),
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes). The 1st row must contain the toxicity endpoint for each sample.",fileInput("pcoronaRaw","Biological data:", accept = "text/csv")),
     #Scaling
     h4("Parameters of analysis"),
     checkboxInput("scaling_descr","Scaling of physicochemical data"),
     checkboxInput("scaling_bio","Scaling of biological data"),
     checkboxInput("DEproteins", "Use of differentially expressed genes or proteins from GSVA analysis"),
     actionButton("run_tr","Training"),
     #Neighboring
     selectInput("calc","Affinity calculation method:",choices = c("Cosine similarity"="cos","Manhattan distance"="manhattan","Euclidean distance"="euclidean"), multiple = FALSE),
     #Calculation base
     radioButtons("RAcorrBase","Prediction base on:", c("Physicochemical data"="PhChem","Biological data"="Bio")),
     #Thresholds
     uiOutput("slider_PhCh"),
     uiOutput("slider_Bio"),
     #Netplot
     h4("Nanoparticle's universe"),
     numericInput("ref","Reference nanoparticle:",value=1,min=1,max=NA),
     tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot")),
     tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot")),
     actionButton("viz1","Visualization")
     ),
                                                                     
     #RESULTS
     mainPanel(uiOutput("image2"),
     h4(textOutput("CorrCoeff")),textOutput("R2"),br(),
     h4(textOutput("TrainingRes")),dataTableOutput("resultsTable"), uiOutput("DownBut1"),
     h4(textOutput("NanoUni")),plotOutput("netAll"))
     )),
          
     #Tab 3
     tabPanel("Read across prediction", sidebarLayout(
     #INPUTS
     sidebarPanel(h4("Toxicity endpoint prediction"),
     #Files
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("descrRaw_p","Physicochemical data:", accept = "text/csv")),
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("pcoronaRaw_p","Biological data:", accept = "text/csv")),
     #Scaling
     h4("Parameters of analysis"),
     checkboxInput("scaling_descr_p","Scaling of physicochemical data"),
     checkboxInput("scaling_bio_p","Scaling of biological data"),
     #DE proteins?
     checkboxInput("DEproteins", "Use of differentially expressed genes or proteins from GSVA analysis"),
     actionButton("run_pred","Prediction"),
     #Netplot
     h4("Nanoparticle's universe"),
     numericInput("ref_p","Reference nanoparticle:",value=1,min=1,max=NA),
     tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot_p")),
     tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot_p")),
     actionButton("viz2","Visualization")
     ),
                                                            
     #RESULTS
     mainPanel(uiOutput("image3"),
     h4(textOutput("PredTab")),dataTableOutput("predTable"), uiOutput("DownBut2"),
     h4(textOutput("NanoUniPred")),plotOutput("netAll_p"))
     )),
                                                          
     #Tab4
     #tabPanel("Manual",(a("Manual",target="_blank",href="manual.pdf")))
     tabPanel("Help",p("To open in new tab press", a("here",target="_blank",href="manual.pdf")),p("Address correspondence to: dimitraDOTvarsouATgmailDOTcom"),tags$iframe(src="manual.pdf", width="900", height="600"))
              
),
tags$footer(span(HTML('<footer>
                      <img src="http://si.ntua.gr/pyrforos-digamma.png", height="80", width="80"</img>
                      </footer>'), align="center"), p(a("Unit of Process Control and Informatics",target="_blank",href="http://www.chemeng.ntua.gr/labs/control_lab/"),br(),"School of Chemical Engineering", br(),"National Technical university of Athens"), align="center")
)
)