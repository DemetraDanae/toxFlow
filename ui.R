#toxFlow GSVA-Read across web tools
#Copyright (C) 2017  Dimitra Danai Varsou

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#____________________________________________________________________________

library("shiny")
library(shinythemes)
library(shinyjs)
library(igraph)
library(network)
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

ui <- (fluidPage(
  
      #For icon on tab!
      list(tags$head(HTML('<link rel="icon", href="http://si.ntua.gr/pyrforos-digamma.png", type="image/png" />'))),
      div(style="padding: 1px 0px; width: '100%'",
          titlePanel(
            title="", windowTitle="toxFlow"
          )
      ),
  
      theme = shinytheme("Flatly"), navbarPage(p(strong("toxFlow"), style = "font-family: 'Agency FB'; Agency FB Bold"),
                                                          
      #Tab 1
      tabPanel("GSVA",
      useShinyjs(),  # Include shinyjs
      sidebarLayout(
      #INPUTS
      div(id="GSVA_in",
      sidebarPanel(h4("Gene set variation analysis"),
      
      icon('anchor'),
      span(textOutput("instr1_1"), style="color:#6495ED"),
      span(textOutput("instr1_2"), style="color:#6495ED"), 
      span(textOutput("instr1_3"), style="color:#6495ED"),br(),             
                   
      tabsetPanel(type = "tabs",
                  #Files
                  tabPanel("Input data",br(),
                       
                        #Files
                        selectInput("Files_gsva","Choose files", c("Import dataset"="Import_gsva",
                                                                      "Use demo dataset"="Files1")),
                        uiOutput("Import_gsva_data"),
                        uiOutput("Import_gsva_class"),
                        #Scaling?
                        checkboxInput("scaling","Scaling of raw data"),
                        selectInput("ids","Accession ID:", c("UNIPROT"="UNIPROT", "REFSEQ"="REFSEQ","ENTREZID"="ENTREZID","SYMBOL"="SYMBOL"))#,
                        ),
                  tabPanel("Parameters of analysis",br(),
                        #Gene Set Collection
                        selectInput("gsc","Gene set collection:", c("C5: GO MF gene sets (from MSigDB)"="c5",
                                                                       "Chemical-GO enriched associations (from Comparative Toxicogenomics Database)"="Tox","Other..."="Other")),
                        uiOutput("Other"),
                        #Bootstraps
                        numericInput("boot","Number of bootstraps:",value=1,min=1,max=NA),
                        #Gene limits
                        h5(strong("Gene set size:")),
                        fluidRow(
                        column(6,uiOutput("slider_geneLimits1")),
                        column(6,uiOutput("slider_geneLimits2"))),
                        #pvalue cut off
                        sliderInput("adjPvalueCutoff","Adjusted p-value cut off:",min=0.001, max=0.05, value=0.001,step=0.001),
                        selectInput("map","GO mapping for acyclic graph:",c("GOMFPARENTS","GOMFCHILDREN")),
                        actionButton("run_gsva","Run analysis"))
                  ))),
                                                                     
     #OUTPUTS
     mainPanel(uiOutput("image1"),
     h4(textOutput("resultsText")),dataTableOutput("resultsTable_gsva"), uiOutput("DownGSVA"), br(),
     h4(textOutput("heatmapText")),plotOutput("heatmap"),uiOutput("DownHeat"),br(),
     #textOutput("prot"),
     h4(textOutput("GOgraphText")),plotOutput("GOGraph"),uiOutput("DownGOgraph"))
     )
     ),
                                                          
     #Tab 2
     tabPanel("Read-across training", 
     sidebarLayout(
     #INPUTS
     div(id="Training_in",
     sidebarPanel(h4("Read-across training of model, using leave-one-out cross-validation method"),
     
     icon('anchor'),
     span(textOutput("instr2_1"), style="color:#6495ED"),
     span(textOutput("instr2_2"), style="color:#6495ED"), 
     span(textOutput("instr2_3"), style="color:#6495ED"),
     span(textOutput("instr2_4"), style="color:#6495ED"),
     span(textOutput("instr2_5"), style="color:#6495ED"),br(), 
                  
     tabsetPanel(type = "tabs", 
          #Files
          tabPanel("Input data",br(),                   
                   selectInput("Files_RA","Choose files", c("Import dataset"="Import_RA",
                                                                                "Use demo dataset"="Files2")),
                          uiOutput("Import_RA_phChem"),
                          uiOutput("Import_RA_bio"),
                          checkboxInput("only_one","Only one file available",value = FALSE),
                   radioButtons("the_file","Available file:", c("Physicochemical file"="only_phchem","Biological file"="only_bio")),                   #Scaling
                   checkboxInput("scaling_descr","Scaling of physicochemical data"),
                   checkboxInput("scaling_bio","Scaling of biological data"),
                   checkboxInput("DEproteins", "Use of differentially expressed genes or proteins from GSVA analysis")#,
                   ), 
          tabPanel("Parameters of analysis",br(), 
                #Attribute filtering
                sliderInput("lvl_phchem","Choose physicochemical attribute filtering level", min=0.05, max=0.7,value=0.05),
                sliderInput("lvl_bio","Choose biological attribute filtering level", min=0.05, max=0.7, value=0.30),
                #Neighboring
                selectInput("calc","Similarity calculation method:",choices = c("Cosine similarity"="cos","Manhattan distance"="manhattan","Euclidean distance"="euclidean"), multiple = FALSE),
                #Thresholds
                uiOutput("slider_PhCh"),
                uiOutput("slider_Bio"),
                #Calculation base
                radioButtons("RAcorrBase","Prediction base on:", c("Physicochemical data"="PhChem","Biological data"="Bio")),
                actionButton("run_tr","Training"), br(),
                #Netplot
                h4("Nanoparticle's universe"),
                uiOutput("find_uni"),
                checkboxInput("size","Nanoparticles' diameter"),
                uiOutput("size"),
                checkboxInput("class","Nanoparticles' classification"),
                uiOutput("class"),br(),
                h5(strong("Physicochemical thresholds:")),
                fluidRow(
                column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot1"))),
                column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot2")))),
                h5(strong("Biological thresholds:")),
                fluidRow(
                column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot1"))),
                column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot2")))),
                actionButton("viz1","Visualization"))
     ))),
                                                                     
     #RESULTS
     mainPanel(uiOutput("image2"),
     h4(textOutput("CorrCoeff")),textOutput("R2"),br(),
     h4(textOutput("TrainingRes")),dataTableOutput("resultsTable"),br(), uiOutput("DownBut1"),
     h4(textOutput("NanoUni")),plotOutput("netAll"))
     )),
          
     #Tab 3
     tabPanel("Read-across prediction", sidebarLayout(
     #INPUTS
     div(id="Prediction_in",
     sidebarPanel(h4("Toxicity endpoint prediction"),
                  
     icon('anchor'),
     span(textOutput("instr3_1"), style="color:#6495ED"),
     span(textOutput("instr3_2"), style="color:#6495ED"), 
     span(textOutput("instr3_3"), style="color:#6495ED"),
     span(textOutput("instr3_4"), style="color:#6495ED"),
     span(textOutput("instr3_5"), style="color:#6495ED"),br(),
      
     #Files
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("descrRaw_p","Physicochemical data:", accept = "text/csv")),
     tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("pcoronaRaw_p","Biological data:", accept = "text/csv")),
     actionButton("run_pred","Prediction"),
     #Netplot
     h4("Nanoparticle's universe"),
     uiOutput("find_uni_pred"),
     checkboxInput("size_p","Nanoparticles' diameter"),
     uiOutput("size_p"),
     checkboxInput("class_p","Nanoparticles' classification"),
     uiOutput("class_p"),
     h5(strong("Physicochemical thresholds:")),
     fluidRow(
       column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot1_p"))),
       column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_PhChplot2_p")))),
     h5(strong("Biological thresholds:")),
     fluidRow(
       column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot1_p"))),
       column(6,tags$div(title="Adjusting color of neighbors in NP universe plot. 1st threshold specifies the closest neighbors and 2nd the middle neighbors.",uiOutput("slider_Bioplot2_p")))),
     actionButton("viz2","Visualization")
     )),
                                                            
     #RESULTS
     mainPanel(uiOutput("image3"),
     h4(textOutput("PredTab")),dataTableOutput("predTable"), uiOutput("DownBut2"),
     h4(textOutput("NanoUniPred")),plotOutput("netAll_p"))
     )),
                                                          
     #Tab4
     tabPanel("Help",p("To open in a new tab click", a("here",target="_blank",href="manual.pdf")),p("Report problems:",a("dimitra.varsou (at) gmail.com",target="_blank",href="mailto:dimitra.varsou@gmail.com")),p("Full code available on:", a("GitHub",target="_blank",href="https://github.com/DemetraDanae/toxFlow")), 
              p("Tired of reading? -Find a brief tutorial on", a("YouTube",target="_blank",href="https://www.youtube.com/watch?v=kGp2PuTiDrg")),tags$iframe(src="manual.pdf", width="900", height="600"))
              
),
tags$footer(span(HTML('<footer>
                      <img src="http://si.ntua.gr/pyrforos-digamma.png", height="80", width="80"</img>
                      </footer>'), align="center"), p(a("Unit of Process Control and Informatics",target="_blank",href="http://www.chemeng.ntua.gr/labs/control_lab/"),br(),"School of Chemical Engineering", br(),"National Technical University of Athens"), align="center")#,br(),"toxFlow Copyright", icon("copyright"), "2017  Dimitra Danai Varsou", align="center")
)
)
