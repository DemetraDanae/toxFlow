#toxFlow GSVA-Read across web tools

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
options(shiny.maxRequestSize=250*1024^2)

#Functions_____________________________________________________________________________
cos.sim <- function(X,n) #X matrix, n number of columns
{
  C<- matrix(nrow=n,ncol=n)
  for (i in 1:n){
    for(j in 1:n) {
      A <-X[,i]
      B <-X[,j]
      C[i,j]<-sum(A*B)/sqrt(sum(A^2)*sum(B^2))
    }
  }
  return("Similarity Matrix"=C)
}

r2.funct.lm<- function(y,y.new){#y==y, y.new=predicted
  #only for lm with intercept
  R2<- cor(y,y.new)^2
  return(R2)
}

scl<-function(X,n,total){
  X_sc<-data.frame()
  for (i in 1:n){
    min<-min(X[i,])
    max<-max(X[i,])
    for (j in 1:total){
      X_sc[i,j]<-(X[i,j]-min)/(max-min)
    }
  }
  row.names(X_sc)<-rownames(X)
  colnames(X_sc)<-colnames(X)
  
  return(X_sc)
}

filt<-function(X_sc,method,lvl){ #Filtering data X_sc
  
  X_t<-t(X_sc) #transpose data matrix
  
  #Calculating correlations between n.cell and factors
  corrNcellDescr<-cor(X_t, y=NULL, use="everything", method = method)
  
  range<-nrow(corrNcellDescr)
  j<-c()
  for (i in 2:range){
    if ((is.na(corrNcellDescr[i,1])==TRUE)||(abs(corrNcellDescr[i,1])<lvl)){
      j[i]<-i
    }
  }
  j<-na.omit(j)
  j<-as.vector(j)
  X_sc<-X_sc[-j,]
  
  X_sc<-X_sc[-1,] #drop 1st row (n.cell)
  X_sc<-na.omit(X_sc)
  return(X_sc)
}

relfun<-function(calc,X_sc){
  
  #Calculating (cosine) similarity between nanoparticles
  if (calc=="cos"){
    n<-ncol(X_sc)
    relX<-cos.sim(X_sc,n)
    rownames(relX)<-colnames(X_sc)
    colnames(relX)<-colnames(X_sc)
    #relX_sc<-relX
  }else{
    #Calculating distance between nanoparticles
    X1<-t(X_sc)
    relX<-dist(X1, method = calc, diag = FALSE)
    relX<-as.matrix(relX)
    relX<-1/(1+relX)
    }
  
  return(relX)
}

nanoID_crossval<-function(X){
  nano_name<-colnames(X)
  nano_number<-c(1:length(nano_name))
  nano_n.cell<-t(X[1,])
  nanoID<-data.frame(nano_number,nano_name,nano_n.cell,stringsAsFactors = TRUE)
  colnames(nanoID)<-c("#", "name", "net.cell")
  return(nanoID)
}

nanoID_pred<-function(X){
  nano_name<-colnames(X)
  nano_number<-c(1:length(nano_name))
  nanoID<-data.frame(nano_number,nano_name,stringsAsFactors = TRUE)
  colnames(nanoID)<-c("#", "name")
  return(nanoID)
}

RAVfun<-function(RAcorrBase,selData,nano_ID){
  #Calculating Read Across Value
  if (RAcorrBase=="PhChem"){
    col<-"PhChem_corr"
  }
  if (RAcorrBase=="Bio") {
    col<-"Bio_corr"
  }
  
  RAV<-0
  range<-nrow(selData)
  denom<-sum(selData[2:range,col]) #denominator
  for (i in 2:range){
    #print(paste0("denominator:",denom,"i",i,"selData",selData[i,col],"nanoId:",nano_ID[selData$"nano_name"[i],3], "selDataNanoName", selData$"nano_name"[i],"\n"))
    
    RAV<-RAV+(selData[i,col]*nano_ID[selData[i,1],3])/denom
  }
  #print(selData)
  return(RAV)
}

selData<-function(PhCh_th,Bio_th,calc,ref,relDescr,relBio,total,nano_ID){
  
  selData<-data.frame(stringsAsFactors = TRUE)
  selData[1,1]<-nano_ID[ref,1] #nano_number
  selData[1,2]<-nano_ID[ref,2] #nano_name
  selData[1,3]<-relDescr[ref,ref] 
  selData[1,4]<-relBio[ref,ref]
  
  #Selection of neighboring data
  i<-2 #selData rows
  NP<-1 #Nanoparticle
  while (NP<=total){
    if (NP!=ref){
      if ((relDescr[NP,ref]>=PhCh_th)&(relBio[NP,ref]>=Bio_th)){  #Thresholds
        selData[i,1]<-nano_ID[NP,1] #nano_number
        selData[i,2]<-nano_ID[NP,2] 
        selData[i,3]<-relDescr[NP,ref] #corr_factor
        selData[i,4]<-relBio[NP,ref] #corr_factor (pcorona)
        i<-i+1
      }
    }
    NP<-NP+1
  }

  colnames(selData)<-c("#","nano_name", "PhChem_corr","Bio_corr")
  
  return(selData)}

#Gene Set Collections
#setwd("./www")

#1st Gene Set Collection
mSigDB1 <- readLines("c5.mf.v5.1.entrez.gmt") #entrez id's
mSigDB1 <- strsplit(mSigDB1, "\t")

range<-length(mSigDB1)

nam<-vector()
for (i in 1:range){
  nam[i]<-paste(unlist(strsplit(tolower(mSigDB1[[i]][1]),'_')),collapse=' ')
}

link<-data.frame()
for (i in 1:range){
  link[i,1]<-nam[i]
  link[i,2]<-mSigDB1[[i]][2] #url
}
colnames(link)<-c("Terms","links")

names(mSigDB1) <- sapply(nam, function(x) x[1])
mSigDB1 <- sapply(mSigDB1, function(x) x[3:length(x)])

#2nd Gene Set Collection
CTD<-load("CTD_terms.RData")
CTD_gsc<-CTD_gsc

go<-read.table("GOids-terms.txt",header=TRUE)

server <- function(input, output) {
  
#Tab 1: GSVA
  
#1. UI outputs__________________________________________________________________________
  
  #Image
  output$image1<-renderUI(if(input$run_gsva==FALSE){
    img(src="https://raw.githubusercontent.com/DemetraDanae/toxFlow/master/toxFlow_logo.png", align="center")
  })
  
  #Files
  output$Import_gsva_data<-renderUI({
    if(input$Files_gsva=="Import_gsva"){ 
      tags$div(title="Requiering .csv file. Columns must contain the NPs (or samples) and rows the markers (or genes).",fileInput("rawData", "Biological data:", accept = "text/csv"))}
  })
  
  output$Import_gsva_class<-renderUI({
    if(input$Files_gsva=="Import_gsva"){ 
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' phenotype or other categorical variable of interest based on which GSVA values will be calculated.",fileInput("classifData", "Data classification:", accept = "text/csv"))}
  })
  
  output$instr1_1<-renderText("Please insert data or use the demo dataset")
  
  output$instr1_2<-renderText("Ready to run analysis")
  
  output$instr1_3<-renderText("Please load a gene set collection")
  
  #Other gsc
  output$Other<-renderUI({
    if(input$gsc=="Other"){
      tags$div(title="The file must be a .csv. The first column must contain GOterms and the second the ENTREZids",fileInput("other", " ", accept = "text/csv"))}
  })
  
  output$slider_geneLimits1<-renderUI(numericInput("geneLimits1"," ",min=1,max=5000,value=1,step=1))
    
  output$slider_geneLimits2<-renderUI(numericInput("geneLimits2"," ",min=input$geneLimits1,max=10000,value=6000,step=1))
  
  output$DownGSVA<-renderUI(if(input$run_gsva){downloadButton("Down_gsva", label = "Download")})
  
  output$DownHeat<-renderUI(if(input$run_gsva){downloadButton("Down_Heat",label = "Download")})
  
  output$DownGOgraph<-renderUI(if(input$run_gsva){downloadButton("Down_Gograph",label = "Download")})
  
#2. Start GSVA______________________________________________________________
  
  #import files
  rawData<-reactive({
    if (input$Files_gsva=="Files1"){
      raw<-read.csv2("bio_gsva.csv",header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    }else{
      raw<-read.csv2(input$rawData$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
      raw<-t(raw)} #new format of files required
    raw})
  
  additionalInfo<-reactive({
    if (input$Files_gsva=="Files1"){
      add<-read.csv2("NP_classification.csv",header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    }else{
      add<-read.csv2(input$classifData$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)}
    add
  })
  
  #Enable/disable run button when files are not uploaded
  observe({
    if ((((is.null(input$rawData) || input$rawData == "") || is.null(input$classifData) || input$classifData == "") & input$Files_gsva!="Files1") || (input$gsc=="Other" & is.null(input$other))){
      shinyjs::disable("run_gsva")
      } else {
      shinyjs::enable("run_gsva")
      }
  })
  
  #Messages if files not uploaded
  observe({
    if ((((is.null(input$rawData) || input$rawData == "") || is.null(input$classifData) || input$classifData == "") & input$Files_gsva!="Files1")){
      shinyjs::show("instr1_1")
      shinyjs::hide("instr1_2")
    } else {
      shinyjs::hide("instr1_1")
      shinyjs::show("instr1_2")
    }
  })
  
  #Messages if gsc not uploaded
  observe({
    if (input$gsc=="Other" & is.null(input$other)){
      shinyjs::show("instr1_3")
      shinyjs::hide("instr1_2")
    } else {
      shinyjs::hide("instr1_3")
      shinyjs::show("instr1_2")
    }
  })
  
  #Disable ID for demo dataset
  observe({
    if (input$Files_gsva=="Files1"){
      shinyjs::disable("ids")
    } else {
      shinyjs::enable("ids")
    }
  })
  
  #double values?
  uniqueData<-reactive({
        X1<-rawData()
  X1})
  
  #Scaling?
  rawData_sc<- reactive({
    if (input$scaling==TRUE){
      rawData_sc<-scl(uniqueData(),nrow(uniqueData()),ncol(uniqueData()))
    } else{
      rawData_sc<-data.frame()
      rawData_sc<-uniqueData()
    }
    rawData_sc
  })
  
  #Find Entrez ids
  X<-reactive({
    protIDs<-row.names(rawData_sc())
    protIDs<-as.character(protIDs)
    
    validate(
      need(sum(protIDs%in%keys(org.Hs.eg.db, input$ids))>0, "None of the keys entered are valid keys for the provided Accession ID. Please use the right Accession ID to see a listing of valid arguments.")
    )
    
    if (input$ids!="ENTREZID"){
      
      annots<-select(org.Hs.eg.db,keys=protIDs,keytype = input$ids, columns=c("ENTREZID"))
      
      X<- merge(rawData_sc(), annots, by.x=0, by.y=input$ids)
      X <- na.omit(X) #exclude NA values
      rownames(X) <- seq(length=nrow(X)) #reset rows' numbering
      X<-X[!duplicated(X[,1]), ] #exclude double values
      rownames(X) <- seq(length=nrow(X))
      rownames(X)<-X$ENTREZID
      X<- subset(X, select = -c(Row.names,ENTREZID))
    }else{
      X<-rawData_sc()
    }
    X})
  
  #geneSet?
  geneSetcol<-reactive({ #same as geneSet
    if(input$gsc=="c5"){
      geneSetcol<-mSigDB1
    } else if (input$gsc=="Tox"){
      geneSetcol<-CTD_gsc
    } else if (input$gsc=="Other"){
      entr<-read.csv2(input$other$datapath,header=TRUE,sep=";",quote = "\"",dec=".", stringsAsFactors = FALSE,strip.white = TRUE)
      
      entr<-entr[complete.cases(entr),]
      singleIDs<-unique(entr[,1])
      #---------------------------------------------------------------------------------
      
      gsc<-list()
      rng_singleIDs<-length(singleIDs)
      gsc[[1]]<-entr[,2][which(entr[,1]%in%singleIDs[1])]
      for (i in 2:rng_singleIDs){
        gsc[[i]]<-entr[,2][which(entr[,1]%in%singleIDs[i])]
        gsc[[i]]<-unique(gsc[[i]])
      }
      
      terms<-vector()
      for (i in 1:rng_singleIDs){
        terms[i]<- paste(unlist(strsplit(singleIDs[i],'-')),collapse=' ')}
      
      names(gsc)<-terms  
      
      geneSetcol<-gsc
    }
    geneSetcol
  })
  
  #ES calculation
  ES<- reactive({
    
    #geneSet?
    if(input$gsc=="c5"){ #same as geneSetcol
      geneSet<-mSigDB1
    } else if (input$gsc=="Tox"){
      geneSet<-CTD_gsc
    } else if (input$gsc=="Other"){
      entr<-read.csv2(input$other$datapath,header=TRUE,sep=";",quote = "\"",dec=".", stringsAsFactors = FALSE,strip.white = TRUE)
      
      entr<-entr[complete.cases(entr),]
      singleIDs<-unique(entr[,1])
      #---------------------------------------------------------------------------------
      
      gsc<-list()
      rng_singleIDs<-length(singleIDs)
      gsc[[1]]<-entr[,2][which(entr[,1]%in%singleIDs[1])]
      for (i in 2:rng_singleIDs){
        gsc[[i]]<-entr[,2][which(entr[,1]%in%singleIDs[i])]
        gsc[[i]]<-unique(gsc[[i]])
      }
      
      terms<-vector()
      for (i in 1:rng_singleIDs){
        terms[i]<- paste(unlist(strsplit(singleIDs[i],'-')),collapse=' ')}
      
      names(gsc)<-terms  
      
      geneSet<-gsc
    }
    
    X1<-sdsgenes <- apply(X(), 1, sd) #same as X_filt
    X1<- X()[-which(sdsgenes<=quantile(sdsgenes,0.25)),]
    
    #set.seed(11293)
    min<-input$geneLimits1
    max<-input$geneLimits2
    #old version of GSVA function
    #gsva_score<-gsva(as.matrix(X1),geneSet,method="gsva",rnaseq=FALSE,abs.ranking=FALSE, min.sz=min,max.sz=max,no.bootstraps=input$boot,bootstrap.percent =.632,mx.diff=TRUE,verbose=TRUE)
    #ES<-gsva_score$es.obs
    #new version of GSVA (less arguments)
    gsva_score<-gsva(as.matrix(X1),geneSet,method="gsva",abs.ranking=FALSE, min.sz=min,max.sz=max,mx.diff=TRUE,verbose=TRUE)
    ES<-gsva_score
    print("ok2")
    ES
  })
  
  X1terms<-reactive({
    
    X1<-sdsgenes <- apply(X(), 1, sd) 
    X1<- X()[-which(sdsgenes<=quantile(sdsgenes,0.25)),]
    
    X1terms<-data.frame()
    range_geneSets<-length(geneSetcol())
    range_genes<-nrow(X1)
    
    k<-1
    i<-1
    while (i<=range_geneSets){
      l<-length(geneSetcol()[[i]])
      if (l==0){
        i<-i+1
      } else {
        for (j in 1:l){
          if (geneSetcol()[[i]][j]%in%rownames(X1)){
            X1terms[k,1]<-geneSetcol()[[i]][j]
            X1terms[k,2]<-names(geneSetcol())[[i]]
            k<-k+1
          }
        }
        i<-i+1
      }}
    
    for (i in 1:range_genes){
      if (rownames(X())[i]%in%X1terms[,1]==FALSE){
        X1terms[k,1]<-rownames(X1)[i]
        X1terms[k,2]<-NA
        k<-k+1
      }
    }
    
    names(X1terms)<-c("genes","terms")
    X1terms
  })  
  
  #Classification
  classif <-reactive({
    model.matrix(~factor(additionalInfo()$Classification))}) 
  
  #pathway-level differential expression analysis
  
  DEgeneSets<-reactive({
    adjPvalueCutoff <-input$adjPvalueCutoff
    design <- classif()
    fit <- lmFit(ES(), design)
    fit <- eBayes(fit)
    allGeneSets <- topTable(fit, coef=colnames(design)[2], number=Inf)
    DEgeneSets <- topTable(fit, coef=colnames(design)[2], number=Inf,p.value=adjPvalueCutoff, adjust="BH")
    res <- decideTests(fit, p.value=adjPvalueCutoff)
    
    DEgeneSets
  })
  
  #DE gene sets included_______________________________________________________________
  
  #DE proteins________________________________________________________________
  DEproteins<-reactive({
    DEgeneSets<-DEgeneSets()
    range<-nrow(DEgeneSets)
    DEgeneSets_terms<-vector()
    for (i in 1:range){
      DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(rownames(DEgeneSets)[i]),'_')),collapse=' ')
    }
    
    X1.3<-list() #Entrezids of DEgeneSets
    X1.3<-X1terms()[,1][which(X1terms()[,2]%in%DEgeneSets_terms)]
    #X1.3<-X1.2[,1][which(X1.2[,2]%in%rownames(DEgeneSets))]
    X1.3<-X1.3[!duplicated(X1.3)] #exclude double values
    
    protIDs<-row.names(rawData_sc())
    protIDs<-as.character(protIDs)
    annots<-select(org.Hs.eg.db,keys=protIDs,keytype = input$ids, columns=c("ENTREZID"))
    DEproteins<-annots$UNIPROT[which(annots$ENTREZID%in%X1.3)]
    
    DEproteins
    #print(DEproteins)
  })
  
  #output$prot<-renderText({DEproteins()})
  
  #Constructing final table_________________________________________________________
  
  results<-eventReactive(input$run_gsva,{
    
    if (input$gsc=="c5"){
      withProgress(message = "Processing...", {
        DEgeneSets<-DEgeneSets()
        
        #final table
        range<-nrow(DEgeneSets)
        DEgeneSets_terms<-vector()
        for (i in 1:range){
          DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(rownames(DEgeneSets)[i]),'_')),collapse=' ')
        }
        
        results<-data.frame(stringsAsFactors = TRUE)
        for (i in 1:range){
          if (DEgeneSets_terms[[i]]%in%go[,2]){
            A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
            A<-as.character(A)
            results[i,1]<-A} #GOMFID
          else{
            results[i,1]<-"Not included"
          }
          results[i,2]<-DEgeneSets_terms[[i]] #term
          results[i,3]<-link[,2][which(link[,1]==DEgeneSets_terms[[i]])] #link
        }
        
        results<-cbind(results,DEgeneSets$adj.P.Val) #pvalue
        
        size<-vector()
        for (i in 1:range){
          y<-geneSetcol()[which(names(geneSetcol())==DEgeneSets_terms[[i]])]
          size[i]<-length(y[[1]])
        }
        
        #Counts
        counts<-data.frame()
        for (i in 1:range){
          counts[i,1]<-DEgeneSets_terms[i]
          if (DEgeneSets_terms[i]%in%X1terms()[,2]==TRUE){
            counts[i,2]<-length(X1terms()[,2][which(X1terms()[,2]==counts[i,1])])
          } else {
            counts[i,2]<-0
          }
        }
        colnames(counts)<-c("Terms", "Counts")
        
        results<-cbind(results,size)
        results<-cbind(results,counts[,2])
        
        colnames(results)<-c("GOMFID","Term","url","p-value","Size","Count")
        
        results_total<-results
        results_total$url<-paste0("<a href='",results_total$url,"'>",results_total$Term,"</a>")
        results_total<-results_total[,-2]
        colnames(results_total)<-c("GOMFID","Term","p-value","Size","Counts")
      })
    } else {
      withProgress(message = "Processing...", {
        DEgeneSets<-DEgeneSets()
        
        range<-nrow(DEgeneSets)
        DEgeneSets_terms<-vector()
        for (i in 1:range){
          A<-as.character(rownames(DEgeneSets)[i])
          DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(A),',')),collapse=' ')
        }
        
        results<-data.frame(stringsAsFactors = TRUE)
        for (i in 1:range){
          if (DEgeneSets_terms[[i]]%in%go[,2]){
            A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
            A<-as.character(A)
            results[i,1]<-A} #GOMFID
          else{
            results[i,1]<-"Not included"
          }
          results[i,2]<-DEgeneSets_terms[[i]] #term
        }
        
        results<-cbind(results,DEgeneSets$adj.P.Val) #pvalue
        
        size<-vector()
        for (i in 1:range){
          y<-geneSetcol()[which(DEgeneSets_terms[[i]]==names(geneSetcol()))]
          size[i]<-length(y[[1]])
        }
        
        #Counts
        counts<-data.frame()
        for (i in 1:range){
          counts[i,1]<-DEgeneSets_terms[i]
          if (DEgeneSets_terms[i]%in%X1terms()[,2]==TRUE){
            counts[i,2]<-length(X1terms()[,2][which(X1terms()[,2]==counts[i,1])])
          } else {
            counts[i,2]<-0
          }
        }
        colnames(counts)<-c("Terms", "Counts")
        
        results<-cbind(results,size)
        results<-cbind(results,counts[,2])
        
        colnames(results)<-c("GOMFID","Term","p-value","Size","Counts")
        
        results_total<-results
        
      })
    }
    final<-datatable(results_total, escape = FALSE)
    final 
    #heat.map()
  })
  
  results.Text<-eventReactive(input$run_gsva,{"Differentially expressed gene sets"})
  
  #Heatmap
  heatmap.Text<-eventReactive(input$run_gsva,{
    "Gene sets-genes heatmap"})
  
  heat.map<-eventReactive(input$run_gsva,{
    
    ES_heat<-ES()
    
    total<-nrow(additionalInfo())
    col<-c()
    for (i in 1:total){
      if (additionalInfo()[i,1]==levels(additionalInfo()[,1])[1]){
        col[i]<-"#66CD00"
      } else{
        col[i]<-"#8B2500"
      }
    }
    
    heatmap(ES_heat,margins=c(5,3.5), scale="none",cexRow = 0.7, cexCol = 0.5, ColSideColors=col)
    legend("left",      # location of the legend on the heatmap plot
           legend = c(levels(additionalInfo()[,1])[1], levels(additionalInfo()[,1])[2]), # category labels
           col = c("#66CD00", "#8B2500"),  # color key
           lty= 1,             # line style
           lwd = 5,
           cex=0.5 # line width
    )
  })
  
  GOgraph.Text<-eventReactive(input$run_gsva,{"GO directed acyclic graph"})
  
  GO.graph<-eventReactive(input$run_gsva,{
    
    DEgeneSets<-DEgeneSets()
    
    range<-nrow(DEgeneSets)
    DEgeneSets_terms<-vector()
    for (i in 1:range){
      DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(rownames(DEgeneSets)[i]),'_')),collapse=' ')
    }
    
    results<-vector()
    for (i in 1:range){
      if (DEgeneSets_terms[[i]]%in%go[,2]){
        A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
        A<-as.character(A)
        results[i]<-A} #GOMFID
      else{
        results[i]<-"Not included"
      }
    }
    
    g<-results[which(results!="Not included")]
    
    if (input$map=="GOMFPARENTS"){
      map<-GOMFPARENTS
    } else{
      map<-GOMFCHILDREN
    }
    g_plot<-GOGraph(g,map)
    plot(g_plot)
  })
  
#3. Outputs___________________________________________________________________________
  output$resultsText<-renderText({
    results.Text()
  })
  
  output$resultsTable_gsva<-renderDataTable({
    results()
  })
  
  output$Down_gsva<- downloadHandler(filename = function() {
    paste("GSVA_analysis", Sys.Date(), '.csv', sep='')},
    
    if(input$gsc=="c5"){
      content = function(file1) {
        DEgeneSets<-DEgeneSets()
        range<-nrow(DEgeneSets)
        DEgeneSets_terms<-vector()
        for (i in 1:range){
          DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(rownames(DEgeneSets)[i]),'_')),collapse=' ')
        }
        
        results<-data.frame(stringsAsFactors = TRUE)
        for (i in 1:range){
          if (DEgeneSets_terms[[i]]%in%go[,2]){
            A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
            A<-as.character(A)
            results[i,1]<-A} #GOMFID
          else{
            results[i,1]<-"Not included"
          }
          results[i,2]<-DEgeneSets_terms[[i]] #term
        }
        
        results<-cbind(results,DEgeneSets$adj.P.Val) #pvalue
        
        size<-vector()
        for (i in 1:range){
          y<-geneSetcol()[which(DEgeneSets_terms[[i]]==names(geneSetcol()))]
          size[i]<-length(y[[1]])
        }
        
        #Counts
        counts<-data.frame()
        for (i in 1:range){
          counts[i,1]<-DEgeneSets_terms[i]
          if (DEgeneSets_terms[i]%in%X1terms()[,2]==TRUE){
            counts[i,2]<-length(X1terms()[,2][which(X1terms()[,2]==counts[i,1])])
          } else {
            counts[i,2]<-0
          }
        }
        colnames(counts)<-c("Terms", "Counts")
        
        results<-cbind(results,size)
        results<-cbind(results,counts[,2])
        
        colnames(results)<-c("GOMFID","Term","p-value","Size","Counts") 
        write.csv2(results,file1)}
    } else{
      content = function(file1) {
        DEgeneSets<-DEgeneSets()
        
        range<-nrow(DEgeneSets)
        DEgeneSets_terms<-vector()
        for (i in 1:range){
          A<-as.character(rownames(DEgeneSets)[i])
          DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(A),',')),collapse=' ')
        }
        
        results<-data.frame(stringsAsFactors = TRUE)
        for (i in 1:range){
          if (DEgeneSets_terms[[i]]%in%go[,2]){
            A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
            A<-as.character(A)
            results[i,1]<-A} #GOMFID
          else{
            results[i,1]<-"Not included"
          }
          results[i,2]<-DEgeneSets_terms[[i]] #term
        }
        
        results<-cbind(results,DEgeneSets$adj.P.Val) #pvalue
        
        size<-vector()
        for (i in 1:range){
          y<-geneSetcol()[which(names(geneSetcol())==DEgeneSets_terms[[i]])]
          size[i]<-length(y[1])
        }
        
        #Counts
        counts<-data.frame()
        for (i in 1:range){
          counts[i,1]<-DEgeneSets_terms[i]
          if (DEgeneSets_terms[i]%in%X1terms()[,2]==TRUE){
            counts[i,2]<-length(X1terms()[,2][which(X1terms()[,2]==counts[i,1])])
          } else {
            counts[i,2]<-0
          }
        }
        colnames(counts)<-c("Terms", "Counts")
        
        results<-cbind(results,size)
        results<-cbind(results,counts[,2])
        
        colnames(results)<-c("GOMFID","Term","p-value","Size","Counts")
        write.csv2(results,file1)
      }
    }
  )
  
  output$heatmapText<-renderText({heatmap.Text()})
  
  output$heatmap<-renderPlot({
    heat.map()
  })
  
  output$Down_Heat<- downloadHandler(filename = function() {
    paste("Heatmap", Sys.Date(), '.png', sep='')
  },
  content = function(file5) {
    png(file5,width=15,height=6,units="in",res=100)
    ES_heat<-ES()
    
    total<-nrow(additionalInfo())
    col<-c()
    for (i in 1:total){
      if (additionalInfo()[i,1]==levels(additionalInfo()[,1])[1]){
        col[i]<-"#66CD00"
      } else{
        col[i]<-"#8B2500"
      }
    }
    
    heatmap(ES_heat,margins=c(5,3.5), scale="none",cexRow = 0.7, cexCol = 0.5, ColSideColors=col)
    legend("left",      # location of the legend on the heatmap plot
           legend = c(levels(additionalInfo()[,1])[1], levels(additionalInfo()[,1])[2]), # category labels
           col = c("#66CD00", "#8B2500"),  # color key
           lty= 1,             # line style
           lwd = 5,
           cex=0.5 # line width
    )
    #print(H)
    dev.off()
  })
  
  output$GOgraphText<-renderText({GOgraph.Text()})
  
  output$GOGraph<-renderPlot({
    GO.graph()
  })
  
  output$Down_Gograph<- downloadHandler(filename = function() {
    paste("AcyclicGraph", Sys.Date(), '.pdf', sep='')
  },
  content = function(file3) {
    pdf(file3,width=38,height=15)
    DEgeneSets<-DEgeneSets()
    
    range<-nrow(DEgeneSets)
    DEgeneSets_terms<-vector()
    for (i in 1:range){
      DEgeneSets_terms[i]<-paste(unlist(strsplit(tolower(rownames(DEgeneSets)[i]),'_')),collapse=' ')
    }
    
    results<-vector()
    for (i in 1:range){
      if (DEgeneSets_terms[[i]]%in%go[,2]){
        A<-go[,1][which(go[,2]==DEgeneSets_terms[[i]])]
        A<-as.character(A)
        results[i]<-A} #GOMFID
      else{
        results[i]<-"Not included"
      }
    }
    
    g<-results[which(results!="Not included")]
    
    if (input$map=="GOMFPARENTS"){
      map<-GOMFPARENTS
    } else{
      map<-GOMFCHILDREN
    }
    g_plot<-GOGraph(g,map)
    plot(g_plot)
    dev.off()
  })

  ############################################################################

#Tab 2: Training Read Across
  
#1. UI outputs____________________________________________
  #Image
  output$image2<-renderUI(if(input$run_tr==FALSE){
    img(src="https://raw.githubusercontent.com/DemetraDanae/toxFlow/master/toxFlow_logo.png", align="center")
  })
  
  #Files
  output$Import_RA_phChem<-renderUI({
    if(input$Files_RA=="Import_RA"){ 
      tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes). The 1st row must contain the toxicity endpoint for each sample.",fileInput("descrRaw","Physicochemical data:", accept = "text/csv"))}
  })
  
  output$Import_RA_bio<-renderUI({
    if(input$Files_RA=="Import_RA"){ 
      tags$div(title="Requiering .csv files. Columns must contain the NPs (or samples) and rows the markers (or genes). The 1st row must contain the toxicity endpoint for each sample.",fileInput("pcoronaRaw","Biological data:", accept = "text/csv"))}
  })
  
  output$instr2_1<-renderText("Please insert data or use the demo dataset")
  
  output$instr2_2<-renderText("Ready to run analysis")
  
  output$instr2_3<-renderText("Please load the corresponding file for visualization")
  
  output$instr2_4<-renderText("Ready for visualization")
  
  output$instr2_5<-renderText("Protein IDs not found in file of GSVAnalysis. Please import another file")
  
  #Physicochemical threshold 1 for netplots
  output$slider_PhChplot1<-renderUI(
    numericInput("PhCh_th1_1"," ",min=0,max=1,value=0.1,step=0.01)
  )
  
  #Physicochemical threshold 2 for netplots
  output$slider_PhChplot2<-renderUI(
    numericInput("PhCh_th1_2"," ",min=input$PhCh_th1_1, max=1,value=0.5,step=0.01)
  )
  
  #Biological threshold 1 for netplots
  output$slider_Bioplot1<-renderUI(
      numericInput("Bio_th1_1","",min=0,max=1,value=0.1,step=0.01)
  )
  
  #Biological threshold 2 for netplots
  output$slider_Bioplot2<-renderUI(
      numericInput("Bio_th1_2"," ",min=input$Bio_th1_1,max=1,value=0.5,step=0.01)
  )
  
  #Button
  output$DownBut1<-renderUI(if(input$run_tr){downloadButton("DownRes", label = "Download")})
  
  output$filtProteins<-renderUI(if(input$run_gsva){
    checkboxInput("DEproteins", "Use of differentially expressed genes or proteins from GSVA analysis")
  })
  
  output$find_uni<-renderUI(
    if (input$run_tr){
      numericInput("ref","Reference nanoparticle:",value=1,min=1,max=total())}
    else{
      numericInput("ref","Reference nanoparticle:",value=1,min=1,max=NA)}
  )
  
  output$size<-renderUI({
    if(input$size){ 
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' diameter.",fileInput("size_file", " ", accept = "text/csv"))}
  })
  
  output$class<-renderUI({
    if(input$class){ 
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' phenotype or other categorical variable of interest.",fileInput("class_file", " ", accept = "text/csv"))}
  })
  
#2. Start cross validation (training)____________________________________________________
 
  #Import files
  observe({
    if (input$Files_RA=="Files2"){
      shinyjs::hide("only_one")
    } else{
      shinyjs::show("only_one")
    }
  })
  
  observe({
    if (input$only_one==FALSE){
      shinyjs::hide("the_file")
    } else{
      shinyjs::show("the_file")
    }
  })
  
  observe({
    if (input$only_one==TRUE){
        if(input$"the_file"=="only_phchem"){
          shinyjs::show("Import_RA_phChem")
          shinyjs::show("scaling_descr")
          shinyjs::show("lvl_phchem")
          shinyjs::show("slider_PhCh")
          shinyjs::show("slider_PhChplot1")
          shinyjs::show("slider_PhChplot2")
          shinyjs::hide("Import_RA_bio")
          shinyjs::hide("scaling_bio")
          shinyjs::hide("lvl_bio")
          shinyjs::hide("slider_Bio")
          shinyjs::hide("slider_Bioplot1")
          shinyjs::hide("slider_Bioplot2")
          shinyjs::hide("RAcorrBase")
        } else {
          shinyjs::hide("Import_RA_phChem")
          shinyjs::hide("scaling_descr")
          shinyjs::hide("lvl_phchem")
          shinyjs::hide("slider_PhCh")
          shinyjs::hide("slider_PhChplot1")
          shinyjs::hide("slider_PhChplot2")
          shinyjs::hide("RAcorrBase")
          shinyjs::show("Import_RA_bio")
          shinyjs::show("scaling_bio")
          shinyjs::show("lvl_bio")
          shinyjs::show("slider_Bio")
          shinyjs::show("slider_Bioplot1")
          shinyjs::show("slider_Bioplot2")}
    }else{
      shinyjs::show("Import_RA_phChem")
      shinyjs::show("scaling_descr")
      shinyjs::show("lvl_phchem")
      shinyjs::show("slider_PhCh")
      shinyjs::show("slider_PhChplot1","slider_PhChplot2")
      shinyjs::show("Import_RA_bio")
      shinyjs::show("scaling_bio")
      shinyjs::show("lvl_bio")
      shinyjs::show("slider_Bio")
      shinyjs::show("slider_Bioplot1","slider_Bioplot2")
      shinyjs::show("RAcorrBase")
    }
  })
  
  descr<-reactive({
     if ((input$"Files_RA"=="Files2")){
      raw1<-read.csv2("phChem_readAcross.csv",header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    } else{
      if ((input$"only_one"==TRUE)&(input$"the_file"=="only_bio")){
        raw1<-pcorona()
      }else {
        file<-read.csv2(input$descrRaw$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
        raw1<-as.data.frame(t(file)) #for new format
        }}
      raw1})
  
  pcorona<-reactive({
     if (input$"Files_RA"=="Files2"){
      raw2<-read.csv2("bio_readAcross.csv",header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    }else{
      if ((input$"only_one"==TRUE)&(input$"the_file"=="only_phchem")){
          raw2<-descr()
      }else{
        file<-read.csv2(input$pcoronaRaw$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
        raw2<-as.data.frame(t(file)) #for new format
        }}
      raw2})

  #Disable run button if files not uploaded 
  observe({
    if (input$Files_RA=="Files2"){
      shinyjs::enable("run_tr")
    } else {
      if (input$only_one==FALSE){
        if ((is.null(input$descrRaw) || input$descrRaw == "") || (is.null(input$pcoronaRaw) || input$pcoronaRaw == "")){
          shinyjs::disable("run_tr")
        } else{
          shinyjs::enable("run_tr")
        }
      }else{
        if ((is.null(input$descrRaw) || input$descrRaw == "") & (is.null(input$pcoronaRaw) || input$pcoronaRaw == "")){
          shinyjs::disable("run_tr")
        } else {
          shinyjs::enable("run_tr")
        }
      }}
    })
  
  #Disable checkbox if gsva not performed previously or if only phchem file available
  observe({
    if ((input$run_gsva == FALSE)||(input$only_one==TRUE & input$the_file=="only_phchem")){
      shinyjs::disable("DEproteins")
    } else {
      shinyjs::enable("DEproteins")}
    })
  
  #Disable visualization if training not performed or if additional files not uploaded
  observe({
    if (input$run_tr == FALSE || (input$size & is.null(input$size_file)) || (input$class & is.null(input$class_file))) {
      shinyjs::disable("viz1")
    } else {
     shinyjs::enable("viz1")
    }
  })
  
  #Messages if files not uploaded
  observe({
    if (input$Files_RA=="Files2"){
      shinyjs::hide("instr2_1")
      shinyjs::show("instr2_2")
      shinyjs::hide("instr2_4")
    } else {
      if (input$only_one==FALSE){
        if ((is.null(input$descrRaw) || input$descrRaw == "") || (is.null(input$pcoronaRaw) || input$pcoronaRaw == "")){
          shinyjs::show("instr2_1")
          shinyjs::hide("instr2_2")
          shinyjs::hide("instr2_4")
        } else{
          shinyjs::hide("instr2_1")
          shinyjs::show("instr2_2")
          shinyjs::hide("instr2_4")
        }
        }else{
          if ((is.null(input$descrRaw) || input$descrRaw == "") & (is.null(input$pcoronaRaw) || input$pcoronaRaw == "")){
            shinyjs::show("instr2_1")
            shinyjs::hide("instr2_2")
            shinyjs::hide("instr2_4")
          } else {
            shinyjs::hide("instr2_1")
            shinyjs::show("instr2_2")
            shinyjs::hide("instr2_4")
        }
      }}
  })
  
  #Messages if files for visualization not uploaded
  observe({
    if ((input$size & is.null(input$size_file)) || (input$class & is.null(input$class_file))){
      shinyjs::show("instr2_3")
      shinyjs::hide("instr2_2")
      shinyjs::hide("instr2_4")
      shinyjs::hide("instr2_5")
    } else {
      shinyjs::show("instr2_4")
      shinyjs::hide("instr2_3")
      shinyjs::hide("instr2_2")
      shinyjs::hide("instr2_5")
    }
  })
  
  rng_descr<-eventReactive(input$run_tr,{nrow(descr())})
  
  rng_pcorona<-eventReactive(input$run_tr,{nrow(pcorona())})
  
  total<-eventReactive(input$run_tr,{
    if (input$Files_RA=="Files2" || input$only_one==FALSE){
      tot<-ncol(descr())
    } else {
    if (input$the_file=="only_phchem"){#((is.null(input$descrRaw)==FALSE & input$descrRaw != "")){
      tot<-ncol(descr())
    } else {
      tot<-ncol(pcorona())
    }}
    tot
   })
  
  nano_ID<-eventReactive(input$run_tr,{
    if (input$Files_RA=="Files2" || input$only_one==FALSE){
      nano_ID<-nanoID_crossval(descr())
      #print(t(descr()[1,]))
      #print(is.data.frame(descr()))
    } else{ 
      if (input$the_file=="only_phchem"){#((is.null(input$descrRaw)==FALSE & input$descrRaw != "")){
        nano_ID<-nanoID_crossval(descr())
      } else {
        nano_ID<-nanoID_crossval(pcorona())
      }}
    nano_ID
    })
  
  descr_sc<- eventReactive(input$run_tr,{
    
    if (input$scaling_descr==TRUE){  
      descr_sc<-scl(descr(),rng_descr(),total())
    } else {
      descr_sc<-data.frame()
      descr_sc<-descr()
    }
    descr_sc})
  
  pcorona_sc<- eventReactive(input$run_tr,{
    
    if (input$scaling_bio==TRUE){
      pcorona_sc<-scl(pcorona(),rng_pcorona(),total())
    } else{
      pcorona_sc<-data.frame()
      pcorona_sc<-pcorona()
    }
    #Differentially expressed proteins
    if (input$DEproteins==TRUE){
      DEproteins<-DEproteins()
      netcell<-pcorona_sc[1,]
      DEproteins_common<-pcorona_sc[which(rownames(pcorona_sc)%in%DEproteins),]
      pcorona_sc<-rbind(netcell,DEproteins_common)
    }
    pcorona_sc
  })
  
  calculations<-reactive({
    
    descr_sc<-descr_sc()
    pcorona_sc<-pcorona_sc()
    
    #Set levels
    if (input$only_one==FALSE){
      level_phch<-input$lvl_phchem
      level_biol<-input$lvl_bio
    } else {
      if (input$the_file=="only_phchem"){
        level_phch<-input$lvl_phchem
        level_biol<-0.005
      } else{
        level_phch<-0.005
        level_biol<-input$lvl_bio
      }
    }
    
    descr_sc<-filt(descr_sc,"pearson",level_phch)
    pcorona_sc<-filt(pcorona_sc,"pearson",level_biol)
    
    relDescr<-relfun(input$calc,descr_sc)
    relBio<-relfun(input$calc,pcorona_sc)
    
    #Set thresholds and prediction base
    if (input$only_one==FALSE){
      thr_phch<-input$PhCh_th
      thr_biol<-input$Bio_th
      base<-input$RAcorrBase
    } else {
      if (input$the_file=="only_phchem"){
        thr_phch<-input$PhCh_th
        base<-"PhChem"
        thr_biol<-0
      } else{
        thr_biol<-input$Bio_th
        base<-"Bio"
        thr_phch<-0
      }
    }
   
    results_full<-data.frame(stringsAsFactors = TRUE)
    results<-data.frame(stringsAsFactors = TRUE)
    
    for (ref in 1:total()){
      selectedData<-selData(thr_phch,thr_biol,input$calc,ref,relDescr,relBio,total(),nano_ID())
      RAV<-RAVfun(base,selectedData,nano_ID())
      results_full[ref,1]<-nano_ID()[ref,1] #nano_number
      results_full[ref,2]<-nano_ID()[ref,2] #nano_name
      results_full[ref,3]<-nano_ID()[ref,3]   #n.cell
      results_full[ref,4]<-RAV
      results_full[ref,5]<-nano_ID()[ref,3]-RAV #residual
    }
    colnames(results_full)<-c("#","Nanoparticle's ID", "Endpoint","Read across value", "Residuals")
    
    results<-na.omit(results_full)
    rownames(results) <- seq(length=nrow(results))
    
    results
  })
  
  Training.Res<-eventReactive(input$run_tr,{"Training results"})
  
#3. Outputs__________________________________________________________  
  
  output$TrainingRes<-renderText({Training.Res()})
  
  output$resultsTable<-renderDataTable({
    display_calc<-calculations()
    display_calc[,4]<-round(x = display_calc[,4],digits = 6) #Display digits of RAV
    display_calc[,5]<-round(x = display_calc[,5],digits = 6)
    display_calc
    }) 
  
  output$DownRes<- downloadHandler(filename = function() {
    paste("Training", Sys.Date(), '.csv', sep='')
  },
  content = function(con) {
    write.csv2(calculations(),con)
  })
  
  Corr.Coeff<-eventReactive(input$run_tr,{"Correlation Coefficient"})
  
  output$CorrCoeff<-renderText({Corr.Coeff()})
  
  output$R2<-renderText({
    
    predict<-calculations()
    
    #R^2 corr_based calculation
    y<-predict[,3]
    yPred<-predict[,4]
    R2<-r2.funct.lm(y,yPred)
    
    R2
    
  })
  
  Nano.Uni<-eventReactive(input$run_tr,{"Nanoparticles' universe"})
  
  output$NanoUni<-renderText({Nano.Uni()})
  
  diametr<-eventReactive(input$size,{
    diametr<-read.csv2(input$size_file$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    diametr})
  
  categ<-eventReactive(input$class,{
    categ<-read.csv2(input$class_file$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    categ})
  
  net<-eventReactive(input$viz1,{
    
    #Set levels
    if (input$only_one==FALSE){
      level_phch<-input$lvl_phchem
      level_biol<-input$lvl_bio
    } else {
      if (input$the_file=="only_phchem"){
        level_phch<-input$lvl_phchem
        level_biol<-0.005
      } else{
        level_phch<-0.005
        level_biol<-input$lvl_bio
      }
    }
    
    descr_sc<-descr_sc()
    pcorona_sc<-pcorona_sc()
    
    descr_sc<-filt(descr_sc,"pearson",level_phch)
    pcorona_sc<-filt(pcorona_sc,"pearson",level_biol)
    
  relDescr<-relfun(input$calc,descr_sc)
  relBio<-relfun(input$calc,pcorona_sc)
  
  #Set thresholds and prediction base
  if (input$only_one==FALSE){
    PhCh1_1<-input$PhCh_th1_1
    PhCh1_2<-input$PhCh_th1_2
    Bio1_1<-input$Bio_th1_1
    Bio1_2<-input$Bio_th1_2
  } else {
    if (input$the_file=="only_phchem"){
      PhCh1_1<-input$PhCh_th1_1
      PhCh1_2<-input$PhCh_th1_2
      Bio1_1<-0
      Bio1_2<-0
    } else{
      Bio1_1<-input$Bio_th1_1
      Bio1_2<-input$Bio_th1_2
      PhCh1_1<-0
      PhCh1_2<-0
    }
  }
  
  selected<-nano_ID()
  range<-nrow(selected)
  id<-c()
  rel<-c()
  for (i in seq(1,range)){
    id[i]<-paste("s0", i, sep="")
  }
  nodes<-cbind.data.frame(id,selected[,1:2])
  colnames(nodes)<-c("id","nano_number","nano_name")
  
  l<-data.frame() #coordinates
  
  r1<-seq(from=8, to=20, by=1.5)
  r2<-seq(from=30, to=43, by=0.8)
  r3<-seq(from=55, to=65, by=0.3)
  
  f<-seq(from=0, to=2*pi, by=0.008*pi)
  
  links<-data.frame(stringsAsFactors = TRUE)
    for (i in seq(1,range)){
      links[i,1]<-id[which(nodes$nano_number==input$ref)]
      links[i,2]<-paste("s0", i, sep="")
      links[i,3]<-"mention"   #type?
      id_ref<-which(nodes$nano_number==input$ref)
      NP<-nodes[i,2]
      if (i==id_ref){ #color_of_nodes
        links[i,4]<-"#FF3030"
        l[i,1]<-0
        l[i,2]<-0}
      else if ((relDescr[input$ref,NP]>=PhCh1_1)&(relBio[input$ref,NP]>=Bio1_1)){
        links[i,4]<-"#008B8B"
        rad<-sample(r1,1)
        theta<-sample(f,1)
        l[i,1]<-rad*cos(theta)
        l[i,2]<-rad*sin(theta)
      }
      else if ((relDescr[input$ref,NP]>=PhCh1_2)&(relDescr[input$ref,NP]<PhCh1_1)&(relBio[input$ref,NP]>=Bio1_2)&(relBio[input$ref,NP]<Bio1_1)) {
        links[i,4]<-"#00CDCD"
        rad<-sample(r2,1)
        theta<-sample(f,1)
        l[i,1]<-rad*cos(theta)
        l[i,2]<-rad*sin(theta)
      }
      else {#if ((relDescr[input$ref,NP]>input$PhCh_th1_2)&(relBio[input$ref,NP]>input$Bio_th1_2)){
        links[i,4]<-"#C6E2FF"
        rad<-sample(r3,1)
        theta<-sample(f,1)
        l[i,1]<-rad*cos(theta)
        l[i,2]<-rad*sin(theta)
      }
  }
  colnames(links)<-c("from", "to", "type","color")#,"groups")
  
  if (input$class){
    col<-c()
    for (i in 1:total()){
      if (categ()[i,1]==levels(categ()[,1])[1]){
        col[i]<-"#66CD00"
      } else{
        col[i]<-"#8B2500"
      }
    }
    links[,4]<-col
  }
  
  net <- graph.data.frame(links, nodes, directed=T)
  #net <- simplify(net, remove.multiple = F, remove.loops = T)
  
  if (input$size){
    sz<-as.matrix(diametr())/2}
  else{
    sz<-15
  }
  
  l<-as.matrix(l)
  plot1<-plot(net, vertex.color=links$color, vertex.label=V(net)$nano_name, vertex.label.cex=0.7,vertex.label.dist =0,
              vertex.size=sz,vertex.frame.color="white",edge.color="white",vertex.label.family="Helvetica",vertex.label.color="black",edge.arrow.size=0.1,
              edge.arrow.width=0.5,edge.label.cex=0.7, margin=c(1,1,1,1),rescale=F,layout = l*0.017)#,
  #title(main=paste(input$ref,nano_ID[input$ref,2]),cex.main=0.5)
  #mark.groups=list(group, c(15:17)), mark.col=c("#C5E5E7","#ECD89A","#00CDCD"), mark.border=NA)
  #par(mai=c(0.5,0.5,0.5,0.5))
  if (input$class){
    legend("left",      # location of the legend on the heatmap plot
           legend = c(levels(categ()[,1])[1], levels(categ()[,1])[2]), # category labels
           col = c("#66CD00", "#8B2500"),  # color key
           lty= 1,             # line style
           lwd = 5,
           cex=0.5 # line width
    )}
  plot1
  })
  
  output$netAll<-renderPlot({
    net()
  })
  
  ############################################################################
  
#Tab 3: Prediction  
  
#1. UI outputs________________________________________________________  
  #Image
  output$image3<-renderUI(if(input$run_pred==FALSE){
    img(src="https://raw.githubusercontent.com/DemetraDanae/toxFlow/master/toxFlow_logo.png", align="center")
  })
  
  output$find_uni_pred<-renderUI(
    if (input$run_pred){
      numericInput("ref_p","Reference nanoparticle:",value=1,min=1,max=ncol(descr_p()))}
    else{
      numericInput("ref_p","Reference nanoparticle:",value=1,min=1,max=NA)}  
  )
  
  output$instr3_1<-renderText("Please insert data")
  
  output$instr3_2<-renderText("Ready to run analysis")
  
  output$instr3_3<-renderText("Please load the corresponding file for visualization")
  
  output$instr3_4<-renderText("Ready for visualization")
  
  output$instr3_5<-renderText("The predictions can be generated after read-across training has been performed")
  
  #Physicochemical threshold 1 for netplot
  output$slider_PhChplot1_p<-renderUI(
      numericInput("PhCh_th1_1p","",min=0,max=1,value=0.1,step=0.01)
  )
  
  #Physicochemical threshold 2 for netplot
  output$slider_PhChplot2_p<-renderUI(
    numericInput("PhCh_th1_2p","",min=input$PhCh_th1_1p,max=1,value=0.5,step=0.01)
  )
  
  #Biological threshold 1 for netplot
  output$slider_Bioplot1_p<-renderUI(
   numericInput("Bio_th1_1p","",min=0,max=1,value=0.1,step=0.01)
  )
  
  #Biological threshold 2 for netplot
  output$slider_Bioplot2_p<-renderUI(
      numericInput("Bio_th1_2p","",min=input$Bio_th1_1p,max=1,value=0.5,step=0.01)
  )
  
  output$size_p<-renderUI({
    if(input$size_p){ 
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' diameter.",fileInput("size_file_p", " ", accept = "text/csv"))}
  })
  
  output$class_p<-renderUI({
    if(input$class_p){ 
      tags$div(title="Requiering .csv file. The 1st column must contain the samples and the 2nd samples' phenotype or other categorical variable of interest.",fileInput("class_file_p", " ", accept = "text/csv"))}
  })
  
  output$DownBut2<-renderUI(
    if(input$run_pred){
    downloadButton("DownPred", label = "Download")})
  
  #2. Start prediction_______________________________________________
  observe({
    if (input$only_one==TRUE){
      if(input$"the_file"=="only_phchem"){
        shinyjs::show("descrRaw_p")
        shinyjs::show("slider_PhChplot1_p")
        shinyjs::show("slider_PhChplot2_p")
        shinyjs::hide("pcoronaRaw_p")
        shinyjs::hide("slider_Bioplot1_p")
        shinyjs::hide("slider_Bioplot2_p")
      } else {
        shinyjs::hide("descrRaw_p")
        shinyjs::hide("slider_PhChplot1_p")
        shinyjs::hide("slider_PhChplot2_p")
        shinyjs::show("pcoronaRaw_p")
        shinyjs::show("slider_Bioplot1_p")
        shinyjs::show("slider_Bioplot2_p")}
    }else{
      shinyjs::show("descrRaw_p")
      shinyjs::show("slider_PhChplot1_p","slider_PhChplot2_p")
      shinyjs::show("pcoronaRaw_p")
      shinyjs::show("slider_Bioplot1_p","slider_Bioplot2_p")
     }
  })
  
  descr_p<-reactive({
    if ((input$"only_one"==TRUE)&(input$"the_file"=="only_bio")){
      raw1<-descr()
    }else {
    raw1<-read.csv2(input$descrRaw_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    raw1<-t(raw1)}#for new format
    raw1})
  
  pcorona_p<-reactive({
    if ((input$"only_one"==TRUE)&(input$"the_file"=="only_phchem")){
      raw2<-pcorona()
    }else{
    raw2<-read.csv2(input$pcoronaRaw_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    raw2<-t(raw2)}#for new format
  raw2})
  
  #Disable run button if files not uploaded and training not performed
  observe({
    if (input$run_tr==FALSE){
      shinyjs::disable("run_pred")
    } else {
      if (input$only_one==FALSE){
        if ((is.null(input$descrRaw_p) || input$descrRaw_p == "") || (is.null(input$pcoronaRaw_p) || input$pcoronaRaw_p == "")){
          shinyjs::disable("run_pred")
        } else{
          shinyjs::enable("run_pred")
        }
      } else{
        if ((is.null(input$descrRaw_p) || input$descrRaw_p == "") & (is.null(input$pcoronaRaw_p) || input$pcoronaRaw_p == "")){
          shinyjs::disable("run_pred")
        } else{
          shinyjs::enable("run_pred")
        }
      }
    }
    })
  
  #Disable visualization if training not performed or if additional files not uploaded
  observe({
    if (input$run_pred == FALSE || (input$size_p & is.null(input$size_file_p)) || (input$class_p & is.null(input$class_file_p))) {
      shinyjs::disable("viz2")
    } else {
      shinyjs::enable("viz2")
    }
  })
  
  #Messages if files not uploaded
  observe({
    if (input$run_tr==FALSE){
      shinyjs::hide("instr3_1")
      shinyjs::hide("instr3_2")
      shinyjs::hide("instr3_3")
      shinyjs::hide("instr3_4")
      shinyjs::show("instr3_5")
    } else {
      if (input$only_one==FALSE){
        if ((is.null(input$descrRaw_p) || input$descrRaw_p == "") || (is.null(input$pcoronaRaw_p) || input$pcoronaRaw_p == "")){
          shinyjs::show("instr3_1")
          shinyjs::hide("instr3_2")
          shinyjs::hide("instr3_3")
          shinyjs::hide("instr3_4")
          shinyjs::hide("instr3_5")
        } else{
          shinyjs::hide("instr3_1")
          shinyjs::show("instr2_2")
          shinyjs::hide("instr3_3")
          shinyjs::hide("instr3_4")
          shinyjs::hide("instr3_5")
        }
      } else{
        if ((is.null(input$descrRaw_p) || input$descrRaw_p == "") & (is.null(input$pcoronaRaw_p) || input$pcoronaRaw_p == "")){
          shinyjs::show("instr3_1")
          shinyjs::hide("instr3_2")
          shinyjs::hide("instr3_3")
          shinyjs::hide("instr3_4")
          shinyjs::hide("instr3_5")
        } else{
          shinyjs::hide("instr3_1")
          shinyjs::show("instr2_2")
          shinyjs::hide("instr3_3")
          shinyjs::hide("instr3_4")
          shinyjs::hide("instr3_5")
        }
      }
    }
  })
  
  #Messages if files for visualization not uploaded
  observe({
    if (input$run_pred==TRUE){
      if((input$size_p & is.null(input$size_file_p)) || (input$class_p & is.null(input$class_file_p))){
      shinyjs::hide("instr3_1")
      shinyjs::hide("instr3_2")
      shinyjs::show("instr3_3")
      shinyjs::hide("instr3_4")
      shinyjs::hide("instr3_5")
      } else {
      shinyjs::hide("instr3_1")
      shinyjs::hide("instr2_2")
      shinyjs::hide("instr3_3")
      shinyjs::show("instr3_4")
      shinyjs::hide("instr3_5")
    }
  }})
  
  unknown<-reactive({ #unknown for prediction
    if (input$only_one==FALSE){
      ncol(descr_p())
    }else{
      if (input$the_file=="only_phchem"){
        ncol(descr_p())
      } else {
        ncol(pcorona_p())
      }
    }
  }) 
  
  predCalc<-eventReactive(input$run_pred,{
    
    withProgress(message = "Processing...", {
      
      prediction_full<-data.frame(stringsAsFactors = TRUE)
      
      #Set levels
      if (input$only_one==FALSE){
        level_phch<-input$lvl_phchem
        level_biol<-input$lvl_bio
      } else {
        if (input$the_file=="only_phchem"){
          level_phch<-input$lvl_phchem
          level_biol<-0.005
        } else{
          level_phch<-0.005
          level_biol<-input$lvl_bio
        }
      }
      
      descr_sc<-filt(descr_sc(),"pearson",level_phch)
      pcorona_sc<-filt(pcorona_sc(),"pearson",level_biol)
      
      for (predRef in 1:unknown()){
        
        incProgress(1/unknown(), detail = paste("Calculating sample:", predRef))
        
        #Connect with filtered base
        if (input$only_one==FALSE){
        descrPrediction<-cbind(descr_p()[,predRef],descr()[-1,])
        colnames(descrPrediction)[1]<-colnames(descr_p())[predRef]
        descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
         
        pcoronaPrediction<-cbind(pcorona_p()[,predRef],pcorona()[-1,])
        colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[predRef]
        pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
        tot<-ncol(descrPrediction)
        } else {
          if (input$the_file=="only_phchem"){
            descrPrediction<-cbind(descr_p()[,predRef],descr()[-1,])
            colnames(descrPrediction)[1]<-colnames(descr_p())[predRef]
            descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
            
            tot<-ncol(descrPrediction)
            
            empty<-rnorm(nrow(pcorona_p())-1,0,1)
            pcoronaPrediction<-cbind(empty,pcorona()[-1,])
            colnames(pcoronaPrediction)[1]<-"empty"
            pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
        } else {
            empty<-rnorm(nrow(descr_p())-1,0,1)
            descrPrediction<-cbind(empty,descr()[-1,])
            colnames(descrPrediction)[1]<-"empty"
            descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
            
            pcoronaPrediction<-cbind(pcorona_p()[,predRef],pcorona()[-1,])
            colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[predRef]
            pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
          
            tot<-ncol(pcoronaPrediction)
          }
        }
        
        rng_descrPrediction<-nrow(descrPrediction)
        rng_pcoronaPrediction<-nrow(pcoronaPrediction)
        
        if (input$scaling_descr==TRUE){  
          descrPrediction_sc<-scl(descrPrediction,rng_descrPrediction,tot)
        } else {
          descrPrediction_sc<-data.frame()
          descrPrediction_sc<-descrPrediction
        }
        
        if (input$scaling_bio==TRUE){
          pcoronaPrediction_sc<-scl(pcoronaPrediction,rng_pcoronaPrediction,tot)
        } else{
          pcoronaPrediction_sc<-data.frame()
          pcoronaPrediction_sc<-pcoronaPrediction
        }
        
        #Differentially expressed proteins
        if (input$DEproteins==TRUE){
          DEproteins<-DEproteins()
          #netcell<-pcorona[1,]
          DEproteins_common<-pcoronaPrediction_sc[which(rownames(pcoronaPrediction_sc)%in%DEproteins),]
          #pcoronaPrediction_sc<-rbind(netcell,DEproteins_common)
          pcoronaPrediction_sc<-DEproteins_common
          }
        
        relDescr_p<-relfun(input$calc,descrPrediction_sc)
        relBio_p<-relfun(input$calc,pcoronaPrediction_sc)
      
        IDs<-data.frame()
        if (input$only_one==FALSE){
          tab<-descrPrediction_sc
        } else {
          if (input$the_file=="only_phchem"){
            tab<-descrPrediction_sc
          }else{
            tab<-pcoronaPrediction_sc
          }
        }
        
        for (i in 1:tot){
          IDs[i,1]<-nano_ID()[colnames(tab)[i],"#"]
          IDs[i,2]<-colnames(tab)[i]
          IDs[i,3]<-nano_ID()[colnames(tab)[i],"net.cell"]
        }
        IDs[1,1]<-0
        IDs[1,3]<-0
        colnames(IDs)<-c("#", "nano_name", "net.cell")
        
       #Set thresholds and prediction base
        if (input$only_one==FALSE){
          thr_phch<-input$PhCh_th
          thr_biol<-input$Bio_th
          base<-input$RAcorrBase
        } else {
          if (input$the_file=="only_phchem"){
            thr_phch<-input$PhCh_th
            base<-"PhChem"
            thr_biol<-min(relBio_p)
          } else{
            thr_biol<-input$Bio_th
            base<-"Bio"
            thr_phch<-min(relDescr_p)
          }
        }
        
        selectedData<-selData(thr_phch,thr_biol,input$calc,1,relDescr_p,relBio_p,tot,IDs)
        
        #Calculating Read Across Value
        if (base=="PhChem"){
          col<-"PhChem_corr"
        }
        if (base=="Bio") {
          col<-"Bio_corr"
        }
        
        if (nrow(selectedData)==1){
          RAV<-"loose restrictions for succesful prediction"
        } else {
        RAV<-0
        range<-nrow(selectedData)
        denom<-sum(selectedData[2:range,col]) #denominator
        for (i in 2:range){
          RAV<-RAV+(selectedData[i,col]*IDs[which(IDs$nano_name==selectedData$nano_name[i]),3])/denom
        }}
        prediction_full[predRef,1]<-IDs[1,1] #nano_number
        prediction_full[predRef,2]<-IDs[1,2] #nano_name
        prediction_full[predRef,3]<-RAV
      }
      
      prediction_full[,1]<-c(1:unknown())
      colnames(prediction_full)<-c("#","Nanoparticle's ID","Read across value")
      
      prediction_full
      })})
 
#3. Outputs________________________________________________________ 
 
  Pred.Tab<-eventReactive(input$run_pred,{"Prediction Table"})
  
  output$PredTab<-renderText({
    Pred.Tab()
  })
  
  output$predTable<-renderDataTable({
    predCalc()
   })
  
  output$DownPred<- downloadHandler(filename = function() {
    paste("Prediction", Sys.Date(), '.csv', sep='')
  },
  content = function(file) {
    write.csv2(predCalc(),file)
  })
  
  diametr_p<-eventReactive(input$size,{
    diametr_p<-read.csv2(input$size_file_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    diametr_p})
  
  categ_p<-eventReactive(input$class,{
    categ_p<-read.csv2(input$class_file_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    categ_p})
  
  netpred<-eventReactive(input$viz2,{
    
    #Set levels
    if (input$only_one==FALSE){
      level_phch<-input$lvl_phchem
      level_biol<-input$lvl_bio
    } else {
      if (input$the_file=="only_phchem"){
        level_phch<-input$lvl_phchem
        level_biol<-0.005
      } else{
        level_phch<-0.005
        level_biol<-input$lvl_bio
      }
    }
    
    descr_sc<-filt(descr_sc(),"pearson",level_phch)
    pcorona_sc<-filt(pcorona_sc(),"pearson",level_biol)
    
    if (input$only_one==FALSE){
      descrPrediction<-cbind(descr_p()[,input$ref_p],descr()[-1,])
      colnames(descrPrediction)[1]<-colnames(descr_p())[input$ref_p]
      descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
      
      pcoronaPrediction<-cbind(pcorona_p()[,input$ref_p],pcorona()[-1,])
      colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[input$ref_p]
      pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
      tot<-ncol(descrPrediction)
    } else {
      if (input$the_file=="only_phchem"){
        descrPrediction<-cbind(descr_p()[,input$ref_p],descr()[-1,])
        colnames(descrPrediction)[1]<-colnames(descr_p())[input$ref_p]
        descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
        
        tot<-ncol(descrPrediction)
        
        empty<-rnorm(nrow(pcorona_p())-1,0,1)
        pcoronaPrediction<-cbind(empty,pcorona()[-1,])
        colnames(pcoronaPrediction)[1]<-"empty"
        pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
      } else {
        empty<-rnorm(nrow(descr_p())-1,0,1)
        descrPrediction<-cbind(empty,descr()[-1,])
        colnames(descrPrediction)[1]<-"empty"
        descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),] #filtered attributes
        
        pcoronaPrediction<-cbind(pcorona_p()[,input$ref_p],pcorona()[-1,])
        colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[input$ref_p]
        pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
        
        tot<-ncol(pcoronaPrediction)
      }
    }
    
    rng_descrPrediction<-nrow(descrPrediction)
    rng_pcoronaPrediction<-nrow(pcoronaPrediction)
    
    if (input$scaling_descr==TRUE){  
      descrPrediction_sc<-scl(descrPrediction,rng_descrPrediction,tot)
    } else {
      descrPrediction_sc<-data.frame()
      descrPrediction_sc<-descrPrediction
    }
    
    if (input$scaling_bio==TRUE){
      pcoronaPrediction_sc<-scl(pcoronaPrediction,rng_pcoronaPrediction,tot)
    } else{
      pcoronaPrediction_sc<-data.frame()
      pcoronaPrediction_sc<-pcoronaPrediction
    }
    
    #Differentially expressed proteins
    if (input$DEproteins==TRUE){
      DEproteins<-DEproteins()
      #netcell<-pcorona[1,]
      DEproteins_common<-pcoronaPrediction_sc[which(rownames(pcoronaPrediction_sc)%in%DEproteins),]
      #pcoronaPrediction_sc<-rbind(netcell,DEproteins_common)
      pcoronaPrediction_sc<-DEproteins_common
    }
    
    relDescr_p<-relfun(input$calc,descrPrediction_sc)
    relBio_p<-relfun(input$calc,pcoronaPrediction_sc)
    
    IDs<-data.frame()
    if (input$only_one==FALSE){
      tab<-descrPrediction_sc
    } else {
      if (input$the_file=="only_phchem"){
        tab<-descrPrediction_sc
      }else{
        tab<-pcoronaPrediction_sc
      }
    }
    
    for (i in 1:tot){
      IDs[i,1]<-nano_ID()[colnames(tab)[i],"#"]
      IDs[i,2]<-colnames(tab)[i]
      IDs[i,3]<-nano_ID()[colnames(tab)[i],"net.cell"]
    }
    IDs[1,1]<-0
    IDs[1,3]<-0
    colnames(IDs)<-c("#", "nano_name", "net.cell")
    
    #Set thresholds
    if (input$only_one==FALSE){
      PhCh1_1p<-input$PhCh_th1_1p
      PhCh1_2p<-input$PhCh_th1_2p
      Bio1_1p<-input$Bio_th1_1p
      Bio1_2p<-input$Bio_th1_2p
    } else {
      if (input$the_file=="only_phchem"){
        PhCh1_1p<-input$PhCh_th1_1p
        PhCh1_2p<-input$PhCh_th1_2p
        Bio1_1p<-min(relBio_p)
        Bio1_2p<-min(relBio_p)*2
      } else{
        Bio1_1p<-input$Bio_th1_1p
        Bio1_2p<-input$Bio_th1_2p
        PhCh1_1p<-min(relDescr_p)
        PhCh1_2p<-min(relDescr_p)*2
      }
    }
    
    selected<-IDs
    range<-nrow(selected)
    id<-c()
    rel<-c()
    for (i in seq(1,range)){
      id[i]<-paste("s0", i, sep="")
    }
    nodes<-cbind.data.frame(id,selected[,1:2])
    colnames(nodes)<-c("id","nano_number","nano_name")
    
    l<-data.frame() #coordinates
    
    r1<-seq(from=8, to=20, by=1.5)
    r2<-seq(from=30, to=43, by=0.8)
    r3<-seq(from=55, to=65, by=0.3)
    
    f<-seq(from=0, to=2*pi, by=0.008*pi)
    
    links<-data.frame(stringsAsFactors = TRUE)
    #reference-unknown
    links[1,1]<-paste("s0", 1, sep="")
    links[1,2]<-id[1]
    links[1,3]<-"mention"   #type?
    links[1,4]<-"#FF3030"
    l[1,1]<-0
    l[1,2]<-0
    #rest NPs of training set
      for (i in seq(2,range)){
        links[i,1]<-id[1]
        links[i,2]<-paste("s0", i, sep="")
        links[i,3]<-"mention"   #type?
        NP<-i
        if ((relDescr_p[1,NP]>=PhCh1_1p)&(relBio_p[1,NP]>=Bio1_1p)){
          links[i,4]<-"#008B8B"
          rad<-sample(r1,1)
          theta<-sample(f,1)
          l[i,1]<-rad*cos(theta)
          l[i,2]<-rad*sin(theta)
        }
        else if ((relDescr_p[1,NP]>=PhCh1_2p)&(relDescr_p[1,NP]<PhCh1_1p)&(relBio_p[1,NP]>=Bio1_2p)&(relBio_p[1,NP]<Bio1_1p)){
          links[i,4]<-"#00CDCD"
          rad<-sample(r2,1)
          theta<-sample(f,1)
          l[i,1]<-rad*cos(theta)
          l[i,2]<-rad*sin(theta)
        }
        else {#if ((relDescr_p[1,NP]>input$PhCh_th1_2p)&(relBio_p[1,NP]>input$Bio_th1_2p)){
          links[i,4]<-"#C6E2FF"
          rad<-sample(r3,1)
          theta<-sample(f,1)
          l[i,1]<-rad*cos(theta)
          l[i,2]<-rad*sin(theta)
        }
      }
    colnames(links)<-c("from", "to", "type","color")#,"groups")
    
    if (input$class_p){
      col<-c()
      for (i in 1:total()){
        if (categ_p()[i,1]==levels(categ_p()[,1])[1]){
          col[i]<-"#66CD00"
        } else{
          col[i]<-"#8B2500"
        }
      }
      links[,4]<-col
    }
    
    net <- graph.data.frame(links, nodes, directed=T)
    
    if (input$size_p){
      sz<-as.matrix(diametr_p())/2}
    else{
      sz<-15
    }
    
    l<-as.matrix(l)
    plot1<-plot(net, vertex.color=links$color, vertex.label=V(net)$nano_name, vertex.label.cex=0.7,vertex.label.dist =0,
                vertex.size=sz,vertex.frame.color="white",edge.color="white",vertex.label.family="Helvetica",vertex.label.color="black",edge.arrow.size=0.1,
                edge.arrow.width=0.5,edge.label.cex=0.7, margin=c(1,1,1,1),rescale=F,layout = l*0.017)#,
    
    if (input$class_p){
      legend("left",      # location of the legend on the heatmap plot
             legend = c(levels(categ()[,1])[1], levels(categ()[,1])[2]), # category labels
             col = c("#66CD00", "#8B2500"),  # color key
             lty= 1,             # line style
             lwd = 5,
             cex=0.5 # line width
      )}
    plot1
  })
  
  Nano.Uni.Pred<-eventReactive(input$run_pred,{"Nanoparticles' universe"})
  
  output$NanoUniPred<-renderText({Nano.Uni.Pred()})
  
  output$netAll_p<-renderPlot({
    netpred()
  })
}