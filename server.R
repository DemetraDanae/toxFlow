#Gene Set Collections
setwd("C:/Users/Dimitra/Desktop/Thesis/Shiny/Final App/GSVA files")
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
  
  #UI outputs__________________________________________________________________________
  #images
  output$image1<-renderUI(if(input$run_gsva==FALSE){
    img(src="logo5.png", align="center")
  })
  
  output$image2<-renderUI(if(input$run_tr==FALSE){
    img(src="logo5.png", align="center")
  })
  
  output$image3<-renderUI(if(input$run_pred==FALSE){
    img(src="logo5.png", align="center")
  })
  
  #For gsva
  tags$div(title="The file must be a .csv. The first column must contain GOterms and the second the ENTREZids",output$Other<-renderUI({
    if(input$gsc=="Other"){
      fileInput("other", " ", accept = "text/csv")}
  }))
  
  output$DownGSVA<-renderUI(if(input$run_gsva){downloadButton("Down_gsva", label = "Download")})
  
  #For cross validation
  output$slider_PhCh<-renderUI(
    if(input$calc=="cos"){
      sliderInput("PhCh_th","Physicochemical threshold:",min=0,max=1,value=0.5,step=0.01)
    } else{
      sliderInput("PhCh_th","Physicochemical threshold:",min=0,max=6,value=1,step=0.05)
    }
  )
  
  output$slider_Bio<-renderUI(
    if(input$calc=="cos"){
      sliderInput("Bio_th","Biological threshold:",min=0,max=1,value=0.5,step=0.01)
    } else{
      sliderInput("Bio_th","Biological threshold:",min=0,max=6,value=1,step=0.05)
    }
  )
  
  output$slider_PhChplot<-renderUI(
    if(input$calc=="cos"){
      sliderInput("PhCh_th1","Physicochemical thresholds:",min=0,max=1,value=c(0.1,0.5),step=0.01)
    } else{
      sliderInput("PhCh_th1","Physicochemical thresholds:",min=0,max=6,value=c(3,5),step=0.05)
    }
  )
  
  output$slider_Bioplot<-renderUI(
    if(input$calc=="cos"){
      sliderInput("Bio_th1","Biological thresholds:",min=0,max=1,value=c(0.1,0.5),step=0.01)
    } else{
      sliderInput("Bio_th1","Biological thresholds:",min=0,max=6,value=c(3,5),step=0.05)
    }
  )
  
  output$DownBut1<-renderUI(if(input$run_tr){downloadButton("DownRes", label = "Download")})
  
  #For prediction
  output$slider_PhChplot_p<-renderUI(
    if(input$calc=="cos"){
      sliderInput("PhCh_th1_p","Physicochemical thresholds:",min=0,max=1,value=c(0.1,0.5),step=0.01)
    } else{
      sliderInput("PhCh_th1_p","Physicochemical thresholds:",min=0,max=6,value=c(3,5),step=0.05)
    }
  )
  
  output$slider_Bioplot_p<-renderUI(
    if(input$calc=="cos"){
      sliderInput("Bio_th1_p","Biological thresholds:",min=0,max=1,value=c(0.1,0.5),step=0.01)
    } else{
      sliderInput("Bio_th1_p","Biological thresholds:",min=0,max=6,value=c(3,5),step=0.05)
    }
  )
  
  output$DownBut2<-renderUI(if(input$run_pred){downloadButton("DownPred", label = "Download")})
  
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
  
  filt<-function(X_sc,method){ #Filtering data X_sc
    
    X_t<-t(X_sc) #transpose data matrix
    
    #Calculating correlations between n.cell and factors
    corrNcellDescr<-cor(X_t, y=NULL, use="everything", method = method)
    
    range<-nrow(corrNcellDescr)
    j<-c()
    for (i in 2:range){
      if ((is.na(corrNcellDescr[i,1])==TRUE)|(corrNcellDescr[i,1]<0)){
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
    }else{
      #Calculating distance between nanoparticles
      X1<-t(X_sc)
      relX<-dist(X1, method = calc, diag = FALSE)
      relX<-as.matrix(relX)
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
      RAV<-RAV+(selData[i,col]*nano_ID[selData$"nano_name"[i],3])/denom
    }
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
    if (calc=="cos"){
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
    } else if ((calc=="manhattan")|(calc=="euclidean")){
      while (NP<=total) {
        if (NP!=ref){
          if ((relDescr[NP,ref]<=PhCh_th)&(relBio[NP,ref]<=Bio_th)){  #Thresholds
            selData[i,1]<-nano_ID[NP,1] #nano_number
            selData[i,2]<-nano_ID[NP,2] 
            selData[i,3]<-relDescr[NP,ref] #corr_factor
            selData[i,4]<-relBio[NP,ref] #corr_factor (pcorona)
            i<-i+1
          }
        }
        NP<-NP+1
      }
    }
    
    colnames(selData)<-c("#","nano_name", "PhChem_corr","Bio_corr")
    
    return(selData)}
  
  netplot<-function(ref,nano_ID,PhCh_th1,PhCh_th2,Bio_th1,Bio_th2,relDescr,relBio){
    #Network diagram
    selected<-nano_ID
    range<-nrow(selected)
    id<-c()
    rel<-c()
    for (i in seq(1,range)){
      id[i]<-paste("s0", i, sep="")
    }
    nodes<-cbind.data.frame(id,selected[,1:2])
    colnames(nodes)<-c("id","nano_number","nano_name")
    
    links<-data.frame(stringsAsFactors = TRUE)
    for (i in seq(1,range)){
      links[i,1]<-id[which(nodes$nano_number==ref)]
      links[i,2]<-paste("s0", i, sep="")
      links[i,3]<-"mention"   #type?
      id_ref<-which(nodes$nano_number==ref)
      NP<-nodes[i,2]
      if (i==id_ref){ #color_of_nodes
        links[i,4]<-"#FFA500" }
      else if ((relDescr[ref,NP]<=PhCh_th2)&(relDescr[ref,NP]>PhCh_th1)&(relBio[ref,NP]<=Bio_th2)&(relBio[ref,NP]>Bio_th1)) {
        links[i,4]<-"#00CDCD"
      }
      else if ((relDescr[ref,NP]>PhCh_th2)&(relBio[ref,NP]>Bio_th2)){
        links[i,4]<-"#008B8B"
      }
      else{
        links[i,4]<-"#00FFFF"
      }
      #if ((relDescr[ref,NP]<=PhCh_th2)&(relDescr[ref,NP]>PhCh_th1)&(relBio[ref,NP]<=Bio_th2)&(relBio[ref,NP]>Bio_th1)){
      #links[i,5]<-1
      # } else if ((relDescr[ref,NP]>PhCh_th2)&(relBio[ref,NP]>Bio_th2)){
      #   links[i,5]<-2
      #} else {
      #   links[i,5]<-3
      #}
    }
    colnames(links)<-c("from", "to", "type","color")#,"groups")
    
    #group1<-links[which(links[,5]==1),2]
    #group2<-links[which(links[,5]==2),2]
    #group3<-links[which(links[,5]==3),2]
    
    net <- graph.data.frame(links, nodes, directed=T)
    net.bg <- barabasi.game(range)
    l <- layout.circle(net.bg)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    #net <- simplify(net, remove.multiple = F, remove.loops = T)
    
    E(net)$width <- 5^E(net)$weight
    
    plot1<-plot(net, vertex.color=links$color, vertex.label=V(net)$nano_name, vertex.label.cex=0.7,vertex.label.dist =0,
                vertex.size=20,vertex.frame.color="white",edge.color="white",vertex.label.family="Helvetica",vertex.label.color="black",edge.arrow.size=0.1,
                edge.arrow.width=0.5,edge.label.cex=0.7,margin=0.005, main=paste(ref,nano_ID[ref,2]), margin=c(5,10,5,10))#,
    #mark.groups=list(group1,group2,group3), mark.col=c("#C5E5E7","#ECD89A","#00CDCD"), mark.border=NA)
    
    return(plot1)}
  
  
  #For GSVA
  #import files
  rawData<-reactive({
    raw<-read.csv2(input$rawData$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    raw})
  
  additionalInfo<-reactive({read.csv2(input$classifData$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)})
  
  #double values?
  uniqueData<-reactive({
    if (input$double==TRUE){  
      X1<- numeric(dim(rawData())[2])
      rownames.X1<- character(0)
      for(i in 1: dim(rawData())[1]){
        r<- colnames(rawData())[i]
        if(length(intersect(r,rownames.X1))==0){ 
          r.uni<- which(colnames(rawData())%in%r)
          X1.uni<- rawData()[r.uni,]
          if(length(r.uni)>1){X1<- rbind(X1,apply(X1.uni,2,mean))}else{X1<- rbind(X1,X1.uni)}
          rownames.X1[i]<- r
        }	
      }
      
      X1<- X1[-1,] 
      rownames.X1<- rownames.X1[!is.na(rownames.X1)]
      rownames(X1)<- rownames.X1
    } else {
      X1<-rawData()
    }
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
    min<-input$geneLimits[1]
    max<-input$geneLimits[2]
    gsva_score<-gsva(as.matrix(X1),geneSet,method="gsva",rnaseq=FALSE,abs.ranking=FALSE, min.sz=min,max.sz=max,no.bootstraps=input$boot,bootstrap.percent =.632,mx.diff=TRUE,verbose=TRUE)
    ES<-gsva_score$es.obs
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
        
        colnames(results)<-c("GOMFID","Term","p-value","Size","Count")
        
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
    
    heatmap(ES_heat,margins=c(6.5,3), scale="none",cexRow = 0.7, cexCol = 0.5, ColSideColors=col)
    legend("topleft",      # location of the legend on the heatmap plot
           legend = c(levels(additionalInfo()[,1])[1], levels(additionalInfo()[,1])[2]), # category labels
           col = c("#66CD00", "#8B2500"),  # color key
           lty= 1,             # line style
           lwd = 5,
           cex=0.5 # line width
    )
    
    #heatmap(ES_heat,margins=c(6.5,3), scale="none",cexRow = 0.7, cexCol = 0.5)
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
  
  #Outputs___________________________________________________________________________
  output$resultsText<-renderText({
    results.Text()
  })
  
  output$resultsTable_gsva<-renderDataTable({
    #X1terms()
    results()
    #DEgeneSets()
  })
  
  output$Down_gsva<- downloadHandler(filename = function() {
    paste("GSVA_analysis", Sys.Date(), '.csv', sep='')
  },
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
    
    colnames(results)<-c("GOMFID","Term","p-value","Size","Count")
    
    write.csv2(results,file1)
  })
  
  output$heatmapText<-renderText({heatmap.Text()})
  
  output$heatmap<-renderPlot({
    heat.map()
  })
  
  output$GOgraphText<-renderText({GOgraph.Text()})
  
  output$GOGraph<-renderPlot({
    GO.graph()
  })
  
  #For cross validation (training)____________________________________________________
  descr<-reactive({
    raw1<-read.csv2(input$descrRaw$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    raw1})
  pcorona<-reactive({
    raw2<-read.csv2(input$pcoronaRaw$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
  raw2})
  
  rng_descr<-eventReactive(input$run_tr,{nrow(descr())})
  rng_pcorona<-eventReactive(input$run_tr,{nrow(pcorona())})
  total<-eventReactive(input$run_tr,{ncol(descr())})
  nano_ID<-eventReactive(input$run_tr,{nanoID_crossval(descr())})
  
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
    
    descr_sc<-filt(descr_sc,"pearson")
    pcorona_sc<-filt(pcorona_sc,"pearson")
    
    relDescr<-relfun(input$calc,descr_sc)
    relBio<-relfun(input$calc,pcorona_sc)
    
    results_full<-data.frame(stringsAsFactors = TRUE)
    results<-data.frame(stringsAsFactors = TRUE)
    
    for (ref in 1:total()){
      selectedData<-selData(input$PhCh_th,input$Bio_th,input$calc,ref,relDescr,relBio,total(),nano_ID())
      RAV<-RAVfun(input$RAcorrBase,selectedData,nano_ID())
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
  
  output$TrainingRes<-renderText({Training.Res()})
  
  output$resultsTable<-renderDataTable({
    calculations()
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
  
 
  net<-eventReactive(input$viz1,{descr_sc<-filt(descr_sc(),"pearson")
  pcorona_sc<-filt(pcorona_sc(),"pearson")
  relDescr<-relfun(input$calc,descr_sc)
  relBio<-relfun(input$calc,pcorona_sc)
  
  netplot(input$ref,nano_ID(),input$PhCh_th1[1],input$PhCh_th1[2],input$Bio_th1[1],input$Bio_th1[2],relDescr,relBio)
  })
  
   output$netAll<-renderPlot({
    net()
  })
  
  #For prediction__________________________________________________________________________
  descr_p<-reactive({
    raw1<-read.csv2(input$descrRaw_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
    raw1})
  pcorona_p<-reactive({raw2<-read.csv2(input$pcoronaRaw_p$datapath,header=TRUE,quote = "\"",dec=".", stringsAsFactors = TRUE,strip.white = TRUE,row.names=1)
  raw2})
  
  unknown<-reactive({ncol(descr_p())}) #unknown for prediction
  
  predCalc<-eventReactive(input$run_pred,{
    
    withProgress(message = "Processing...", {
      
      prediction_full<-data.frame(stringsAsFactors = TRUE)
      
      descr_sc<-filt(descr_sc(),"pearson")
      pcorona_sc<-filt(pcorona_sc(),"pearson")
      
      for (predRef in 1:unknown()){
        
        incProgress(1/unknown(), detail = paste("Calculating sample:", predRef))
        
        #Connect with filtered base
        descrPrediction<-cbind(descr_p()[,predRef],descr()[-1,])
        colnames(descrPrediction)[1]<-colnames(descr_p())[predRef]
        descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),]
        
        pcoronaPrediction<-cbind(pcorona_p()[,predRef],pcorona()[-1,])
        colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[predRef]
        pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
        
        tot<-ncol(descrPrediction)
        rng_descrPrediction<-nrow(descrPrediction)
        rng_pcoronaPrediction<-nrow(pcoronaPrediction)
        
        if (input$scaling_descr_p==TRUE){  
          descrPrediction_sc<-scl(descrPrediction,rng_descrPrediction,tot)
        } else {
          descrPrediction_sc<-data.frame()
          descrPrediction_sc<-descrPrediction
        }
        
        if (input$scaling_bio_p==TRUE){
          pcoronaPrediction_sc<-scl(pcoronaPrediction,rng_pcoronaPrediction,tot)
        } else{
          pcoronaPrediction_sc<-data.frame()
          pcoronaPrediction_sc<-pcoronaPrediction
        }
        
        #Differentially expressed proteins
        if (input$DEproteins==TRUE){
          DEproteins<-DEproteins()
          netcell<-pcoronaPrediction_sc[1,]
          DEproteins_common<-pcoronaPrediction_sc[which(rownames(pcoronaPrediction_sc)%in%DEproteins),]
          pcoronaPrediction_sc<-rbind(netcell,DEproteins_common)
          }
        
        relDescr_p<-relfun(input$calc,descrPrediction_sc)
        relBio_p<-relfun(input$calc,pcoronaPrediction_sc)
        
        IDs<-data.frame()
        for (i in 1:tot){
          IDs[i,1]<-nano_ID()[colnames(descrPrediction_sc)[i],"#"]
          IDs[i,2]<-colnames(descrPrediction_sc)[i]
          IDs[i,3]<-nano_ID()[colnames(descrPrediction_sc)[i],"net.cell"]
        }
        IDs[1,1]<-0
        IDs[1,3]<-0
        colnames(IDs)<-c("#", "nano_name", "net.cell")
        
        selectedData<-selData(input$PhCh_th,input$Bio_th,input$calc,1,relDescr_p,relBio_p,tot,IDs)
        #RAV<-RAVfun(RAcorrBase,selectedData,IDs)
        
        #Calculating Read Across Value
        if (input$RAcorrBase=="PhChem"){
          col<-"PhChem_corr"
        }
        if (input$RAcorrBase=="Bio") {
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
        #if (is.na(RAV==TRUE)){
        #  prediction_full[predRef,3]<-"loose restrictions for successful prediction"
        #} else{
          prediction_full[predRef,3]<-RAV
        #}
      }
      
      prediction_full[,1]<-c(1:unknown())
      colnames(prediction_full)<-c("#","Nanoparticle's ID","Read across value")
      
      prediction_full
    })})
  
  Pred.Tab<-eventReactive(input$run_pred,{"Prediction Table"})
  
  output$PredTab<-renderText({Pred.Tab()})
  
  output$predTable<-renderDataTable({
    predCalc()
  })
  
  output$DownPred<- downloadHandler(filename = function() {
    paste("Prediction", Sys.Date(), '.csv', sep='')
  },
  content = function(file) {
    write.csv2(predCalc(),file)
  })
  
  netpred<-eventReactive(input$viz2,{
    descr_sc<-filt(descr_sc(),"pearson")
    pcorona_sc<-filt(pcorona_sc(),"pearson")
    
    #Connect with filtered base
    descrPrediction<-cbind(descr_p()[,input$ref_p],descr()[-1,])
    colnames(descrPrediction)[1]<-colnames(descr_p())[input$ref_p]
    descrPrediction<-descrPrediction[which(row.names(descrPrediction)%in%row.names(descr_sc)),]
    
    pcoronaPrediction<-cbind(pcorona_p()[,input$ref_p],pcorona()[-1,])
    colnames(pcoronaPrediction)[1]<-colnames(pcorona_p())[input$ref_p]
    pcoronaPrediction<-pcoronaPrediction[which(row.names(pcoronaPrediction)%in%row.names(pcorona_sc)),]
    
    tot<-ncol(descrPrediction)
    rng_descrPrediction<-nrow(descrPrediction)
    rng_pcoronaPrediction<-nrow(pcoronaPrediction)
    
    if (input$scaling_descr_p==TRUE){  
      descrPrediction_sc<-scl(descrPrediction,rng_descrPrediction,tot)
    } else {
      descrPrediction_sc<-data.frame()
      descrPrediction_sc<-descrPrediction
    }
    
    if (input$scaling_bio_p==TRUE){
      pcoronaPrediction_sc<-scl(pcoronaPrediction,rng_pcoronaPrediction,tot)
    } else{
      pcoronaPrediction_sc<-data.frame()
      pcoronaPrediction_sc<-pcoronaPrediction
    }
    
    #Differentially expressed proteins
    if (input$DEproteins==TRUE){
      DEproteins<-DEproteins()
      netcell<-pcoronaPrediction_sc[1,]
      DEproteins_common<-pcoronaPrediction_sc[which(rownames(pcoronaPrediction_sc)%in%DEproteins),]
      pcoronaPrediction_sc<-rbind(netcell,DEproteins_common)
    }
    
    relDescr_p<-relfun(input$calc,descrPrediction_sc)
    relBio_p<-relfun(input$calc,pcoronaPrediction_sc)
    
    IDs<-data.frame()
    for (i in 1:tot){
      IDs[i,1]<-i
      IDs[i,2]<-colnames(descrPrediction_sc)[i]
      IDs[i,3]<-nano_ID()[colnames(descrPrediction_sc)[i],"net.cell"]
    }
    IDs[1,3]<-0
    colnames(IDs)<-c("#", "nano_name", "net.cell")
    
    A2<-IDs[1,2]
    IDs[1,2]<-IDs[which(rownames(IDs)==input$ref_p),2]
    IDs[which(rownames(IDs)==input$ref_p),2]<-A2
    
    netplot(input$ref_p,IDs,input$PhCh_th1_p[1],input$PhCh_th1_p[2],input$Bio_th1_p[1],input$Bio_th1_p[2],relDescr_p,relBio_p)
  })
  
  Nano.Uni.Pred<-eventReactive(input$run_pred,{"Nanoparticles' universe"})
  
  output$NanoUniPred<-renderText({Nano.Uni.Pred()})
  
  output$netAll_p<-renderPlot({
    netpred()
  })
}
