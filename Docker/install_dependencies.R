#Universal Bioconductor package installation function
install.bioc <- function(pkg){
	vers <- getRversion()
	if (vers >= "3.6"){
		if (!requireNamespace("BiocManager", quietly = TRUE)) {
			install.packages("BiocManager")
		}
		BiocManager::install(pkg)
	}else{
		if (!requireNamespace("BiocInstaller", quietly = TRUE)){
			source("https://bioconductor.org/biocLite.R")
			biocLite(pkg, suppressUpdates=TRUE)
		}else{
			BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
		}
	}
}

#Install Bioconductor dependencies
bioc_pkgs <- c(
  'limma', 'GSEAbase','GSVA', 'GSVAdata','org.Hs.eg.db','GSEAlm','GO.db','GOstats', 'Rgraphviz'
)
bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
if(length(bioc_pkgs.inst)>0){
	print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
	for(pkg in bioc_pkgs.inst){
		print(paste0("Installing Package:'", pkg, "'..."))
		install.bioc(pkg)
		print("Installed!!!")
	}
}

#Install CRAN dependencies
cran_pkgs <- c(
  'igraph','shiny', 'shinyjs', 'DT', 'shinyFiles', 'rhandsontable', 'base', 'shinythemes', 'network', 'extrafont', 'lattice', 'utils', 'RColorBrewer', 'readr', 'shinyWidgets', 'rgl', 'shinyBS', 'shinyLP', 'rapportools'
)
cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(cran_pkgs.inst)>0){
	print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
	for(pkg in cran_pkgs.inst){
		print(paste0("Installing Package:'", pkg, "'..."))
		install.packages(
			pkg, repo="http://cran.rstudio.org", dependencies=TRUE
		)
		print("Installed!!!")
	}
}
