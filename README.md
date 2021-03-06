# toxFlow
R code for toxFlow application. toxFlow is a Shiny application that was build under R version 3.2.3. In this repository ui.R and server.R files are uploaded, as well as additional files and a manual.

toxFlow is an application of web tools for Gene Set Variation Analysis (GSVA) and toxicity prediction using read across technique. toxFlow is developed as a project for my Diploma Thesis, in School of Chemical Engineering (National Technical University of Athens, March-September 2016). Additional information about the employed methods and an detailed case study, can be found in the corresponding <a href="https://pubs.acs.org/doi/pdfplus/10.1021/acs.jcim.7b00160">publication</a>; Varsou DD, Tsiliki G, Nymark P, Kohonen P, Grafström R, and Sarimveis H, <i> toxFlow: A Web-Based Application for Read-Across Toxicity Prediction Using Omics and Physicochemical Data</i>, J. Chem. Inf. Model.

This application is hosted in the server of Unit of Process Control & Informatics of School of Chemical Engineering (NTUA), in the following link <a href="https://toxflow.jaqpot.org/"> https://toxflow.jaqpot.org/ </a>. 

You can find a video tutorial on <a href="https://www.youtube.com/watch?v=kGp2PuTiDrg"> YouTube</a>

doi: <a href="https://doi.org/10.5281/zenodo.836713"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.595814.svg" alt="DOI"></a>

You can find the toxflow docker image here: https://hub.docker.com/r/demetradanae/toxflow

# License
This application is released under <a href="https://www.gnu.org/licenses/gpl.html"> GNU General Public License v.3</a>. 
```html
toxFlow GSVA-Read across web tools

Copyright (C) 2017 Dimitra-Danai Varsou

This program is free software: you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
You should have received a copy of the GNU General Public License along with this program.  
If not, see here: http://www.gnu.org/licenses/.

```

# Cite as
If you find toxFlow useful, please cite me! :)

Dimitra Danai Varsou. (2020, March 19). DemetraDanae/toxFlow: toxFlow Release 1.0.0 (Version v1.0.0). Supplement to:
https://github.com/DemetraDanae/toxFlow/releases/tag/v.1.0.0 Zenodo. https://zenodo.org/badge/latestdoi/68043137

# Updates in toxFlow
Feb 19, 2020

-Deletion of "Help" tab and addition of landing page

-Update of "GSVA" part according to bioconductor/gsva v.1.26 package. Visible changes: Number of bootstraps parameter has been removed.

-Userguide update
