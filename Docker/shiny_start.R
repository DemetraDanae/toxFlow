library(shiny)
library(shinyjs)

runApp(
	appDir="/home/docker/toxflow", 
	port=3838, 
	launch.browser=FALSE, 
	host="0.0.0.0"
)
