FROM rocker/r-base:latest

RUN apt-get update && apt-get install -y \
	openjdk-8-jdk \
	libcurl4-openssl-dev \
	libxml2-dev \
	libv8-dev 

RUN ["R", "CMD", "javareconf"]

USER docker

COPY install_dependencies.R /home/docker/toxflow/install_dependencies.R

RUN ["Rscript", "/home/docker/toxflow/install_dependencies.R"]

COPY ./*.R /home/docker/toxflow/
COPY ./www /home/docker/toxflow/

COPY ./Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

ENTRYPOINT ["Rscript", "/home/docker/toxflow/shiny_start.R"]
