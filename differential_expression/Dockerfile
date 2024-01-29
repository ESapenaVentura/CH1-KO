FROM rocker/r-ver:4.3.0


WORKDIR /app

# Install necessary system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        libpng-dev \
        libglpk40


# Install dependencies
RUN R -e "install.packages(c('Seurat', 'dplyr', 'nichenetr', 'ggplot2', 'RANN', 'tidyverse', 'optparse'))"

# The code to run when container is started:
COPY src/ /app