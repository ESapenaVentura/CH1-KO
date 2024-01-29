# CH1-KO
Challenge 1 - Quantify the size and variation of knock out effect

## Pre-requisites
- Nextflow
- Docker

And that's about it!

## How to run

### General

### Running docker image
1. Build the docker image
   ```
   docker build differential_expression/ "dif-exp:latest"
   ```
   Please note that this step may take some time - While building the container, it will download and install all the
   necessary R libraries with dependencies
2. Run the docker image
   ```shell
   docker run -it dif-exp:latest
   ```
   This will land you in the R command line, where you should be able to run any code with the dependencies already installed.