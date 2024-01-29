# First part is just to parse arguments from command line
library(optparse)
library(dplyr)
library(tidyverse)
library(ggplot2)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
              help="Counts file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read file
matrix = readRDS(opt$file)


plot = D1[[]]%>%
         filter(genotype%in%genotype.day[["D1"]])%>%
         ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
         geom_boxplot()+
         ggtitle("D1")+
         theme_bw()+
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))