#!/usr/bin/env Rscript

# First part is just to parse arguments from command line
library(optparse)
library(dplyr)
library(tidyverse)

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

##### Create output file name using input filename #####

input_filename = tail(strsplit(opt$file, "/")[[1]], n=1)
input_filename_no_extension = tools::file_path_sans_ext(input_filename)

output_filename = paste(file.path(opt$out, input_filename_no_extension), "counts.rds", sep="_")

##### Read Data #####

# select top 10 genotypes with the most number of cells for each cell type

matrix = readRDS(opt$file)
meta = as.data.frame(matrix[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,output_filename)