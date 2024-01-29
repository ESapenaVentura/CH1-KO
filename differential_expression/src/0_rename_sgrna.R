# First part is just to parse arguments from command line
library(optparse)
library(dplyr)
library(tidyverse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="Dataset file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Read file
matrix = readRDS(opt$file)

# rename sgrna
tmp = strsplit(matrix$sgrna, split = "_")
matrix$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})

saveRDS(matrix, file=opt$file)