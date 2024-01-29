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

# Read file
matrix = readRDS(opt$file)

# Output file name
output_filename = file.path(opt$out, "D18.ctproption.pdf")

# cell type proportions for D18
meta = matrix[[]]%>%filter(celltype!='Unknown')
df <- meta%>%
  group_by(genotype, celltype) %>%
  summarize(n=n())
df1<-meta%>%group_by(genotype) %>%
  summarize(n.genotype=n())
df=df%>%left_join(df1)
df=df%>%mutate(proportion = n/n.genotype)
fig = df%>%ggplot(aes(fill=celltype, y=proportion, x=genotype)) +
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),legend.position = "top")+
  guides(fill = guide_legend(title =''))
pdf(output_filename, width = 6, height = 4)
  fig
dev.off()