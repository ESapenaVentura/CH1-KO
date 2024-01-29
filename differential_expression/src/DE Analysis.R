library(dplyr)
library(Seurat)
library(nichenetr)
library(ggplot2)
library(RANN)
library(tidyverse)

##### Read Data #####

# select top 10 genotypes with the most number of cells for each cell type at each day

D1 = readRDS("Sample_G_merged_annotated.rds")
D3 = readRDS("Sample_H_annotated.rds")
D7 = readRDS("Sample_I_annotated.rds")
D11 = readRDS("Sample_J_annotated.rds")
D18 = readRDS("Sample_L_merged_annotated.rds")

meta = as.data.frame(D1[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,"Analysis/D1_ctgt_counts.RDS")

meta = as.data.frame(D3[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,"Analysis/D3_ctgt_counts.RDS")

meta = as.data.frame(D7[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,"Analysis/D7_ctgt_counts.RDS")

meta = as.data.frame(D11[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,"Analysis/D11_ctgt_counts.RDS")


meta = as.data.frame(D18[[]])
n.df = as.data.frame(meta%>%dplyr::group_by(celltype,genotype)%>%dplyr::summarise(n()))
saveRDS(n.df,"Analysis/D18_ctgt_counts.RDS")


# cell type proportions for D18
meta = D18[[]]%>%filter(celltype!='Unknown')
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
pdf("Analysis/D18.ctproption.pdf", width = 6, height = 4)
  fig
dev.off()



# rename sgrna
tmp = strsplit(D1$sgrna, split = "_")
D1$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})

tmp = strsplit(D3$sgrna, split = "_")
D3$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})

tmp = strsplit(D7$sgrna, split = "_")
D7$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})

tmp = strsplit(D11$sgrna, split = "_")
D11$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})

tmp = strsplit(D18$sgrna, split = "_")
D18$sgrna2=sapply(tmp,function(x){
  paste0(x[2],"-",x[1])
})


###### check sequencing depth ######

# similar sequencing depth across celltype-genotype combinations in most cases

pdf("Analysis/nCount_RNA.pdf")
D1[[]]%>%
  filter(genotype%in%genotype.day[["D1"]])%>%
  ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
  geom_boxplot()+
  ggtitle("D1")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
D3[[]]%>%
  filter(genotype%in%genotype.day[["D3"]])%>%
  ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
  geom_boxplot()+
  ggtitle("D3")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
D7[[]]%>%
  filter(genotype%in%genotype.day[["D7"]])%>%
  ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
  geom_boxplot()+
  ggtitle("D7")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
D11[[]]%>%
  filter(genotype%in%genotype.day[["D11"]])%>%
  ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
  geom_boxplot()+
  ggtitle("D11")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
D18[[]]%>%
  filter(genotype%in%genotype.day[["D18"]])%>%
  ggplot(aes(x=sgrna2,y=nCount_RNA,color = genotype))+
  geom_boxplot()+
  ggtitle("D18")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


##### Cell-Level WilCoxon test  #####



cell_DE<-function(obj,ct.list,geno2test,min_cell=30){
  DE.list = list()
  for(ct in ct.list){
    print(ct)
    DE.list[[ct]] = list()
    for(g in geno2test){
      if(!paste(ct,g,sep='_')%in%levels(as.factor(obj$celltype.genotype))){
        DE.list[[ct]][[g]]=NA
      }else{
        n1 = (obj[[]]%>%count(celltype.genotype)%>%filter(celltype.genotype==paste(ct,g,sep='_')))$n
        n0 = (obj[[]]%>%count(celltype.genotype)%>%filter(celltype.genotype==paste(ct,"WT",sep='_')))$n
        if(n1<min_cell){
          DE.list[[ct]][[g]]=NA
        }else{
          res = FindMarkers(object = obj, 
                            ident.1 =  paste(ct,g,sep='_'),
                            ident.2 = paste(ct,"WT",sep='_'),
                            test.use = "wilcox",
                            min.pct=max(0.1,10/(min(n1,n0))))
          DE.list[[ct]][[g]]<-res%>%dplyr::arrange(p_val)
        }
      }
    }
  }
  return(DE.list)
}

get_prop_sig<-function(DE,max.p=0.05){ # get proportion of significant DE genes
  prop_sig=lapply(DE,function(df.list){
    sapply(df.list,function(df){
      if(all(is.na(df))){
        return(NA)
      }
      return(sum(df$p_val<max.p=0.05)/nrow(df))
    })
  })
  return(prop_sig)
}

wilcoxon.path = "Analysis/DE/wilcoxon/"
geno2test = setdiff(levels(as.factor(D1$genotype)),"WT")

###### D1 ######

obj = D1
ct.list=c("ESC")
obj$celltype.genotype <- paste(obj$celltype, obj$genotype, sep = "_")
Idents(obj) <- "celltype.genotype"
DE = cell_DE(obj,ct.list,geno2test)
saveRDS(DE,"Analysis/DE/wilcoxon/D1_DE.RDS")


# check proportions of DE vs sample size
prop.sig=get_prop_sig(DE)
obj.sub = obj%>%subset(celltype%in%ct.list)
n.counts = obj.sub[[]]%>%count(celltype,genotype)
cor(prop.sig[[ct]],(n.counts%>%filter(celltype ==ct&genotype!="WT"))$n,use = 'complete.obs') # 0.5858891



###### D3 ######

obj = D3
summary(as.factor(obj$celltype))
ct.list=c("DE_1","DE_2")
obj$celltype.genotype <- paste(obj$celltype, obj$genotype, sep = "_")
Idents(obj) <- "celltype.genotype"
DE = cell_DE(obj,ct.list,geno2test)
saveRDS(DE,"Analysis/DE/wilcoxon/D3_DE.RDS")


# check proportions of DE vs sample size
prop.sig=get_prop_sig(DE)
obj.sub = obj%>%subset(celltype%in%ct.list)
n.counts = obj.sub[[]]%>%count(celltype,genotype)
for(ct in ct.list){
  print(cor(prop.sig[[ct]][(n.counts%>%filter(celltype ==ct&genotype!="WT"))$genotype],
            (n.counts%>%filter(celltype ==ct&genotype!="WT"))$n,use = 'complete.obs'))
}
# [1] 0.5858891
# [1] 0.2361605

###### D7 ######

obj = D7
summary(as.factor(obj$celltype))
ct.list=c("Posterior foregut_1","Posterior foregut_2")
obj$celltype.genotype <- paste(obj$celltype, obj$genotype, sep = "_")
Idents(obj) <- "celltype.genotype"
DE = cell_DE(obj,ct.list,geno2test)
saveRDS(DE,"Analysis/DE/wilcoxon/D7_DE.RDS")

# check proportions of DE vs sample size
prop.sig=get_prop_sig(DE)
obj.sub = obj%>%subset(celltype%in%ct.list)
n.counts = obj.sub[[]]%>%count(celltype,genotype)
for(ct in ct.list){
  print(cor(prop.sig[[ct]][(n.counts%>%filter(celltype ==ct&genotype!="WT"))$genotype],
            (n.counts%>%filter(celltype ==ct&genotype!="WT"))$n,use = 'complete.obs'))
}
# [1] 0.4367513
# [1] 0.576018



###### D11 ######

obj = D11
summary(as.factor(obj$celltype))
ct.list=c('Pancreatic progenitor_1','Pancreatic progenitor_2')
obj$celltype.genotype <- paste(obj$celltype, obj$genotype, sep = "_")
Idents(obj) <- "celltype.genotype"
DE = cell_DE(obj,ct.list,geno2test)
saveRDS(DE,"Analysis/DE/wilcoxon/D11_DE.RDS")


# check proportions of DE vs sample size
prop.sig=get_prop_sig(DE)
obj.sub = obj%>%subset(celltype%in%ct.list)
n.counts = obj.sub[[]]%>%count(celltype,genotype)
for(ct in ct.list){
  print(cor(prop.sig[[ct]][(n.counts%>%filter(celltype ==ct&genotype!="WT"))$genotype],
            (n.counts%>%filter(celltype ==ct&genotype!="WT"))$n,use = 'complete.obs'))
}
# [1] 0.4077133
# [1] 0.3207598


###### D18 ######

obj = D18
summary(as.factor(obj$celltype))
ct.list=c("SC-enterchromaffin","Nonendocrine","SC-alpha","SC-beta")
obj$celltype.genotype <- paste(obj$celltype, obj$genotype, sep = "_")
Idents(obj) <- "celltype.genotype"
DE = cell_DE(obj,ct.list,geno2test)
saveRDS(DE,"Analysis/DE/wilcoxon/D18_DE.RDS")
DE = readRDS("Analysis/DE/wilcoxon/D18_DE.RDS")

# check proportions of DE vs sample size
prop.sig=get_prop_sig(DE)
obj.sub = obj%>%subset(celltype%in%ct.list)
n.counts = obj.sub[[]]%>%count(celltype,genotype)
for(ct in ct.list){
  print(cor(prop.sig[[ct]][(n.counts%>%filter(celltype ==ct&genotype!="WT"))$genotype],
            (n.counts%>%filter(celltype ==ct&genotype!="WT"))$n,use = 'complete.obs'))
}
# [1] 0.7345467
# [1] 0.7765717
# [1] 0.7021777
# [1] 0.7843352


# p-value histograms
for( day in c("D1","D3","D7","D11","D18")){
  print(day)
  DE = readRDS(paste(wilcoxon.path,day,"_DE.RDS",sep=''))
  pdf(paste("Analysis/DE/wilcoxon/",day,"_DE.pdf"))
  for( i in 1:length(DE)){
    ct = names(DE)[i]
    for(g in names(DE[[i]])){
      if(!all(is.na(DE[[i]][[g]]))){
        df = DE[[i]][[g]]
        print(hist(df$p_val,main = paste(ct,g,sep=' : ')))
      }
    }
  }
  dev.off()
}



##### Pseudobulk test (DESeq2) #####

# use Seurat `FindMarkers`
# results can be reproduced using DESeq2 package


bulk_DE<-function(obj,ct.list,min.pct=0.5){
  # create pseudobulk object
  obj$celltype = as.factor(obj$celltype)
  levels(obj$celltype)=gsub("_","-",levels(obj$celltype)) # rename celltype levels
  
  if(length(levels(obj$celltype))>1){ 
    obj$pseudo.sample.id = as.factor(paste(obj$genotype,obj$sgrna2,obj$celltype,sep='_'))
    pseudo.obj <- AggregateExpression(obj, assays = "RNA", return.seurat = T, 
                                      group.by = c("genotype", "sgrna2", "celltype"))
  }else{ # only one cell type
    obj$pseudo.sample.id = as.factor(paste(obj$genotype,obj$sgrna2,sep='_'))
    pseudo.obj <- AggregateExpression(obj, assays = "RNA", return.seurat = T, group.by = c("genotype", "sgrna2"))
  }
  
  # keep pseudobulk samples that contain > 30 cells
  pseudo.ids.keep=(obj[[]]%>%count(pseudo.sample.id)%>%filter(n>30))$pseudo.sample.id
  pseudo.obj$pseudo.sample.id = rownames(pseudo.obj[[]])
  pseudo.obj=pseudo.obj%>%subset(pseudo.sample.id%in%pseudo.ids.keep)
  # create celltype-geonotype combination id
  if(length(levels(obj$celltype))>1){
    pseudo.obj$celltype.genotype <- paste(pseudo.obj$celltype, pseudo.obj$genotype, sep = "_")
    Idents(pseudo.obj) <- "celltype.genotype"
  }else{
    Idents(pseudo.obj) <- "genotype"
  }
  
  # DE analysis for each celltype-genotype combination
  DE.list = list()
  if(length(levels(obj$celltype))>1){
    for(ct in ct.list){ # >1 cell types in data
      print(ct)
      DE.list[[ct]] = list()
      
      for(g in geno2test){
        if(!paste(ct,g,sep='_')%in%levels(as.factor(pseudo.obj$celltype.genotype))){
          DE.list[[ct]][[g]]=NA
        }else{
          n1 = (pseudo.obj[[]]%>%count(celltype.genotype)%>%filter(celltype.genotype==paste(ct,g,sep='_')))$n
          if(n1<2){
            DE.list[[ct]][[g]]=NA
          }else{
            res = FindMarkers(object = pseudo.obj, 
                              ident.1 =  paste(ct,g,sep='_'),
                              ident.2 = paste(ct,"WT",sep='_'),
                              test.use = "DESeq2",
                              min.cells.group=2,
                              min.pct=min.pct)
            DE.list[[ct]][[g]]<-res%>%dplyr::arrange(p_val)
          }
        }
      }
     }
    }else{ # only 1 cell type in data
      ct = ct.list[1]
      DE.list[[ct]] = list()
      for(g in geno2test){
        if(!g%in%levels(as.factor(pseudo.obj$genotype))){
          DE.list[[ct]][[g]]=NA
        }else{
          n1 = (pseudo.obj[[]]%>%count(genotype)%>%filter(genotype==g))$n
          if(n1<2){
            DE.list[[ct]][[g]]=NA
          }else{
            res <- FindMarkers(object = pseudo.obj , 
                               ident.1 = g,
                               ident.2 = "WT",
                               test.use = "DESeq2",
                               min.cells.group=1,
                               min.pct=0.5)
            DE.list[[ct]][[g]] <-res%>%dplyr::arrange(p_val)
          }
        }
      }
    }
  return(DE.list)
}

geno2test = setdiff(levels(as.factor(D1$genotype)),"WT")

obj = D1
ct.list=c("ESC")
pseudo_DE = bulk_DE(obj,ct.list)
saveRDS(pseudo_DE,"Analysis/DE/pseudo/pseudo_D1_DE.RDS")


obj = D3
ct.list=c("DE-1","DE-2")
pseudo_DE = bulk_DE(obj,ct.list)
saveRDS(pseudo_DE,"Analysis/DE/pseudo/pseudo_D3_DE.RDS")


obj = D7
ct.list=c("Posterior foregut-1","Posterior foregut-2")
pseudo_DE = bulk_DE(obj,ct.list)
saveRDS(pseudo_DE,"Analysis/DE/pseudo/pseudo_D7_DE.RDS")


obj = D11
pseudo_DE = bulk_DE(obj,ct.list)
saveRDS(pseudo_DE,"Analysis/DE/pseudo/pseudo_D11_DE.RDS")



obj = D18
ct.list=c("SC-enterchromaffin","Nonendocrine","SC-alpha","SC-beta")
pseudo_DE = bulk_DE(obj,ct.list)
saveRDS(pseudo_DE,"Analysis/DE/pseudo/pseudo_D18_DE.RDS")




# p-value histograms
for( day in c("D1","D3","D7","D11","D18")){
  print(day)
  pseudo_DE = readRDS(paste(pseudo.path,day,"_DE.RDS",sep=''))
  pdf(paste("Analysis/DE/pseudo/",day,".p.val.pdf"))
  for( i in 1:length(pseudo_DE)){
    ct = names(pseudo_DE)[i]
    for(g in names(pseudo_DE[[i]])){
      if(!all(is.na(pseudo_DE[[i]][[g]]))){
        df = pseudo_DE[[i]][[g]]
        print(hist(df$p_val,main = paste(ct,g,sep=' : ')))
      }
    }
  }
  dev.off()
}

# p-values for genes also tested in cell-level Wilcoxon test

for( day in c("D1","D3","D7","D11","D18")){
  print(day)
  pseudo_DE = readRDS(paste(pseudo.path,day,"_DE.RDS",sep=''))
  DE = readRDS(paste(wilcoxon.path,day,"_DE.RDS",sep=''))
  names(pseudo_DE)=names(DE)
  pdf(paste("Analysis/DE/pseudo/filteredgenes/",day,".p.val.pdf"))
  for(ct in names(pseudo_DE)){
    print(ct)
    for(g in names(pseudo_DE[[ct]])){
      print(g)
      if(!all(is.na(pseudo_DE[[ct]][[g]]))&!all(is.na(DE[[ct]][[g]]))){
        pseudo_res = pseudo_DE[[ct]][[g]]
        res = DE[[ct]][[g]]
        print(dim(res))
        print(dim(pseudo_res))
        pseudo_res$genes = rownames(pseudo_res);res$genes = rownames(res)
        genes = intersect(rownames(pseudo_res),rownames(res))
        pseudo_res=pseudo_res[genes,]
        print(hist(pseudo_res$p_val,main = paste(ct,g,sep=' : ')))
      }
    }
  }
  dev.off()
}


##### compare p-values from pseudobulk and cell-level DE analysis #####

pseudo.path = "Analysis/DE/pseudo/pseudo_"
wilcoxon.path = "Analysis/DE/wilcoxon/"


# plot -log2(p.val) from Wilconxon and pseudobulk DE analysis

for( day in c("D1","D3","D7","D11","D18")){
  print(day)
  pseudo_DE = readRDS(paste(pseudo.path,day,"_DE.RDS",sep=''))
  DE = readRDS(paste(wilcoxon.path,day,"_DE.RDS",sep=''))
  names(pseudo_DE)=names(DE)
  pdf(paste("Analysis/DE/comparison/",day,".p.val.pdf"))
  for(ct in names(pseudo_DE)){
    print(ct)
    for(g in names(pseudo_DE[[ct]])){
      print(g)
      if(!all(is.na(pseudo_DE[[ct]][[g]]))&!all(is.na(DE[[ct]][[g]]))){
        pseudo_res = pseudo_DE[[ct]][[g]]
        res = DE[[ct]][[g]]
        print(dim(res))
        print(dim(pseudo_res))
        pseudo_res$genes = rownames(pseudo_res);res$genes = rownames(res)
        genes = intersect(rownames(pseudo_res),rownames(res))
        pseudo_res=pseudo_res[genes,];res = res[genes,]
        c = cor(-log(res$p_val,2),-log(pseudo_res$p_val,2),method = "spearman",use = "complete.obs")
        text = paste("cor=",round(c,3),sep='')
        df = data.frame(wilcoxon = -log(res$p_val,2),pseudobulk = -log(pseudo_res$p_val,2),gene =res$genes )
        fig=ggplot(data = df, aes(x = wilcoxon, y = pseudobulk)) +
          ggpointdensity::geom_pointdensity()+
          theme_minimal()+
          ggtitle( paste(ct,g,text,sep = ";"))+
          guides(fill = guide_legend(title ='number of neighbors'))+
          ggrepel::geom_text_repel(data = subset(df, pseudobulk > quantile(df$pseudobulk,0.999,na.rm=T) | wilcoxon > quantile(df$wilcoxon,0.999,na.rm=T)), 
                                   aes(label = gene))
        print(ggrastr::rasterize(fig))
      }
    }
  }
  dev.off()
}



###### record correlations ######

cor.mat.list = list()
for( day in c("D1","D3","D7","D11","D18")){
  pseudo_DE = readRDS(paste(pseudo.path,day,"_DE.RDS",sep=''))
  DE = readRDS(paste(wilcoxon.path,day,"_DE.RDS",sep=''))
  names(pseudo_DE)=names(DE)
  cor.mat = matrix(NA,nrow = length(pseudo_DE),ncol = length(pseudo_DE[[1]]))
  rownames(cor.mat) = names(pseudo_DE); colnames(cor.mat) = names(pseudo_DE[[1]])
  for(i in 1:length(pseudo_DE)){
    ct = names(pseudo_DE)[i]
    for(j in 1:length(pseudo_DE[[i]])){
      g = names(pseudo_DE[[ct]])[j]
      if(!all(is.na(pseudo_DE[[ct]][[g]]))&!all(is.na(DE[[ct]][[g]]))){
        pseudo_res = pseudo_DE[[ct]][[g]]
        res = DE[[ct]][[g]]
        pseudo_res$genes = rownames(pseudo_res);res$genes = rownames(res)
        genes = intersect(rownames(pseudo_res),rownames(res))
        pseudo_res=pseudo_res[genes,];res = res[genes,]
        c = cor(-log(res$p_val,2),-log(pseudo_res$p_val,2),method = "spearman",use = "complete.obs")
        cor.mat[i,j] = c
      }
    }
  }
  cor.mat.list[[day]] = cor.mat
}
# spread to data frame
df.lst = lapply(c("D1","D3","D7","D11","D18"),function(day){
  df = reshape::melt(cor.mat.list[[day]])
  df$day = day
  return(df)
})
df = Reduce(rbind,df.lst)
colnames(df)[1:3] = c("celltype","genotype","cor")
df$day = factor(df$day,levels=c("D1","D3","D7","D11","D18"))

# boxplots by cell types
pdf("Analysis/DE/comparison/cor.pdf",width = 5, height = 3)
df%>%ggplot(aes(x=celltype,y=cor,color = day))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=0.8,position = position_jitter(width = 0.01, height = 0.01))+
  theme_minimal()+
  ylab("Spearman Correlation")+
  xlab("")+
  guides(color = guide_legend(title =''))+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8))
dev.off()







