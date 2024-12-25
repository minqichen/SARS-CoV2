setwd('/media/user/sdh/cmq_wj')
.libPaths("/media/user/sdf/R.lib/library/wj/reg")
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(monocle)
library(devtools)
library(reticulate)
library(data.table)
library(pheatmap)
library(Rcpp)
library(harmony)
use_python("/home/mayuanchen/anaconda3/bin/python",required = T)
py_config()
py_module_available("umap")
library(PythonInR)

lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds') 


for(i in levels(lung@active.ident)){
  
  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@active.ident)])
  lung_se<- subset(lung,cells=clust.name)
  
  # saveRDS(lung_se,paste0('/media/user/sdh/cmq_wj/G_rds/',i,'.rds'))
  
  # all.markers_Nfc<-FindAllMarkers(lung_se,only.pos = F,min.pct = 0.01,logfc.threshold = 0)
  # saveRDS(all.markers_Nfc,paste0('/media/user/sdh/cmq_wj/G_allmarker/',i,'_all.markers.rds'))
  Idents(lung_se)<-'orig.ident'
  G65<-FindMarkers(
    lung_se,
    ident.1 = 'G6',
    ident.2 = 'G5',
    only.pos = F,min.pct = 0.1,logfc.threshold = 0)
  write.csv(G65,paste0('/media/user/sdh/cmq_wj/G_allmarker/',i,'_G6_G5.CSV'))
  
  G87<-FindMarkers(
    lung_se,
    ident.1 = 'G8',
    ident.2 = 'G7',
    only.pos = F,min.pct = 0.1,logfc.threshold = 0)
  write.csv(G87,paste0('/media/user/sdh/cmq_wj/G_allmarker/',i,'_G8_G7.CSV'))
  
  G75<-FindMarkers(
    lung_se,
    ident.1 = 'G7',
    ident.2 = 'G5',
    only.pos = F,min.pct = 0.1,logfc.threshold = 0)
  write.csv(G75,paste0('/media/user/sdh/cmq_wj/G_allmarker/',i,'_G7_G5.CSV'))
  
  
  G86<-FindMarkers(
    lung_se,
    ident.1 = 'G8',
    ident.2 = 'G6',
    only.pos = F,min.pct = 0.1,logfc.threshold = 0)
  write.csv(G75,paste0('/media/user/sdh/cmq_wj/G_allmarker/',i,'_G8_G6.CSV'))
}

##   volcano
setwd('/media/user/sdh/cmq_wj/G_allmarker')
library(org.Mm.eg.db)
library(clusterProfiler)
path<-getwd()
dir(path)
aa<-dir(path)[grep(".CSV", dir(path))]
aa<-unlist( strsplit(aa, ".CSV", fixed= T))

for (i in aa) {
  gene<-read.csv(paste0(i,'.CSV'))
  gene$thres<-  as.factor(ifelse(gene[,6] < 0.05 & abs(gene[,3]) > 0.58, 
                                 ifelse(gene[,3]> 0.58,'Up','Down'),'NoSignifi'))
  
  ggplot(data = gene, aes(x=avg_logFC, y=-(log10 (gene[,6])),colour= thres )) +
    geom_point() +
    scale_color_manual(values=c("blue", "grey","red"))+
    labs(x="avg_logFC",y="-log10 P adj",title=i)+ 
    # annotate("text", x = 10, y = 5, label =length(gene$thres[gene$thres=='Up']) ,
    #          colour= "red", hjust = -0.08, size = 5 )+
    # annotate("text", x = 10, y = -5, label =length(gene$thres[gene$thres=='Down']) ,
    #          colour= "blue", hjust = -0.08, size = 5 )+
    theme_bw()+#去除背景
    theme(panel.grid.minor = element_blank(),panel.grid.major =element_blank())#去除网格线
  ggsave(paste0('../volcano/',i,'_volcano.pdf'),width = 6,height = 5)
}
