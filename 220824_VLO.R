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
FeaturePlot(lung,'Csf2',split.by = 'orig.ident',reduction = 'harmony_umap')
DimPlot(lung,label = T)
VlnPlot(lung,'Csf2')


G5<-readRDS('/media/user/sdh/cmq_wj/G_rds/G5.rds')
G6<-readRDS('/media/user/sdh/cmq_wj/G_rds/G6.rds')
G7<-readRDS('/media/user/sdh/cmq_wj/G_rds/G7.rds')
G8<-readRDS('/media/user/sdh/cmq_wj/G_rds/G8.rds')




lung<-AddMetaData(lung,c(paste0(G5@active.ident,'_G5'),
                        paste0(G6@active.ident,'_G6') ,
                         paste0(G7@active.ident,'_G7'),
                         paste0(G8@active.ident,'_G8') ), col.name = 'vlnplot_group') #添加分组为metadata
Idents(lung)<-'vlnplot_group'



