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
library(ComplexHeatmap)

lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds') 
# colors <- c("#003a70","#0085c3","#7ab800","#f2af00","#dc5034",
#             "#ce1126","#b7295a","#5482ab","#6e2585","#71c6c1",
#             "#862633","#009bbb","#444444","#8b8c8d")

colors <- c("#003a70","#0085c3","#7ab800","#f2af00","#dc5034",
            "#ce1126","#b7295a","#5482ab","#6e2585")


lung <- ScaleData(lung, verbose = FALSE)

DoHeatmap(lung,features = top50$gene)
ggsave('heatmap_top50markergen.pdf',width=10,height=9)



AM<-readRDS("/media/user/sdh/cmq_wj/G_rds/Macrophage_alveolar.rds")
AM<-FindVariableFeatures(AM, selection.method = "vst", nfeatures = 32285)
AM <- ScaleData(AM, verbose = FALSE)
Idents(AM)<-AM@meta.data[["orig.ident"]]

NK<-readRDS("/media/user/sdh/cmq_wj/G_rds/NK_cell.rds")
NK<-FindVariableFeatures(NK, selection.method = "vst", nfeatures = 32285)
NK <- ScaleData(NK, verbose = FALSE)
Idents(NK)<-NK@meta.data[["orig.ident"]]

T_cell<-readRDS("/media/user/sdh/cmq_wj/G_rds/T_cell.rds")
T_cell<-FindVariableFeatures(T_cell, selection.method = "vst", nfeatures = 32285)
T_cell <- ScaleData(T_cell, verbose = FALSE)
Idents(T_cell)<-T_cell@meta.data[["orig.ident"]]

B_cell<-readRDS("/media/user/sdh/cmq_wj/G_rds/B_cell.rds")
B_cell<-FindVariableFeatures(B_cell, selection.method = "vst", nfeatures = 32285)
B_cell <- ScaleData(B_cell, verbose = FALSE)
Idents(B_cell)<-B_cell@meta.data[["orig.ident"]]

Monocyte_CCL3<-readRDS("/media/user/sdh/cmq_wj/G_rds/Monocyte_CCL3.rds")
Monocyte_CCL3<-FindVariableFeatures(Monocyte_CCL3, selection.method = "vst", nfeatures = 32285)
Monocyte_CCL3 <- ScaleData(Monocyte_CCL3, verbose = FALSE)
Idents(Monocyte_CCL3)<-Monocyte_CCL3@meta.data[["orig.ident"]]

Monocyte_Ly6c_low<-readRDS("/media/user/sdh/cmq_wj/G_rds/Monocyte_Ly6c_low.rds")
Monocyte_Ly6c_low<-FindVariableFeatures(Monocyte_Ly6c_low, selection.method = "vst", nfeatures = 32285)
Monocyte_Ly6c_low <- ScaleData(Monocyte_Ly6c_low, verbose = FALSE)
Idents(Monocyte_Ly6c_low)<-Monocyte_Ly6c_low@meta.data[["orig.ident"]]

Macrophage_Ki67<-readRDS("/media/user/sdh/cmq_wj/G_rds/Macrophage_Ki67.rds")
Macrophage_Ki67<-FindVariableFeatures(Macrophage_Ki67, selection.method = "vst", nfeatures = 32285)
Macrophage_Ki67 <- ScaleData(Macrophage_Ki67, verbose = FALSE)
Idents(Macrophage_Ki67)<-Macrophage_Ki67@meta.data[["orig.ident"]]

Monocyte_Ly6c_high<-readRDS("/media/user/sdh/cmq_wj/G_rds/Monocyte_Ly6c_high.rds")
Monocyte_Ly6c_high<-FindVariableFeatures(Monocyte_Ly6c_high, selection.method = "vst", nfeatures = 32285)
Monocyte_Ly6c_high <- ScaleData(Monocyte_Ly6c_high, verbose = FALSE)
Idents(Monocyte_Ly6c_high)<-Monocyte_Ly6c_high@meta.data[["orig.ident"]]

neu_seurat<-readRDS('/media/user/sdh/cmq_wj/G_rds/Neutrophil.rds')
neu_seurat<-FindVariableFeatures(neu_seurat, selection.method = "vst", nfeatures = 32285)
neu_seurat <- ScaleData(neu_seurat, verbose = FALSE)
Idents(neu_seurat)<-neu_seurat@meta.data[["orig.ident"]]

seurat_data<-list(AM,NK,T_cell,B_cell,Monocyte_CCL3,Monocyte_Ly6c_low,
                  Macrophage_Ki67,Monocyte_Ly6c_high,neu_seurat)
seurat_name<-c('AM','NK','T_cell','B_cell','Monocyte_CCL3','Monocyte_Ly6c_low',
                  'Macrophage_Ki67','Monocyte_Ly6c_high','Neutrophil')

AM_efferocytosis_receptors <-c('Adgrb1',
      'Havcr1',
      'Havcr2',
      'Timd4',
      'Stab1',
      'Stab2',
      'Ager',
      'Lrp1',
      'Mertk',
      'Axl',
      'Tyro3',
      'Itgb3',
      'Itgb5', 
      'Cd36',
      'Cd300lf',
      'Icam5',
      'Calr',
      'Cd14',
      'Pecam1')
Idents(AM)<-AM@meta.data[["orig.ident"]]
DoHeatmap(AM,features = AM_efferocytosis_receptors)
ggsave('/media/user/sdh/cmq_wj/heatmap_1212/AM_efferocytosis_receptors_heatmap.pdf',
       width=8,height=6)
VlnPlot(AM,features = AM_efferocytosis_receptors)
ggsave('/media/user/sdh/cmq_wj/heatmap_1212/AM_efferocytosis_receptors.pdf',width=10,height=12)

AM_donot_eatme<-c('Sirpa',
 'Pecam1',
'Siglecg',
'Lilrb1')
DoHeatmap(AM,features = AM_donot_eatme)
VlnPlot(AM,features = AM_donot_eatme)
ggsave('/media/user/sdh/cmq_wj/heatmap_1212/AM_donot_eatme.pdf',width=5,height=2.5)
 


all_1<-c('Bnip3','Gsdme','Gzma','Gzmb','Lamp1','Nkg7','Srgn','Ube4b')
for (i in c(1:9)){
DoHeatmap(seurat_data[[i]],features = all_1)
   ggsave(paste0('/media/user/sdh/cmq_wj/heatmap_1212/',seurat_name[i],
                 '_granzyme_mediated_programmed_cell_death.pdf'),
          width = 7,height = 8)
}





cell_death<-c('Bnip3',
              'Gsdme',
              'Gzma',
              'Gzmb',
              'Lamp1',
              'Nkg7',
              'Srgn',
              'Ube4b')
DoHeatmap(AM,features = cell_death)
VlnPlot(AM,features = cell_death)

################33
neu_seurat<-readRDS('/media/user/sdh/cmq_wj/G_rds/Neutrophil.rds')
neu_seurat<-FindVariableFeatures(neu_seurat, selection.method = "vst", nfeatures = 32285)
neu_seurat <- ScaleData(neu_seurat, verbose = FALSE)
mat<-GetAssayData(neu_seurat,slot='scale.data')
Idents(neu_seurat)<-neu_seurat@meta.data[["orig.ident"]]
cluster_info<-sort(neu_seurat@meta.data[["orig.ident"]])
neu_1<-c('Cd47','Pecam1','Cd24a')
neu_2<-c('Panx1',
         'Mbtps1',
         'Cx3cl1',
         'Gpr132',
         'Gzma')
mat_1<-as.matrix(mat[,names(cluster_info)])
# subset(neu_1,neu_1 %in% rownames(mat_1))
mat_2<-mat[subset(neu_2,neu_2 %in% rownames(mat_1)),]
Heatmap(mat_2,
        cluster_columns = F,
        show_column_names = F,
        column_split=cluster_info)

VlnPlot(neu_seurat,neu_1)
ggsave('/media/user/sdh/cmq_wj/heatmap_1212/neu_dont_eat_me.pdf',width = 10,height = 2.5)
VlnPlot(neu_seurat,neu_2)
ggsave('/media/user/sdh/cmq_wj/heatmap_1212/neu_find_me.pdf',width = 10,height = 5)
DoHeatmap(neu_seurat,features = neu_1)
DoHeatmap(neu_seurat,features = neu_2)
dir('/media/user/sdh/cmq_wj/heatmap_gene')

pdf('/media/user/sdh/cmq_wj/heatmap_1212/neu_celldeath.pdf',width=10,height = 10)
for (i in dir('/media/user/sdh/cmq_wj/heatmap_gene')){
   neu<-read.table(paste0('/media/user/sdh/cmq_wj/heatmap_gene/',i))
   neu<-as.character(neu$V1)
   oo<-DoHeatmap(neu_seurat,features = neu)
   print(oo)
}
dev.off()


pdf('/media/user/sdh/cmq_wj/heatmap_1212/neu_celldeath_1.pdf',width=10,height = 10)
for (i in dir('/media/user/sdh/cmq_wj/heatmap_gene')){
   neu<-read.table(paste0('/media/user/sdh/cmq_wj/heatmap_gene/',i))
   neu<-as.character(neu$V1)
   mat_1<-as.matrix(mat[,names(cluster_info)])
   # subset(neu_1,neu_1 %in% rownames(mat_1))
   mat_2<-mat[subset(neu,neu %in% rownames(mat_1)),]
   oo<-Heatmap(mat_2,
           cluster_columns = F,
           show_column_names = F,
           column_split=cluster_info)
   print(oo)
}
dev.off()




pdf('/media/user/sdh/cmq_wj/heatmap_1212/neu_celldeath_2.pdf',width=10,height = 30)
for (i in dir('/media/user/sdh/cmq_wj/heatmap_gene')){
   neu<-read.table(paste0('/media/user/sdh/cmq_wj/heatmap_gene/',i))
   pdf(paste0('/media/user/sdh/cmq_wj/heatmap_1212/',i,'neu_celldeath_2.pdf'),width=10,height = length(neu$V1)/10+3 )
   neu<-as.character(neu$V1)
   oo<-DotPlot(neu_seurat, features = neu)+coord_flip()+
      theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
      labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
      scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
   print(oo)
   dev.off()
}





DotPlot(scedata, features = marker)+coord_flip()+
   theme_bw()+
   theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
   labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
   scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

########### 20230116
neu_seurat<-readRDS('/media/user/sdh/cmq_wj/G_rds/Neutrophil.rds')
neu_seurat<-FindVariableFeatures(neu_seurat, selection.method = "vst", nfeatures = 32285)
neu_seurat <- ScaleData(neu_seurat, verbose = FALSE)
mat<-GetAssayData(neu_seurat,slot='scale.data')
Idents(neu_seurat)<-neu_seurat@meta.data[["orig.ident"]]
cluster_info<-sort(neu_seurat@meta.data[["orig.ident"]])

mat_1<-as.matrix(mat[,names(cluster_info)])


pdf('/media/user/sdh/cmq_wj/heatmap_1212/neu_230116.pdf',width=12,height = 12)
for (i in dir('/media/user/sdh/cmq_wj/neuu')){
   neu<-read.table(paste0('/media/user/sdh/cmq_wj/neuu/',i))
   neu<-as.character(neu$V1)
   mat_1<-as.matrix(mat[,names(cluster_info)])
   # subset(neu_1,neu_1 %in% rownames(mat_1))
   mat_2<-mat[subset(neu,neu %in% rownames(mat_1)),]
   oo<-Heatmap(mat_2,
               cluster_columns = F,
               show_column_names = F,
               column_split=cluster_info,
               name=gsub('.txt','',i),
               heatmap_width  = unit(20, "cm"),
               heatmap_height  = unit(20, "cm")
               )
   print(oo)
}
dev.off()

pdf('/media/user/sdh/cmq_wj/heatmap_1212/neu_230117_Dotplot.pdf',width=15,height = 12)
for (i in dir('/media/user/sdh/cmq_wj/neuu')){
   neu<-read.table(paste0('/media/user/sdh/cmq_wj/neuu/',i))
   neu<-as.character(neu$V1)
   aa<-DotPlot(neu_seurat,features =  neu) 
   print(aa)
}
dev.off()
