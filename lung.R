
.libPaths('/media/user/sdf/R.lib/library/cmq')
library(ggsci)
library(ggpubr)
library(farver)

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
.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Rcpp)
library(harmony)
use_python("/home/mayuanchen/anaconda3/bin/python",required = T)
py_config()
py_module_available("umap")
library(PythonInR)

lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds') 
# colors <- c("#003a70","#0085c3","#7ab800","#f2af00","#dc5034",
#             "#ce1126","#b7295a","#5482ab","#6e2585","#71c6c1",
#             "#862633","#009bbb","#444444","#8b8c8d")

colors <- c("#003a70","#0085c3","#7ab800","#f2af00","#dc5034",
            "#ce1126","#b7295a","#5482ab","#6e2585")
# DimPlot(lung,split.by = 'orig.ident',cols = colors)
# ggsave('UMAP_group.pdf',width = 20,height = 5)

# all.markers<-FindAllMarkers(lung,only.pos = T,min.pct = 0.1,logfc.threshold = 0.2)
# saveRDS(all.markers,'all.markers.rds')
all.markers<-readRDS('all.markers.rds')

lung <- ScaleData(lung, verbose = FALSE)
top50 <- all.markers%>%group_by(cluster)%>%top_n(n=50,wt=avg_logFC) 

DoHeatmap(lung,features = top50$gene)
ggsave('heatmap_top50markergen.pdf',width=10,height=9)

top5 <- all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_logFC) 
DotPlot(lung, features = top5$gene)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 90,hjust = 1,vjust=1))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
ggsave('bubble_top5markergene.pdf',width = 8,height = 8)

DimPlot(lung,cols = colors,label = T)
ggsave('UMAP.pdf',width = 7.5,height = 5)

lung.name <- c(rownames(lung@meta.data)
                    [which(lung@meta.data[["COV2"]] >=1) ])
DimPlot(lung,reduction ='harmony_umap',
        cells.highlight = lung.name,
        split.by = 'orig.ident',cols=c('grey'),
        cols.highlight =c("red"))
ggsave('/media/user/sdh/cmq_wj/COV2_harmony_umap.pdf',width = 13,height = 3)

all.markers_Nfc<-FindAllMarkers(lung,only.pos = F,min.pct = 0.1,logfc.threshold = 0)
# saveRDS(all.markers,'all.markers.rds')

FeaturePlot(lung,'Axl',split.by = 'orig.ident',reduction ='harmony_umap')
ggsave('/media/user/sdh/cmq_wj/Axl_harmony_umap.pdf',width = 13,height = 3)

FeaturePlot(lung,features='COV2',split.by = 'orig.ident',
            reduction ='harmony_umap',cols = c("grey", "red"))

VlnPlot(lung,features='COV2',group.by  = 'orig.ident')
