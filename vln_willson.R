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
library(Rcpp)
library(harmony)
use_python("/home/mayuanchen/anaconda3/bin/python",required = T)
py_config()
py_module_available("umap")
library(PythonInR)

source('/media/user/sdh/cmq_wj/others/vlnplot/violin_gene_exp_function_0.2_jitter.R')

AM<-readRDS("/media/user/sdh/cmq_wj/G_rds/Macrophage_alveolar.rds")
AM<-FindVariableFeatures(AM, selection.method = "vst", nfeatures = 32285)
AM <- ScaleData(AM, verbose = FALSE)
Idents(AM)<-AM@meta.data[["orig.ident"]]


name<-as.data.frame( AM@assays$RNA@counts)
kk_cluster <-as.matrix( AM@meta.data$orig.ident)# kk_cluster是细胞分组信息
rownames(kk_cluster)<-as.character(colnames(name)) 
colnames(kk_cluster)<-'clusters'
mat<-GetAssayData(AM,slot='scale.data')

mat<-as.data.frame( AM@assays$RNA@counts)


mat_1<-as.matrix(mat[,names(cluster_info)])
gene_sel<-c('Havcr2','Calr','Lrp1', 'Itgb5', 'Cd300lf', 'Cd14', 'Cd36','Mertk','Alx')
mat_2<-mat[gene_sel,]
mat_2<-subset(mat,rownames(mat) %in% gene_sel)



kk<-as.data.frame(mat_2)  # kk是表达矩阵 
pdf('VlnPlot_1213_wilcox.pdf',width = 6,height = 4)
for( i in c(1:6)){

  kk1<-as.matrix(kk[i,]) 

violin_gene_exp(gene_sel[i],kk1,kk_cluster)

}
dev.off()
# violin_gene_exp<-function(gene, rpkm, conditions=conditions,  test=TRUE)

pdf('VlnPlot_1214_wilcox.pdf',width = 6,height = 4)
gene_sel<-c('Havcr2','Calr','Lrp1', 'Itgb5', 'Cd300lf', 'Cd14', 'Cd36','Mertk','Alx')
for( i in c(1:9)){
gene<-gene_sel[i]
exp<-as.numeric(kk[gene,])
gene_exp<-data.frame(
  cell=colnames(kk),
  clusters=kk_cluster,
  gene=exp
)

test<-compare_means(gene~clusters, data=gene_exp, method="wilcox.test")
#print(test)
#print(compare_means(gene ~ clusters, data = gene_exp, , method = "wilcox.test", ref.group = ".all."))
my_comparisons=list(c("G7", "G8"),c('G5','G6'))
p<-ggboxplot(gene_exp, x="clusters", y="gene", color="white")+
  # geom_violin(scale="width", width=0.7, adjust=.5,aes(fill=clusters)) +
  stat_summary(fun.y=mean, geom="point", shape=21, size=3,  stroke=1, fill="white")+
  geom_boxplot(fill=rep(c("#add8e6","#ee5e49"),2), outlier.shape = NA, width = 0.5)+ #箱子的颜色
  
  geom_jitter(size=0.3)+  	  
  # geom_hline(yintercept=mean(gene_exp$gene), linetype=2)+
  # geom_hline(yintercept=mean(gene_exp[494:542,]$gene), linetype=2)+  #虚线设为 T-HSC的中值
  
  # scale_fill_aaas()+
  # scale_fill_manual(values = c('#bb0021','#008280','#008b45','#631879','#5f559b','#a20056','#3b4992','#ee0000'))+
  scale_fill_manual(values =rep(c("#add8e6","#ee5e49"),2))+
  theme_bw()+
  ggtitle(gene)+
  expand_limits(y=c(0, max(gene_exp$gene)*1.1))+
  # scale_x_discrete(labels = c('T-HSC','LT-HSC','CLP','CMP','Tcm(LN)','Tcm(SP)','CD4T','CD8T')) +
  stat_compare_means(comparisons=my_comparisons,method = "wilcox.test", label.y = max(gene_exp$gene)+0.5)+      # Add global p-value
  
  # stat_compare_means(comparisons=my_comparisons)+
  # stat_compare_means(
  #   label="p.signif",  #p.format
  #   method="kruskal.test",
  #   ref.group="G7",    #  .all.
  #   label.y=max(gene_exp$gene)*1.05,
  #   size=6
  # )+#Pairwise comparison against all
  theme(
    plot.title=element_text(size=18, face="bold.italic", hjust=0.5),
    axis.text=element_text(size=6),
    #axis.title=element_text(size=16),
    axis.title=element_blank(),
    legend.text=element_text(size=16),
    legend.title=element_blank(),
    aspect.ratio=0.5,
    legend.position="none",
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  )
print(p)
}
dev.off()
