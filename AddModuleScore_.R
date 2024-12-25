setwd('/media/user/sdh/cmq_wj/AddModuleScore')
.libPaths("/media/user/sdf/R.lib/library/wj/reg")
library(Seurat)
library(ggplot2)
colors <- c('#7ab800',"#6e2585","#b7295a","#003a70",
            "#5482ab","#dc5034","#0085c3","#f2af00","#ce1126")
lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds') 
path<-dir('/media/user/sdh/cmq_wj/AddModuleScore')
for (i1 in path){
he<-read.csv(i1)

gene <- list(he$V1)
lung <- AddModuleScore(  #对每个模块的基因集打分
  object = lung,
  features = gene,
  ctrl = 100, #默认值是100
  name = strsplit(i1,'_GOGENE.CSV')[[1]][1])
}

VlnPlot(lung, features=paste0(strsplit(i1,'_GOGENE.CSV')[[1]][1],'1'))

for (i in path){
 
VlnPlot(lung, features=paste0(strsplit(i,'_GOGENE.CSV')[[1]][1],'1'),
                                     cols=colors,pt.size = 0) +
    geom_boxplot(width=.2,col="black",fill="white")+  
    NoLegend()
 ggsave(paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',
               strsplit(i,'_GOGENE.CSV')[[1]][1],'.pdf'),width = 4.5,height = 4)
 
 # FeaturePlot(lung,features=paste0(strsplit(i,'_GOGENE.CSV')[[1]][1],'1'),label = T)
 # 
 # ggsave(paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',strsplit(i,'_GOGENE.CSV')[[1]][1],'_tsne.pdf'),
 #        width = 6,height = 5.2)
 # 
}


###########  G5 G6 G7 G8

for(i in levels(lung@meta.data[["orig.ident"]])){

  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@meta.data[["orig.ident"]])])
  lung_se<- subset(lung,cells=clust.name)
 Idents(lung_se)<-'active.ident'

 path<-dir('/media/user/sdh/cmq_wj/AddModuleScore')
 for (i1 in path){
   he<-read.csv(i1)
  
  gene <- list(he$V1)
   lung_se <- AddModuleScore(  #对每个模块的基因集打分
     object = lung_se,
     features = gene,
       ctrl = 100, #默认值是100
    name = strsplit(i1,'_GOGENE.CSV')[[1]][1])
  
  for (uu in c('Lrp1','Itgb5','Cd300lf','Calr','Cd14','Ager')){
    pdf(paste0('/media/user/sdh/cmq_wj/AM/',i,'_',uu,'.vln.pdf'),
        width = 4.5,height = 4)
    AA<-VlnPlot(lung_se, features=uu,cols=colors,pt.size = 0)+
      geom_boxplot(width=.1,col="black",fill="white")+  
      NoLegend()
    print(AA)
    dev.off()
  }
  }
  

  for (ii in path){
    VlnPlot(lung_se, features=paste0(strsplit(ii,'_GOGENE.CSV')[[1]][1],'1'),
            cols=colors,pt.size = 0)+
      geom_boxplot(width=.2,col="black",fill="white")+  
      NoLegend()
    ggsave(paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',
                  i,'_',strsplit(ii,'_GOGENE.CSV')[[1]][1],'-1.pdf'),
           width = 4.5,height = 4)
    
    # FeaturePlot(lung_se,features=paste0(strsplit(ii,'_GOGENE.CSV')[[1]][1],'1'),label = T,
    #             cols=c('#c7d9ff','grey','#ff4242')
    #             )
    # 
    # ggsave(paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',i,'_',strsplit(ii,'_GOGENE.CSV')[[1]][1],'_tsne.pdf'),
    #        width = 6,height = 5.2)
    
    
    
  }
}


for (uu in c('Lrp1','Itgb5','Cd300lf','Calr','Cd14','Ager')){
  pdf(paste0('/media/user/sdh/cmq_wj/AM/',uu,'.vln.pdf'),
      width = 4.5,height = 4)
  AA<-VlnPlot(lung, features=uu,cols=colors,pt.size = 0，y.max =0.4)+
  geom_boxplot(width=.1,col="black",fill="white")+  
  NoLegend()
  print(AA)
  dev.off()
}

ggsave('/media/user/sdh/cmq_wj/AM/AXL.VLN.PDF',
       width = 4.5,height = 4)




#########3   M1  M2
for(i in levels(lung@meta.data[["orig.ident"]])){
  
  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@meta.data[["orig.ident"]])])
  lung_se<- subset(lung,cells=clust.name)
  
  lung_se <- AddModuleScore(  #对每个模块的基因集打分
  object = lung_se,
  features =list(c('Cxcl9',
                   'Cxcl10',
                   'Cxcl11',
                   'Nos2',
                   'Tnf',
                   'Cd86'))  ,
  ctrl = 100, #默认值是100
  name = 'M1_marker')

  lung_se <- AddModuleScore(  #对每个模块的基因集打分
  object = lung_se,
  features =list(c('Chil3',
                   'Arg1',
                   'Retnla',
                   'Mrc1',
                   'Tgm2'))  ,
  ctrl = 100, #默认值是100
  name = 'M2_marker')


# colors <- c('#7ab800',"#6e2585","#b7295a","#003a70",
#             "#5482ab","#dc5034","#0085c3","#f2af00","#ce1126")

for( ii in c('M1_marker1','M2_marker1')){
  VlnPlot(lung_se, features=ii,
          cols=colors,pt.size = 0,y.max =3,idents=c('Monocyte_CCL3','Monocyte_Ly6c_low','Macrophage_Ki67','Monocyte_Ly6c_high','Macrophage_alveolar') )+
    geom_boxplot(width=.1,col="black",fill="white")+  
    NoLegend()
  ggsave(paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',i,'_',ii,'-2.pdf'),width = 3,height = 4)
}
}




for(i in levels(lung@meta.data[["orig.ident"]])){
  
  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@meta.data[["orig.ident"]])])
  lung_se<- subset(lung,cells=clust.name)
  # Idents(lung_se)<-'active.ident'
  
  path<-dir('/media/user/sdh/cmq_wj/AddModuleScore')
  for (i1 in path){
    he<-read.csv(i1)
    
    gene <- list(he$V1)
    lung_se <- AddModuleScore(  #对每个模块的基因集打分
      object = lung_se,
      features = gene,
      ctrl = 100, #默认值是100
      name = strsplit(i1,'_GOGENE.CSV')[[1]][1]) 
    UU<-VlnPlot(lung_se, features=paste0(strsplit(i1 ,'_GOGENE.CSV')[[1]][1],'1'),
            cols=colors,pt.size = 0,y.max=0.4)+
      geom_boxplot(width=.2,col="black",fill="white")+  
      NoLegend()
    ggsave(UU,paste0('/media/user/sdh/cmq_wj/AddModuleScore_polt/',
                  i,'_',strsplit(i1 ,'_GOGENE.CSV')[[1]][1],'-1.pdf'),
           width = 4.5,height = 4)
  }

  
  }
  