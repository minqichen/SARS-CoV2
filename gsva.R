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
.libPaths("/media/user/sdf/R.lib/library/cmq")
library(msigdbr)
library(GSVA)
library(pheatmap)
lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds') 

# Idents(lung) <- "cell_type" 
expr <- AverageExpression(lung, assays = "RNA", slot = "data")[[1]] 
expr <- expr[rowSums(expr)>0,] #选取非零基因 
expr <- as.matrix(expr) 
head(expr)

genesets <- msigdbr(species = "Mus musculus", category= "C5")


genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame() 
genesets <- split(genesets$gene_symbol, genesets$gs_name)
# genesets1<-subset(genesets,grep('CELL_CYCLE',genesets))
# write.csv(genesets,'c5_GOBP_GENE2GSVA.csv')



gsva.res <- gsva(expr, genesets, method="ssgsea", parallel.sz=8) 
saveRDS(gsva.res, "C5.gsva.res.rds") 
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "C5.gsva_res.csv", row.names = F)





for(i in levels(lung@active.ident)){
  
  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@active.ident)])
  lung_se<- subset(lung,cells=clust.name)
  Idents(lung_se)<-'orig.ident'
  expr <- AverageExpression(lung_se, assays = "RNA", slot = "data")[[1]] 
  expr <- expr[rowSums(expr)>0,] #选取非零基因 
  expr <- as.matrix(expr) 
  head(expr)
  
  gsva.res <- gsva(expr, genesets, method="ssgsea", parallel.sz=8) 
  saveRDS(gsva.res, paste0('./GSVA/',i,"_C5.gsva.res.rds")) 
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.csv(gsva.df,paste0('./GSVA/',i, "_C5.gsva_res.csv"), row.names = F)
}
  
###############   cell cycle
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("CELL_CYCLE",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]
  
pdf(paste0('./GSVA_cycle_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_cols = F)
dev.off()
}  
###############  death 
dir.create('./GSVA_DEATH_plot/')
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("CELL_DEATH",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]

  pdf(paste0('./GSVA_DEATH_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = F)
  dev.off()
}  

###############  cytokine
dir.create('./GSVA_cytokin_plot/')
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("CYTOKIN",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]
  
  pdf(paste0('./GSVA_cytokin_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = F)
  dev.off()
}  
###############  INFLAMMA
dir.create('./GSVA_INFLAMMA_plot/')
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("INFLAMMA",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]
  
  pdf(paste0('./GSVA_INFLAMMA_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = F)
  dev.off()
} 
###############  INFLAMMA
dir.create('./GSVA_CHEMOKIN_plot/')
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("CHEMOKIN",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]
  
  pdf(paste0('./GSVA_CHEMOKIN_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = F)
  dev.off()
}  
###############  activate 
dir.create('./GSVA_ACTIVATE_plot/')
for(i in levels(lung@active.ident)){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl("ACTIVATE",Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]
  
  pdf(paste0('./GSVA_ACTIVATE_plot/',i, ".pdf"),width = 15,height = 15)    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_cols = F)
  dev.off()
}  


###############   DNA REPAIR
.libPaths("/media/user/sdf/R.lib/library/cmq")
library(pheatmap)
setwd('/media/user/sdh/cmq_wj')
# NAME1<-c('DNA_DAMAGE','DNA_REPAIR','CELLULAR_SENESCENCE','VIRUS','APOPTOTIC')
NAME1<-c('CELLULAR_SENESCENCE',
         
         'GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE',
         'GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS',
         'GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE',
         
         'GOBP_DNA_REPAIR',
         'GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR',
         
         'GOBP_RESPONSE_TO_VIRUS',
         'GOBP_CELLULAR_RESPONSE_TO_VIRUS',
         
         'GOBP_CELL_CYCLE',
         'GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE',
         
         
         'GOBP_APOPTOTIC_SIGNALING_PATHWAY',
         'GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY'
         
         )
CELL_TYPE<-c("NK_cell" ,"T_cell" ,"B_cell" ,"Neutrophil" ,"Monocyte_CCL3" ,"Monocyte_Ly6c_low" , 
              "Macrophage_Ki67" ,"Monocyte_Ly6c_high" , "Macrophage_alveolar")

for (NAME in NAME1) {
  if(!dir.exists(paste0('./GSVA_',NAME))){
  dir.create(paste0('./GSVA_',NAME))}
for(i in CELL_TYPE){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  GSVA_SE<-subset(GSVA,grepl(NAME,Genesets))
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-1]

  pdf(paste0('./GSVA_',NAME,'/',NAME,'_',i, ".pdf"),width = 15,height = length(GSVA_SE[,1])/5 +3 )    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_cols = F)
  dev.off()
}  
}

#############  select
NAME1<-c(
         'GOBP_CELLULAR_SENESCENCE',
         'GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE',
         
         'GOBP_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE',
         'GOBP_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS',
         'GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE',
         
         'GOBP_DNA_REPAIR',
         'GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR',
         
         'GOBP_RESPONSE_TO_VIRUS',
         'GOBP_CELLULAR_RESPONSE_TO_VIRUS',
         
         'GOBP_CELL_CYCLE',
         'GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE',
         
         
         'GOBP_APOPTOTIC_SIGNALING_PATHWAY',
         'GOBP_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY'
         
)
CELL_TYPE<-c("NK_cell" ,"T_cell" ,"B_cell" ,"Neutrophil" ,"Monocyte_CCL3" ,"Monocyte_Ly6c_low" , 
             "Macrophage_Ki67" ,"Monocyte_Ly6c_high" , "Macrophage_alveolar")
NAME2<-cbind(NAME1,c(11:23));colnames(NAME2)<-c('Genesets','rank')

for(i in CELL_TYPE){
  GSVA<-read.csv(paste0('./GSVA/',i, "_C5.gsva_res.csv"))
  colnames(GSVA)
  
  # GSVA_SE<-GSVA[GSVA$Genesets %in% NAME1,]
  GSVA_SE<-merge(GSVA,NAME2,by='Genesets')
  
  # View(file[order(file$Code,file$Score),])
  
  GSVA_SE<-GSVA_SE[order(GSVA_SE$rank),]

  
  row.names(GSVA_SE)<-GSVA_SE$Genesets;GSVA_SE<-GSVA_SE[,-c(1,6)]
  # if(i == 'Monocyte_CCL3'){GSVA_SE<-GSVA_SE[,c(1:2)]}

  pdf(paste0('/media/user/sdh/cmq_wj/GSVA_se/',i, ".pdf"),width = 15,height = length(GSVA_SE[,1])/5 +3 )    
  pheatmap(GSVA_SE, show_colnames = T,
           scale = "row",angle_col = "45",
           cellwidth = 35, cellheight = 12, 
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           cluster_cols = F,cluster_rows = F)
  dev.off()
}  

