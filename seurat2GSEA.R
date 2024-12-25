setwd('/media/user/sdh/cmq_wj/GSEA')
if(F){
# source activate cmq
# export R_LIBS=/media/user/sdf/R.lib/library/wj/reg  
# cd /media/user/sdh/cmq_wj
# nohup Rscript seurat2GSEA.R &
  }
# .libPaths("/media/user/sdf/R.lib/library/wj/reg")
# export R_LIBS=/media/user/sdf/R.lib/library/wj/reg
library(Seurat)
lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds')

# spilt cell group
for(i in levels(lung@active.ident)){
  
  clust.name <- c(rownames(lung@meta.data)[grep(i,lung@active.ident)])
  
  if(i=='NK_cell'){NK_cell<- subset(lung,cells=clust.name)}
  if(i=='T_cell'){T_cell<- subset(lung,cells=clust.name)}
  if(i=='B_cell'){B_cell<- subset(lung,cells=clust.name)}
  if(i=='Neutrophil'){Neutrophil<- subset(lung,cells=clust.name)}
  if(i=='Monocyte_CCL3'){Monocyte_CCL3<- subset(lung,cells=clust.name)}  
  if(i=='Monocyte_Ly6c_low'){Monocyte_Ly6c_low<- subset(lung,cells=clust.name)}  
  if(i=='Macrophage_Ki67'){Macrophage_Ki67<- subset(lung,cells=clust.name)}  
  if(i=='Monocyte_Ly6c_high'){Monocyte_Ly6c_high<- subset(lung,cells=clust.name)}
  if(i=='Macrophage_alveolar'){Macrophage_alveolar<- subset(lung,cells=clust.name)}
  
}
# seuo<-NK_cell
# test1<-'G6';test2<-'G5'
# name<-'NK_cell'
# group<-'orig.ident'

Seurat2GSEA<-  function(seuo,group,test1,test2,name){
  Idents(seuo) <- group
  seuo <- subset(seuo,idents=c(test1,test2))  
  seuo@assays$RNA@counts-> counts   # 是用count还是data？ 
  expr_data <- as.matrix(counts)  # 数据量太大怎么办？可以取 pseudocell 啊，见参考。
  count_data<-as.data.frame(rbind(c('NAME','Description',colnames(expr_data)),
                              cbind(rownames(expr_data),'na',expr_data))) 
  con1<-file(paste(name,test1,test2,'expr.gct.txt',sep = "_"),open='w')
  write('#1.2',con1)
  write(paste(length(rownames(expr_data)),'\t',length(colnames(expr_data)),sep = ''),con1)
  tep<-''
  for (hh in 1:length(count_data[,1])){ #每行
    tep<-count_data[hh,1]
    for (uu in 2:length(count_data[1,])){ #每列
      tep<-paste(tep,count_data[hh,uu],sep = '\t')

    }
    write(tep,con1)
    }

  close(con1)
  # 整理GSEA需要的表达谱格式
  write.table(rbind(c('NAME','Description',colnames(expr_data)),
                    cbind(rownames(expr_data),'na',expr_data)),
              file=paste(name,test1,test2,'expr.gct.txt',sep = "_"),quote=F,sep='\t',col.names=F,row.names=F)

  # 整理GSEA需要的分组格式
  pheno<-as.character(seuo@meta.data[,group])  # 注意顺序
  con<-file(paste(name,test1,test2,'pheno.cls',sep = "_"),open='w')
  write(paste(length(pheno),'2','1',sep = '\t'),con)
  write(paste('#',test1,test2,sep = '\t'),con)
  classes<-pheno[1]
  for (i in 2:length(pheno)){
    classes<-paste(classes,pheno[i],sep = '\t')
  }
  write(classes,con)
  close(con)
}


group<-'orig.ident'
for(name in levels(lung@active.ident) ){
  if(name ==  "NK_cell"){
    Seurat2GSEA( NK_cell,group,'G6','G5',name)
    Seurat2GSEA( NK_cell,group,'G7','G5',name)
    Seurat2GSEA( NK_cell,group,'G8','G6',name)
    Seurat2GSEA( NK_cell,group,'G8','G7',name)
  }
  if(name ==  "T_cell"){
    Seurat2GSEA( T_cell,group,'G6','G5',name)
    Seurat2GSEA( T_cell,group,'G7','G5',name)
    Seurat2GSEA( T_cell,group,'G8','G6',name)
    Seurat2GSEA( T_cell,group,'G8','G7',name)
  }
  
  if(name ==  "B_cell"){
    Seurat2GSEA( B_cell,group,'G6','G5',name)
    Seurat2GSEA( B_cell,group,'G7','G5',name)
    Seurat2GSEA( B_cell,group,'G8','G6',name)
    Seurat2GSEA( B_cell,group,'G8','G7',name)
  }
  
  if(name ==  "Neutrophil"){
    Seurat2GSEA( Neutrophil,group,'G6','G5',name)
    Seurat2GSEA( Neutrophil,group,'G7','G5',name)
    Seurat2GSEA( Neutrophil,group,'G8','G6',name)
    Seurat2GSEA( Neutrophil,group,'G8','G7',name)
  }
  
  if(name ==  "Monocyte_CCL3"){
    Seurat2GSEA( Monocyte_CCL3,group,'G6','G5',name)
    Seurat2GSEA( Monocyte_CCL3,group,'G7','G5',name)
    Seurat2GSEA( Monocyte_CCL3,group,'G8','G6',name)
    Seurat2GSEA( Monocyte_CCL3,group,'G8','G7',name)
  }
  
  if(name ==  "Monocyte_Ly6c_low"){
    Seurat2GSEA( Monocyte_Ly6c_low,group,'G6','G5',name)
    Seurat2GSEA( Monocyte_Ly6c_low,group,'G7','G5',name)
    Seurat2GSEA( Monocyte_Ly6c_low,group,'G8','G6',name)
    Seurat2GSEA( Monocyte_Ly6c_low,group,'G8','G7',name)
  }
  if(name ==  "Macrophage_Ki67"){
    Seurat2GSEA( Macrophage_Ki67,group,'G6','G5',name)
    Seurat2GSEA( Macrophage_Ki67,group,'G7','G5',name)
    Seurat2GSEA( Macrophage_Ki67,group,'G8','G6',name)
    Seurat2GSEA( Macrophage_Ki67,group,'G8','G7',name)
  }
  
  if(name ==  "Monocyte_Ly6c_high"){
    Seurat2GSEA( Monocyte_Ly6c_high,group,'G6','G5',name)
    Seurat2GSEA( Monocyte_Ly6c_high,group,'G7','G5',name)
    Seurat2GSEA( Monocyte_Ly6c_high,group,'G8','G6',name)
    Seurat2GSEA( Monocyte_Ly6c_high,group,'G8','G7',name)
  }
  
  if(name ==  "Macrophage_alveolar"){
    Seurat2GSEA( Macrophage_alveolar,group,'G6','G5',name)
    Seurat2GSEA( Macrophage_alveolar,group,'G7','G5',name)
    Seurat2GSEA( Macrophage_alveolar,group,'G8','G6',name)
    Seurat2GSEA( Macrophage_alveolar,group,'G8','G7',name)
  }
}

