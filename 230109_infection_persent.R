G7_seurat<-readRDS('/media/user/sdh/cmq_wj/G_rds/G7.rds')
G8_seurat<-readRDS('/media/user/sdh/cmq_wj/G_rds/G8.rds')


aa1<-as.data.frame(G7_seurat@meta.data[["COV2"]]) 
aa2<-as.data.frame(G7_seurat@active.ident) 
aa<-cbind(aa1,aa2)

aa3<-as.data.frame(G8_seurat@meta.data[["COV2"]]) 
aa4<-as.data.frame(G8_seurat@active.ident) 
aa5<-cbind(aa3,aa4)


G7<-as.data.frame(unique(aa$`G7_seurat@active.ident`))
colnames(G7)<-c('cluster')
for (i in c(1:9)) {
  G7$G7_infe[i]<-sum(subset(aa,aa$`G7_seurat@active.ident`==G7[i,1])$`G7_seurat@meta.data[["COV2"]]`>0)
  G7$G7_uninfe[i]<-sum(subset(aa,aa$`G7_seurat@active.ident`==G7[i,1])$`G7_seurat@meta.data[["COV2"]]`<= 0)
  
  G7$G8_infe[i]<-sum(subset(aa5,aa5$`G8_seurat@active.ident`==G7[i,1])$`G8_seurat@meta.data[["COV2"]]`>0)
  G7$G8_uninfe[i]<-sum(subset(aa5,aa5$`G8_seurat@active.ident`==G7[i,1])$`G8_seurat@meta.data[["COV2"]]`<= 0)
}
G7$G7_pre<-G7$G7_infe/(G7$G7_infe+G7$G7_uninfe)*100
G7$G8_pre<-G7$G8_infe/(G7$G8_infe+G7$G8_uninfe)*100
write.csv(G7,'G7_G8_infected_cell_num.csv')
