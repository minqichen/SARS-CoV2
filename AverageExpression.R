# AverageExpression
neu<-readRDS('/media/user/sdh/cmq_wj/G_rds/Neutrophil.rds')
Idents(neu)<-neu@meta.data[["orig.ident"]]
neu.averages <- AverageExpression(neu, slot = 'data', return.seurat = T)
head(neu.averages@assays$RNA@data)
write.csv(neu.averages@assays$RNA@data, file = "NEU.averages_assays_RNA_data.csv")
