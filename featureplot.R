lung<-readRDS('/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds')

FeaturePlot(lung,'Lcn2')
FeaturePlot(lung,'Fcgr3') #CD16
FeaturePlot(lung,'Fcgr2b') #CD32
FeaturePlot(lung,'Itgam') #Cd11b
FeaturePlot(lung,'Ly6c1') #Cd11b
FeaturePlot(lung,'Ly6c2') 
FeaturePlot(lung,'Fcgr1')  #CD64 has emerged as a marker for monocyte-derived (mo) DCs
FeaturePlot(lung,'Cd80') 
FeaturePlot(lung,'Cd86') 
FeaturePlot(lung,'Cd40') 
FeaturePlot(lung,'Cd18') 
DimPlot(lung,label = T)



pdf('ccl3_maeker_ren.cell.pdf',width=10,height=20)
VlnPlot(lung,c('Il1b', 'Ccl3', 'Cxcl8',
'Ccl3l1', 'Plaur', 'Cxcl4',
'S100a8', 'S100a9', 'S100a12', 'Fcer1g',
'Spi1', 'Fosb', 'Fos', 'Klf4', 'Klf6',
'Cebpd', 'Jund', 'Zeb2', 'Mafb',
'Cebpb','Atf3', 'Klf10',
'Junb', 'Ets2', 'Znf385a', 'Egr1',
'Creb5', 'Bach1', 'Hif1a'),pt.size = 0,ncol=4)
dev.off()
