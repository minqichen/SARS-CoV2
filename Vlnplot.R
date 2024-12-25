c('Nr4a1', 'Tnf','Tnfrsf1b','Ifng',
                              'Tnfrsf1a','Anxa1' , 'Fpr2', 'Fpr1')
FeaturePlot(lung,reduction ='harmony_umap',label = T,ncol = 4,
            split.by =  'orig.ident',
            features ='Anxa1' )
ggsave('FeaturePlot_CELLCHAT.pdf',height = 10,width = 18)

VlnPlot(lung,group.by = 'orig.ident', pt.size = 0.1,
        idents =c("Monocyte_Ly6c_high"  ),
        features = c('Tnf'),ncol=4)


for( id in c("Neutrophil","Monocyte_CCL3","Monocyte_Ly6c_high")){ 
  pdf(paste0('/media/user/sdh/cmq_wj/vln/',id,'_vlnplot.pdf') ,
        width = 5,height = 3)

for (tt in c('Nr4a1', 'Tnf','Anxa1')){   
    aa<-VlnPlot(lung,group.by = 'orig.ident', pt.size = 0.1,
            idents =id,
            features = tt)
    print(aa)
  }
  dev.off()
}


for( id in levels(mg.cc@idents$joint)){
  pdf(paste0('/media/user/sdh/cmq_wj/vln_221121/',id,'_vlnplot.pdf') ,
      width = 15,height = 9)
for (tt in c('Cd68')){  
    
    aa<-VlnPlot(lung,group.by = 'orig.ident', pt.size = 0.1,
                idents =id,
                features = c('Cd68','Itgam', 
                             'Cd80', 'Cd86',
                             'Fcgr1','Fcgr4','Fcgr2b',
                             'Cd163','Mrc1',
                             'Tnf','Il2','Il10'))
    print(aa)
  }
  dev.off()
}
c('Fcna','Spp1','Fabp4')


for( id in levels(mg.cc@idents$joint)){
  pdf(paste0('/media/user/sdh/cmq_wj/vln_221121/',id,'_vlnplot_1127.pdf') ,
      width = 15,height = 9)
  for (tt in c('Cd68')){  
    
    aa<-VlnPlot(lung,group.by = 'orig.ident', pt.size = 0.1,
                idents =id,
                features = c('Fcna','Spp1','Fabp4'))
    print(aa)
  }
  dev.off()
}
#########
for (ii in c('Cd68','Itgam', 
            'Cd80', 'Cd86',
            'Fcgr1','Fcgr4','Fcgr2b',
            'Cd163','Mrc1',
            'Tnf','Il2','Il10') ){
pdf(paste0('/media/user/sdh/cmq_wj/vln_221121/',ii,'_FeaturePlot.pdf') ,
    width = 20,height = 5)
aa<-FeaturePlot(lung,reduction ='harmony_umap',label = T,ncol = 2,
            split.by =  'orig.ident',
            features =ii)
print(aa)
dev.off()
}


for (ii in c('Fcna','Spp1','Fabp4') ){
  pdf(paste0('/media/user/sdh/cmq_wj/vln_221121/',ii,'_FeaturePlot.pdf') ,
      width = 20,height = 5)
  aa<-FeaturePlot(lung,reduction ='harmony_umap',label = T,ncol = 2,
                  split.by =  'orig.ident',
                  features =ii)
  print(aa)
  dev.off()
}


[1] "NK_cell"             "T_cell"             
[3] "B_cell"              "Neutrophil"         
[5] "Monocyte_CCL3"       "Monocyte_Ly6c_low"  
[7] "Macrophage_Ki67"     "Monocyte_Ly6c_high" 
[9] "Macrophage_alveolar"