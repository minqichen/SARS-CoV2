celltypes<-c("NK_cell" ,"T_cell", "B_cell","Neutrophil" ,        
             "Monocyte_CCL3","Monocyte_Ly6c_low","Macrophage_Ki67",
             "Monocyte_Ly6c_high","Macrophage_alveolar")
color.heatmap = "Reds"
################################################################################
###########    G7_G8 
G7_G8<-c('TWEAK','IFN-II','TNF','COMPLEMENT','FASLG','PARs',
         'MIF','SEMA3', 'OSM','ANNEXIN','CXCL','IL1')
G7_G8<-c('FASLG','CCL','IL16','IFN-II',
         'VISFATIN','TNF', 'GALECTIN','APRIL',
         'CSF','LIGHT','TGFb','PROS','IGF',
         'GDF','GALECTIN', 'CSF','CCL',
         'TWEAK','VISFATIN','APRIL',
         'LIGHT','IL16','GDF',
         'PROS','IGF','TGFb'
)
for (signaling in G7_G8){
  title.name=paste('G7_G8',signaling,sep = '_')
  
  if (signaling %in% cc.g7@netP[["pathways"]]){
    aa1<-cc.g7@netP[["prob"]][,,signaling];colnames(aa1)<-paste0('G7_',colnames(aa1))}
  
  if (signaling %in% cc.g8@netP[["pathways"]]){
    aa2<-cc.g8@netP[["prob"]][,,signaling];colnames(aa2)<-paste0('G8_',colnames(aa2))}
  if (!signaling %in% cc.g7@netP[["pathways"]]){
    aa1 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa1)<-paste0('G7_',celltypes)
    rownames(aa1)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g8@netP[["pathways"]]){
    aa2 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa2)<-paste0('G8_',celltypes)
    rownames(aa2)<-paste0(celltypes)
  }
  
  aa_uu<-cbind(aa1,aa2)
  # pheatmap::pheatmap(aa_uu,cluster_rows = F,cluster_cols = F)
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/de_heatmap/',title.name,'_HEATMAP.pdf'),
      width = 6,height = 4)
  pheatmap::pheatmap(aa_uu,scale = 'none',
                     cluster_rows = F,cluster_cols = F,
                     color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                     main=title.name,
                     gaps_col = c(9), border=FALSE)
  dev.off()
}


################################################################################
###########    G6_G8 
G6_G8<-c('GDF','SEMA3','GAS','IGF','MIF','IL2','ANNEXIN',
         'ANGPTL','COMPLEMENT','CCL','LIGHT','TNF','GALECTIN'
)
G6_G8<-c('MIF', 'GDF','FASLG','PROS','SPP1','TGFb',
         'APRIL', 'SPP1','FASLG','PROS','APRIL',
         'IFN-II','CSF','OSM','ANGPTL','CXCL',
         'TWEAK', 'TNF', 'IL1','VISFATIN','ANNEXIN'
)
  
for (signaling in G6_G8){
  title.name=paste('G6_G8',signaling,sep = '_')
  
  if (signaling %in% cc.g6@netP[["pathways"]]){
    aa1<-cc.g6@netP[["prob"]][,,signaling];colnames(aa1)<-paste0('G6_',colnames(aa1))}
  
  if (signaling %in% cc.g8@netP[["pathways"]]){
    aa2<-cc.g8@netP[["prob"]][,,signaling];colnames(aa2)<-paste0('G8_',colnames(aa2))}
  if (!signaling %in% cc.g6@netP[["pathways"]]){
    aa1 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa1)<-paste0('G6_',celltypes)
    rownames(aa1)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g8@netP[["pathways"]]){
    aa2 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa2)<-paste0('G8_',celltypes)
    rownames(aa2)<-paste0(celltypes)
  }
  
  aa_uu<-cbind(aa1,aa2)
  # pheatmap::pheatmap(aa_uu,cluster_rows = F,cluster_cols = F)
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/de_heatmap/',title.name,'_HEATMAP.pdf'),
      width = 6,height = 4)
  pheatmap::pheatmap(aa_uu,scale = 'none',
                     cluster_rows = F,cluster_cols = F,
                     color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                     main=title.name,
                     gaps_col = c(9), border=FALSE)
  dev.off()
}

################################################################################
###########    G5_G6 
G5_G6<-c('CD40','GALECTIN','LIGHT','GDF','IL16','TGFb','GRN',
         'CCL','ANNEXIN','ANGPTL','APRIL','PROS','FLT3',
         'IL2','MIF','GAS','IGF','SEMA3')
G5_G6<-c('PROS',
         'VISFATIN',
         'GDF',
         'VEGF',
         'COMPLEMENT',
         'ANGPTL',
         'CXCL',
         'OSM',
         'IL1',
         'CSF','SPP1',
         'FASLG','TGFb','GRN',
         'VISFATIN',
         'CSF',
         'OSM',
         'IL1',
         'CXCL',
         'IFN-II',
         'TNF','SPP1',
         'FASLG'
         
)
for (signaling in G5_G6){
  title.name=paste('G5_G6',signaling,sep = '_')
  
  if (signaling %in% cc.g5@netP[["pathways"]]){
    aa1<-cc.g5@netP[["prob"]][,,signaling];colnames(aa1)<-paste0('G6_',colnames(aa1))}
  
  if (signaling %in% cc.g6@netP[["pathways"]]){
    aa2<-cc.g6@netP[["prob"]][,,signaling];colnames(aa2)<-paste0('G8_',colnames(aa2))}
  if (!signaling %in% cc.g5@netP[["pathways"]]){
    aa1 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa1)<-paste0('G6_',celltypes)
    rownames(aa1)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g6@netP[["pathways"]]){
    aa2 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa2)<-paste0('G8_',celltypes)
    rownames(aa2)<-paste0(celltypes)
  }
  
  aa_uu<-cbind(aa1,aa2)
  # pheatmap::pheatmap(aa_uu,cluster_rows = F,cluster_cols = F)
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/de_heatmap/',title.name,'_HEATMAP.pdf'),
      width = 6,height = 4)
  pheatmap::pheatmap(aa_uu,scale = 'none',
                     cluster_rows = F,cluster_cols = F,
                     color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                     main=title.name,
                     gaps_col = c(9), border=FALSE)
  dev.off()
}

######################################################
###########    G5_G7 
G5_G7<-c('CD40','GDF','PROS','IL16','ANNEXIN','FLT3',
         'PARs','TWEAK','FASLG','CSF','TGFb',
         'CCL','COMPLEMENT', 'IFN-II','TNF'
         
)
G5_G7<-c('VEGF','OSM','APRIL','VEGF',
         'GRN','TGFb','APRIL',
         'OSM','ANGPTL','IL1',
         'GALECTIN','LIGHT',
         'VISFATIN','IGF','GALECTIN',
         'IL1','LIGHT','CSF','FASLG',
         'IGF','IFN-II'
         
)

for (signaling in G5_G7){
  title.name=paste('G5_G7',signaling,sep = '_')
  
  if (signaling %in% cc.g5@netP[["pathways"]]){
    aa1<-cc.g5@netP[["prob"]][,,signaling];colnames(aa1)<-paste0('G6_',colnames(aa1))}
  
  if (signaling %in% cc.g7@netP[["pathways"]]){
    aa2<-cc.g7@netP[["prob"]][,,signaling];colnames(aa2)<-paste0('G8_',colnames(aa2))}
  if (!signaling %in% cc.g5@netP[["pathways"]]){
    aa1 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa1)<-paste0('G6_',celltypes)
    rownames(aa1)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g7@netP[["pathways"]]){
    aa2 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa2)<-paste0('G8_',celltypes)
    rownames(aa2)<-paste0(celltypes)
  }
  
  aa_uu<-cbind(aa1,aa2)
  # pheatmap::pheatmap(aa_uu,cluster_rows = F,cluster_cols = F)
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/de_heatmap/',title.name,'_HEATMAP.pdf'),
      width = 6,height = 4)
  pheatmap::pheatmap(aa_uu,scale = 'none',
                     cluster_rows = F,cluster_cols = F,
                     color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                     main=title.name,
                     gaps_col = c(9), border=FALSE)
  dev.off()
}



######################################################
###########    G5_G6_G7_G8 
G5_G6_G7_G8<-c('ANNEXIN', 'IFN-II','TNF')
G5_G6_G7_G8<-c('PARs')
for (signaling in G5_G6_G7_G8){
  title.name=paste('G5_G6_G7_G8',signaling,sep = '_')
  
  if (signaling %in% cc.g5@netP[["pathways"]]){
    aa1<-cc.g5@netP[["prob"]][,,signaling];colnames(aa1)<-paste0('G5_',colnames(aa1))}
  
  if (signaling %in% cc.g6@netP[["pathways"]]){
    aa2<-cc.g6@netP[["prob"]][,,signaling];colnames(aa2)<-paste0('G6_',colnames(aa2))}
  
  if (signaling %in% cc.g7@netP[["pathways"]]){
    aa3<-cc.g7@netP[["prob"]][,,signaling];colnames(aa3)<-paste0('G7_',colnames(aa3))}
  
  if (signaling %in% cc.g8@netP[["pathways"]]){
    aa4<-cc.g8@netP[["prob"]][,,signaling];colnames(aa4)<-paste0('G8_',colnames(aa4))}
  
  
  if (!signaling %in% cc.g5@netP[["pathways"]]){
    aa1 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa1)<-paste0('G5_',celltypes)
    rownames(aa1)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g6@netP[["pathways"]]){
    aa2 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa2)<-paste0('G6_',celltypes)
    rownames(aa2)<-paste0(celltypes)
  }
  
  if (!signaling %in% cc.g7@netP[["pathways"]]){
    aa3 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa3)<-paste0('G7_',celltypes)
    rownames(aa3)<-paste0(celltypes)
  }
  if (!signaling %in% cc.g8@netP[["pathways"]]){
    aa4 <-matrix(0,nrow = 9, ncol = 9)
    colnames(aa4)<-paste0('G8_',celltypes)
    rownames(aa4)<-paste0(celltypes)
  }
  
  aa_uu<-cbind(aa1,aa2,aa3,aa4)
  # pheatmap::pheatmap(aa_uu,cluster_rows = F,cluster_cols = F)
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/de_heatmap/',title.name,'_HEATMAP.pdf'),
      width = 12,height = 4)
  pheatmap::pheatmap(aa_uu,scale = 'none',
                     cluster_rows = F,cluster_cols = F,
                     color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                     main=title.name,
                     gaps_col = c(9,18,27), border=FALSE)
  dev.off()
}

###########   5_6_7_8 BUBBLE
levels(mg.cc@idents$joint)

# Neutrophil
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_Neutrophil.pdf'),
    width = 18,height = 6)
netVisual_bubble(mg.cc, sources.use = c(4), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45)
dev.off()
# TNF

pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_TNF.pdf'),
    width = 18,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5,7:9), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,signaling  = 'TNF')
dev.off()
# IFN-II
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_IFN-II.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(1,7), targets.use = c(4:9),  
                 comparison = c(1,2,3,4), angle.x = 45,signaling  = 'IFN-II')
dev.off()
# ANNEXIN
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_ANNEXIN.pdf'),
    width = 18,height = 3)
netVisual_bubble(mg.cc, sources.use = c(1,4:9), targets.use = c(4:9),  
                 comparison = c(1,2,3,4), angle.x = 45,signaling  = 'ANNEXIN')
dev.off()
# select
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220720.pdf'),
    width = 5,height = 4)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(4:6,8:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
dev.off()


pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220720_1.pdf'),
    width = 5,height = 3)
netVisual_bubble(mg.cc, sources.use = c(4), targets.use = c(5:6,8),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
dev.off()
#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220720_2.pdf'),
    width = 5,height = 3)
netVisual_bubble(mg.cc, sources.use = c(6), targets.use = c(4,5,9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
dev.off()
#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220720_3.pdf'),
    width = 5,height = 3)
netVisual_bubble(mg.cc, sources.use = c(8), targets.use = c(6,7),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
dev.off()

#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220720_4.pdf'),
    width = 5,height = 3)
netVisual_bubble(mg.cc, sources.use = c(4,5,8), targets.use = c(6,9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
dev.off()
#
netVisual_bubble(mg.cc, sources.use = c(4,5,8), targets.use = c(9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('TNF','ANNEXIN'))
## 220728
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_1.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('IL1'))
dev.off()
## 220728
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_2.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('COMPLEMENT'))
dev.off()
#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_3.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('COMPLEMENT','IL1'))
dev.off()

#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_4.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('LIGHT','VISFATIN'))
dev.off()

#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_5.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('CCL'))
dev.off()
#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_6.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('GALECTIN','CSF'))
dev.off()

#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_220728_7.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('COMPLEMENT','CSF','VISFATIN','TNF'))
dev.off()



#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_22824_1.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('COMPLEMENT','CSF','VISFATIN','TNF'))
dev.off()


#
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_select_1127.pdf'),
    width = 8,height = 3)
netVisual_bubble(mg.cc, sources.use = c(1:9), targets.use = c(1:9),  
                 comparison = c(1,2,3,4), angle.x = 45,
                 signaling  = c('PARs'))
dev.off()

netAnalysis_contribution(cc.g5, signaling = 'GALECTIN')

