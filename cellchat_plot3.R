levels(mg.cc@idents$joint)
gn1<-'Monocyte_Ly6c_high'
gn1<-'Monocyte_CCL3'
gn1<-'Monocyte_Ly6c_low'
gn1<-'T_cell'

gg1 <-netAnalysis_signalingChanges_scatter(mg.cc, idents.use = gn1,
                                           comparison = c(1, 2), 
                                           signaling.exclude = "MIF")
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G5_G6_',gn1,'_signalingChanges.pdf'),width = 7.5,height = 6)
patchwork::wrap_plots(plots = list(gg1))
dev.off()

gg2 <-netAnalysis_signalingChanges_scatter(mg.cc, idents.use = gn1,
                                           comparison = c(3, 4), 
                                           signaling.exclude = "MIF")
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G7_G8_',gn1,'_signalingChanges.pdf'),width = 7.5,height = 6)
patchwork::wrap_plots(plots = list(gg2))
dev.off()

####### 放大
library(ggrepel)
library(ggplot2)
# rr<-gg1[["data"]]
rr<-gg2[["data"]]

RR1<-rr[which(rr$outgoing < 0.01),]
RR1<-RR1[which(-0.01 < RR1$outgoing),]

RR1<-RR1[which( RR1$incoming<0.01 ),]
RR1<-RR1[which( RR1$incoming> -0.01 ),]
ggplot(RR1,aes(x=outgoing,y=incoming,
               shape=specificity.out.in ,fill=specificity,colour=specificity))+
  geom_point(size=3)+
  scale_shape_manual(values = c(21,22,24,23))+
  scale_color_manual(values=c( "grey10","#F8766D","#00BFC4"))+
  scale_fill_manual(values=alpha(c( "grey10","#F8766D","#00BFC4"),0.5))+
  # geom_label_repel(aes(label=labels))+
  geom_text_repel(aes(outgoing, incoming, label=labels,arrow = T))+
  theme_classic(base_size = 16)+ 
  geom_hline(aes(yintercept=0.0), linetype="dashed",colour='grey')+ 
  geom_vline(aes(xintercept=0.0), linetype="dashed",colour='grey')+
  labs(x="Differential outgoing interaction strength",
       y="Differential incoming interaction strength")

# ggsave(paste0('/media/user/sdh/cmq_wj/cellChat/G5_G6_',gn1,'_signalingChanges_gg.pdf'),
#        width =8,
#        height = 5)
ggsave(paste0('/media/user/sdh/cmq_wj/cellChat/G7_G8_',gn1,'_signalingChanges_gg.pdf'),
       width =8,
       height = 5)

#######  specific pathway heatmap
for (pathways.show in rownames(subset(rr,rr$specificity =='g5 specific'))){
# pathways.show<-'CD40'
NH1<-netVisual_heatmap(cc.g5, signaling = pathways.show, color.heatmap = "Reds")
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_heatmap_',pathways.show,'_g5.pdf'))
print(NH1)
dev.off()
}


for (pathways.show in rownames(subset(rr,rr$specificity =='g6 specific'))){
  # pathways.show<-'CD40'
  NH2<-netVisual_heatmap(cc.g6, signaling = pathways.show, color.heatmap = "Reds")
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_heatmap_',pathways.show,'_g6.pdf'))
  print(NH2);dev.off()
}

gn1<-'Monocyte_Ly6c_high'
gn1<-'Monocyte_CCL3'
gn1<-'Monocyte_Ly6c_low'
gg2 <-netAnalysis_signalingChanges_scatter(mg.cc, idents.use = gn1,
                                           comparison = c(3, 4), 
                                           signaling.exclude = "MIF")
rr<-gg2[["data"]]


for (pathways.show in rownames(subset(rr,rr$specificity =='g7 specific'))){
  # pathways.show<-'CD40'
  NH2<-netVisual_heatmap(cc.g7, signaling = pathways.show, color.heatmap = "Reds")
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_heatmap_',pathways.show,'_g7.pdf'))
  print(NH2);dev.off()
}

for (pathways.show in rownames(subset(rr,rr$specificity =='g8 specific'))){
  # pathways.show<-'CD40'
  NH2<-netVisual_heatmap(cc.g8, signaling = pathways.show, color.heatmap = "Reds")
  pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_heatmap_',pathways.show,'_g8.pdf'))
  print(NH2);dev.off()
}




pdf('/media/user/sdh/cmq_wj/cellChat/G5-8_chord.pdf')
pathway.union <- union(cc.g5@netP$pathways, cc.g8@netP$pathways)
netVisual_aggregate(cc.g5, signaling = pathway.union, layout = "chord")
pathway.union <- union(cc.g6@netP$pathways, cc.g8@netP$pathways)
netVisual_aggregate(cc.g6, signaling = pathway.union, layout = "chord")
pathway.union <- union(cc.g7@netP$pathways, cc.g8@netP$pathways)
netVisual_aggregate(cc.g7, signaling = pathway.union, layout = "chord")
pathway.union <- union(cc.g8@netP$pathways, cc.g8@netP$pathways)
netVisual_aggregate(cc.g8, signaling = pathway.union, layout = "chord")
dev.off() 






