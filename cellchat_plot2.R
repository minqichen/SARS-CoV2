cc.g5@netP$pathways
pathways.show<-"GALECTIN" 
# G5
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G5_',pathways.show,'_plot.pdf'))
netVisual_heatmap(cc.g5, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cc.g5, signaling = pathways.show)
netVisual_aggregate(cc.g5, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cc.g5, signaling = pathways.show, layout = "chord")
dev.off()

# G6
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G6_',pathways.show,'_plot.pdf'))
netVisual_heatmap(cc.g6, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cc.g6, signaling = pathways.show)
netVisual_aggregate(cc.g6, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cc.g6, signaling = pathways.show, layout = "chord")
dev.off()

# G7
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G6_',pathways.show,'_plot.pdf'))
netVisual_heatmap(cc.g7, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cc.g7, signaling = pathways.show)
netVisual_aggregate(cc.g7, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cc.g7, signaling = pathways.show, layout = "chord")
dev.off()

# G8
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G6_',pathways.show,'_plot.pdf'))
netVisual_heatmap(cc.g8, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cc.g8, signaling = pathways.show)
netVisual_aggregate(cc.g8, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cc.g8, signaling = pathways.show, layout = "chord")
dev.off()

#################   
object.list<-xxxx
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
pdf('/media/user/sdh/cmq_wj/cellChat/incom_outcom_perTimeGroup_plot.pdf',width = 6,height = 6)
patchwork::wrap_plots(plots = gg)
dev.off()

###
gg1 <-netAnalysis_signalingChanges_scatter(mg.cc, idents.use = "Monocyte_CCL3", 
                                            signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(mg.cc, idents.use = "T_cell", 
                                            signaling.exclude = c("MIF"))
pdf('/media/user/sdh/cmq_wj/cellChat/G5_G6_Monocyte_CCL3_signalingChanges.pdf',width = 15,height = 6)
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()

rr<-gg1[["data"]]
library(ggrepel)
# RR1<-subset(rr,-0.03 < rr$outgoing )
RR1<-rr[which(rr$outgoing < 0.03),]
RR1<-RR1[which(-0.03 < RR1$outgoing),]

RR1<-RR1[which( RR1$incoming<0.03 ),]

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

ggsave('/media/user/sdh/cmq_wj/cellChat/G5_G6_Monocyte_CCL3_signalingChanges_gg.pdf',
       width =8,
       height = 5)



##############
.libPaths("/media/user/sdf/R.lib/library/cmq");library(ComplexHeatmap)
i=1
# combining all the identified signaling pathways from different datasets 
pathway.union <- unique(c(object.list[[i]]@netP$pathways,
                       object.list[[i+1]]@netP$pathways,
                       object.list[[i+2]]@netP$pathways,
                       object.list[[i+3]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height =7)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height =7)
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height =7)

pdf('/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap.pdf',width = 16,height = 5.5)
draw(ht1 + ht2 + ht3 + ht4, ht_gap = unit(0.5, "cm"))


iht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7, color.heatmap = "GnBu")
iht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7, color.heatmap = "GnBu")
iht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 7, color.heatmap = "GnBu")
iht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 7, color.heatmap = "GnBu")

draw(iht1 + iht2+iht3 + iht4, ht_gap = unit(0.5, "cm"))

aht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 7, color.heatmap = "OrRd")
aht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7, color.heatmap = "OrRd")
aht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 5, height = 7, color.heatmap = "OrRd")
aht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "all", 
                                         signaling = pathway.union, 
                                         title = names(object.list)[i+3], 
                                         width = 5, height = 7, color.heatmap = "OrRd")


draw(aht1 + aht2+aht3 + aht4, ht_gap = unit(0.5, "cm"))
dev.off()

########################################################################
########  netAnalysis_signalingRole_heatmap select right blox
rr1<-as.data.frame(ht1@row_names_param[["labels"]]) 
colnames(rr1)<-'labels'
rr1$count<- ht1@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
rr1$group<-'G5'

rr2<-as.data.frame(ht2@row_names_param[["labels"]]) 
colnames(rr2)<-'labels'
rr2$count<- ht2@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
rr2$group<-'G6'

rr3<-as.data.frame(ht3@row_names_param[["labels"]]) 
colnames(rr3)<-'labels'
rr3$count<- ht3@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
rr3$group<-'G7'

rr4<-as.data.frame(ht4@row_names_param[["labels"]]) 
colnames(rr4)<-'labels'
rr4$count<- ht4@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
rr4$group<-'G8'

rr_all<-rbind(rr1,rr2,rr3,rr4)
library(ggforce)

ggplot(rr_all,aes(x=labels,y=count,fill=group))+
  geom_bar(position = 'dodge',stat='identity',colour='black')+
  scale_fill_brewer(palette = 'Accent')+
  facet_zoom(ylim = c(0, 8))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank())

ggsave(filename = '/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap_right_bar_outcoming.pdf',
       width = 20,height = 3)
##overall
arr1<-as.data.frame(aht1@row_names_param[["labels"]]) 
colnames(arr1)<-'labels'
arr1$count<- aht1@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
arr1$group<-'G5'

arr2<-as.data.frame(aht2@row_names_param[["labels"]]) 
colnames(arr2)<-'labels'
arr2$count<- aht2@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
arr2$group<-'G6'

arr3<-as.data.frame(aht3@row_names_param[["labels"]]) 
colnames(arr3)<-'labels'
arr3$count<- aht3@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
arr3$group<-'G7'

arr4<-as.data.frame(aht4@row_names_param[["labels"]]) 
colnames(arr4)<-'labels'
arr4$count<- aht4@right_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
arr4$group<-'G8'

arr_all<-rbind(arr1,arr2,arr3,arr4)
library(ggforce)

ggplot(arr_all,aes(x=labels,y=count,fill=group))+
  geom_bar(position = 'dodge',stat='identity',colour='black')+
  scale_fill_brewer(palette = 'Accent')+
  facet_zoom(ylim = c(0, 3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank())

ggsave(filename = '/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap_rigaht_bar_overall.pdf',
       width = 20,height = 3)
#########################################################################
########  netAnalysis_signalingRole_heatmap select top blox
tt1<-as.data.frame(ht1@column_names_param[["labels"]]) 
colnames(tt1)<-'labels'
tt1$count<- ht1@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
tt1$group<-'G5'

tt2<-as.data.frame(ht2@column_names_param[["labels"]]) 
colnames(tt2)<-'labels'
tt2$count<- ht2@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
tt2$group<-'G6'

tt3<-as.data.frame(ht3@column_names_param[["labels"]]) 
colnames(tt3)<-'labels'
tt3$count<- ht3@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
tt3$group<-'G7'

tt4<-as.data.frame(ht4@column_names_param[["labels"]]) 
colnames(tt4)<-'labels'
tt4$count<- ht4@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
tt4$group<-'G8'

tt_all<-rbind(tt1,tt2,tt3,tt4)

ggplot(tt_all,aes(x=labels,y=count,fill=group))+
  geom_bar(position = 'dodge',stat='identity',colour='black')+
  scale_fill_brewer(palette = 'Accent')+
  # facet_zoom(ylim = c(0, 3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank())+
  ggtitle ('outcoming')

ggsave(filename = '/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap_top_bar_outcoming.pdf',
       width = 7,height = 3)


##incoming
itt1<-as.data.frame(iht1@column_names_param[["labels"]]) 
colnames(itt1)<-'labels'
itt1$count<- iht1@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
itt1$group<-'G5'

itt2<-as.data.frame(iht2@column_names_param[["labels"]]) 
colnames(itt2)<-'labels'
itt2$count<- iht2@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
itt2$group<-'G6'

itt3<-as.data.frame(iht3@column_names_param[["labels"]]) 
colnames(itt3)<-'labels'
itt3$count<- iht3@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
itt3$group<-'G7'

itt4<-as.data.frame(iht4@column_names_param[["labels"]]) 
colnames(itt4)<-'labels'
itt4$count<- iht4@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
itt4$group<-'G8'

itt_all<-rbind(itt1,itt2,itt3,itt4)

ggplot(itt_all,aes(x=labels,y=count,fill=group))+
  geom_bar(position = 'dodge',stat='identity',colour='black')+
  scale_fill_brewer(palette = 'Accent')+
  # facet_zoom(ylim = c(0, 3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank())+
  ggtitle ('incoming')

ggsave(filename = '/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap_top_bar_incoming.pdf',
       width = 7,height = 3)

##overall
att1<-as.data.frame(aht1@column_names_param[["labels"]]) 
colnames(att1)<-'labels'
att1$count<- aht1@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
att1$group<-'G5'

att2<-as.data.frame(aht2@column_names_param[["labels"]]) 
colnames(att2)<-'labels'
att2$count<- aht2@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
att2$group<-'G6'

att3<-as.data.frame(aht3@column_names_param[["labels"]]) 
colnames(att3)<-'labels'
att3$count<- aht3@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
att3$group<-'G7'

att4<-as.data.frame(aht4@column_names_param[["labels"]]) 
colnames(att4)<-'labels'
att4$count<- aht4@top_annotation@anno_list[["Strength"]]@fun@var_env[["value"]]
att4$group<-'G8'

att_all<-rbind(att1,att2,att3,att4)

ggplot(att_all,aes(x=labels,y=count,fill=group))+
  geom_bar(position = 'dodge',stat='identity',colour='black')+
  scale_fill_brewer(palette = 'Accent')+
  # facet_zoom(ylim = c(0, 3))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_blank())+
  ggtitle ('overall')

ggsave(filename = '/media/user/sdh/cmq_wj/cellChat/signalingRole_heatmap_top_bar_overall.pdf',
       width = 7,height = 3)

############  merge all heatmap
ihe1<-as.data.frame(iht1@matrix) 
ihe2<-as.data.frame(iht2@matrix)
ihe3<-as.data.frame(iht3@matrix)
ihe4<-as.data.frame(iht4@matrix)







############################
###############
library(ggrepel)
# RR1<-subset(rr,-0.03 < rr$outgoing )
RR1<-rr[which(rr$outgoing < 0.03),]
RR1<-RR1[which(-0.03 < RR1$outgoing),]

RR1<-RR1[which( RR1$incoming<0.03 ),]

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

ggsave('/media/user/sdh/cmq_wj/cellChat/G5_G6_Monocyte_CCL3_signalingChanges_gg.pdf',
       width =8,
       height = 5)

pathways.show<-'CD40'
NH1<-netVisual_heatmap(cc.g5, signaling = pathways.show, color.heatmap = "Reds")
NH2<-netVisual_heatmap(cc.g6, signaling = pathways.show, color.heatmap = "Reds")
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_heatmap_',pathways.show,'_g5_g6.pdf'))
par(mfrow = c(1,2));NH1;NH2;dev.off()
##########
levels(mg.cc@idents$joint)
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_Monocyte_CCL3_ALL.pdf'))
netVisual_bubble(mg.cc, sources.use = 5, targets.use = c(1:4,6:9),  comparison = c(1,2,3,4), angle.x = 45)
dev.off()
netVisual_bubble(mg.cc, sources.use = 5, targets.use = c(2,3),  comparison = c(1, 2), angle.x = 45)


pdf(paste0('/media/user/sdh/cmq_wj/cellChat/netVisual_bubble_1.pdf'))
netVisual_bubble(mg.cc, sources.use = c(5,6,8), targets.use = c(2),  comparison = c(1,2,3,4), angle.x = 45)
dev.off()
