.libPaths("/media/user/sde/xiejiayi_database/R_py_script/R/x86_64-pc-linux-gnu-library/3.5")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(future)
.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)
library(ggplot2)


cc.g5 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G5.rds")
cc.g6 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G6.rds")
cc.g7 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G7.rds")
cc.g8 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G8.rds")


# CellChat: compare two group: G7 vs G8
## merge G5, G6, G7 and G8

cc.g5 <- netAnalysis_computeCentrality(cc.g5,slot.name = "netP")
cc.g6 <- netAnalysis_computeCentrality(cc.g6,slot.name = "netP")
cc.g7 <- netAnalysis_computeCentrality(cc.g7,slot.name = "netP")
cc.g8 <- netAnalysis_computeCentrality(cc.g8,slot.name = "netP")
xxxx <- list(g5 = cc.g5, g6 = cc.g6, g7 = cc.g7, g8 = cc.g8)
mg.cc <- xxxx
object.list<-xxxx

mg.cc
mg.cc <- mergeCellChat(mg.cc,add.names = names(mg.cc),cell.prefix = T)




## general visualization

gg1 <- compareInteractions(mg.cc,show.legend = F, group = c(1,2), measure = "count")
gg1
gg2 <- compareInteractions(mg.cc,show.legend = F, group = c(1,2), measure = "weight")
gg2

pdf(paste0('/media/user/sdh/cmq_wj/cellChat/',name,'_netVisual_diffInteraction.pdf'),width = 6,height = 6)
netVisual_diffInteraction(mg.cc, weight.scale = T,comparison = comparison,
                          title.name =paste(name,'Differential number of interactions',sep=' ') )
netVisual_diffInteraction(mg.cc, weight.scale = T, measure = "weight",comparison = comparison,
                          title.name =paste(name,'Differential interaction strength',sep=' ') )
dev.off()

par(mfrow = c(1,2))
h1 <- netVisual_heatmap(mg.cc,title.name =paste(name,'Differential number of interactions',sep=' '),comparison = comparison )
h1
h2 <- netVisual_heatmap(mg.cc, measure = "weight",
                        title.name =paste(name,'Differential interaction strength',sep=' ') ,comparison = comparison)
h2
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/',name,'_netVisual_heatmap.pdf'),width = 12,height = 6)
h1 + h2
dev.off()


## conserved and specific pathway
pp1 <- rankNet(mg.cc, mode = "comparison", comparison = c(3,4), stacked = T, do.stat = T)
pp1
pp2 <- rankNet(mg.cc, mode = "comparison", comparison = c(2,4),stacked = F, do.stat = T)
pp2

pdf('/media/user/sdh/cmq_wj/cellChat/rankNet.pdf',width = 5,height = 5)
rankNet(mg.cc, mode = "comparison", comparison = c(1,2), stacked = T, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(3,4), stacked = T, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(1,3), stacked = T, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(2,4), stacked = T, do.stat = T)
dev.off()

pdf('/media/user/sdh/cmq_wj/cellChat/rankNet_stacked_F.pdf',width = 5,height = 5)
rankNet(mg.cc, mode = "comparison", comparison = c(1,2), stacked = F, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(3,4), stacked = F, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(1,3), stacked = F, do.stat = T)
rankNet(mg.cc, mode = "comparison", comparison = c(2,4), stacked = F, do.stat = T)
dev.off()
###################################
## indicated cell type visualization
par(mfrow = c(1,2))
s.cell <- c("Macrophage_alveolar","Monocyte_Ly6c_low","Monocyte_CCL3","Neutrophil")
cnt1 <- mg.cc@net[[1]]$count[s.cell,s.cell]
cnt2 <- mg.cc@net[[2]]$count[s.cell,s.cell]
weight.max <- max(max(cnt1),max(cnt2))
netVisual_circle(cnt1, weight.scale = T, label.edge = T, edge.weight.max = weight.max,
                 edge.width.max = 12, title.name = "Number of interaction: G7")

netVisual_circle(cnt2, weight.scale = T, label.edge = T, edge.weight.max = weight.max,
                 edge.width.max = 12, title.name = "Number of interaction: G8")


netVisual_diffInteraction(mg.cc, weight.scale = T,
                          measure = "count.merged", label.edge = T)

## conserved and specific pathway



pathway.show <- "IL16"
pathway.show <- "CD40"
pathway.show <- "LIGHT"
pathway.show <- "GALECTIN"
pathway.show <- "PARs"
pathway.show <- "SEMA3"
pathway.show <- "IL1"
pathway.show <- "IFN-II"
pathway.show <- "OSM"
pathway.show <- "ANNEXIN"
pathway.show <- "IGF"
pathway.show <- "GAS"
pathway.show <- "SEMA3"
pathway.show <- "IL2"

pathway.show <- 'Lgals9'
mx1 <- getMaxWeight(xxxx[1],slot.name = "netP", attribute = pathway.show)
mx2 <- getMaxWeight(xxxx[2],slot.name = "netP", attribute = pathway.show)
mx3 <- getMaxWeight(xxxx[3],slot.name = "netP", attribute = pathway.show)
mx4 <- getMaxWeight(xxxx[4],slot.name = "netP", attribute = pathway.show)
mx <- max(mx1,mx2,mx3,mx4)
mx <- max(mx1,mx2)
mx <- mx1
mx <- mx2
par(mfrow = c(2,2),xpd = T)
spec.pp1 <- netVisual_aggregate(xxxx[[1]],signaling = pathway.show, 
                                layout = "circle",edge.weight.max = as.vector(mx),
                                edge.width.max = 10, 
                                signaling.name = paste0(pathway.show," in ",names(xxxx[1])))
spec.pp2 <- netVisual_aggregate(xxxx[[2]],signaling = pathway.show, 
                                layout = "circle",edge.weight.max = as.vector(mx),
                                edge.width.max = 10, 
                                signaling.name = paste0(pathway.show," in ",names(xxxx[2])))
spec.pp3 <- netVisual_aggregate(xxxx[[3]],signaling = pathway.show, 
                                layout = "circle",edge.weight.max = as.vector(mx),
                                edge.width.max = 10, 
                                signaling.name = paste0(pathway.show," in ",names(xxxx[3])))
spec.pp4 <- netVisual_aggregate(xxxx[[4]],signaling = pathway.show, 
                                layout = "circle",edge.weight.max = as.vector(mx),
                                edge.width.max = 10, 
                                signaling.name = paste0(pathway.show," in ",names(xxxx[4])))
## compare receptor-ligand pairs by bubbles
levels(mg.cc@idents$joint)
bub <- netVisual_bubble(mg.cc, sources.use = c(1,3,4,5,6,7,8,9), targets.use = c(2),comparison = c(1,2),angle.x = 45) 
bub <- netVisual_bubble(mg.cc, sources.use = c(2,3,4,5,6,7,8,9), targets.use = c(1),comparison = c(1,2),angle.x = 45) 
bub <- netVisual_bubble(mg.cc, sources.use = c(1,2,3,4,5,6,7,8), targets.use = c(9),comparison = c(1,2),angle.x = 45) 


bub <- netVisual_bubble(mg.cc, sources.use = c(5), targets.use = c(1,2,3,4,6,7,8,9),comparison = c(1,2,3,4),angle.x = 45) 
bub
