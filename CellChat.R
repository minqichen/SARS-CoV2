.libPaths("/media/user/sde/xiejiayi_database/R_py_script/R/x86_64-pc-linux-gnu-library/3.5")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(future)
.libPaths("/media/user/sdf/R.lib/library/wj/sc")
library(Seurat)

sc2 <- readRDS("/media/user/sdh/wj_kx/analysis/RDS/lung_immune_cell_annotated.rds")
ct_all_seu_rn_nwid<-subset(sc2,cells = names(sc2$orig.ident[sc2$orig.ident == "G8"]))
data.input<-ct_all_seu_rn_nwid@assays$RNA@data

# 1. create cellchat object
data.input1<-ct_all_seu_rn_nwid@assays$RNA@data
identity<-as.data.frame(ct_all_seu_rn_nwid@active.ident)
colnames(identity)<-"labels"
cellchat<-createCellChat(data.input1)
cellchat<-addMeta(cellchat,meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize = as.numeric(table(cellchat@idents))

# 2. input receptor-ligand database
CellChatDB <- CellChatDB.mouse
#unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#ligand<-unique(CellChatDB.use[["interaction"]][["ligand"]])
#receptor<-unique(CellChatDB.use[["interaction"]][["receptor"]])
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat@DB <- CellChatDB.use

# 3. preprocessing
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat,PPI.mouse)

# 4. compute interaction
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)




cc.g5 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G5.rds")
cc.g6 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G6.rds")
cc.g7 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G7.rds")
cc.g8 <- readRDS("/media/user/sdh/wj_kx/analysis/CellChat/cellchat.G8.rds")



netVisual_circle(cellchat@net$count,vertex.weight = groupSize, 
                 weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight,vertex.weight = groupSize, 
                 weight.scale = T, label.edge = F, title.name = "Weight of interactions")
cellchat@netP$pathways
pathway.show = "CCL"
pathway.show = "TNF"
pathway.show = "IL16"
pathway.show = "IFN-II"
pathway.show = "CXCL"
pathway.show = "TGFb"
pathway.show = "IL2"
pathway.show = "IL1"
pathway.show = "SEMA3"
netVisual_heatmap(cellchat,signaling = pathway.show,color.heatmap = "Reds")

# compare two groups

netAnalysis_dot(cellchat, pattern = "outgoing")
mynetVisual_aggregate(cellchat, signaling = c("MIF"), layout = "circle", 
                      vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7)

df.net <- subsetCommunication(cellchat,slot.name = "net")
df.netP <- subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.net,"/media/user/sdd/cellchat_secreted.csv")

df.dot<-read.csv("/media/user/sdd/cellchat_secreted_select.csv")
for (i in 1:nrow(df.dot)){
  df.dot[i,8]<-paste0(df.dot[i,1],"_",df.dot[i,2])
}
colnames(df.dot)[8]<-"pair"

dotp<-ggplot(df.dot,aes(x=pair,y=interaction_name_2))+
  geom_point(aes(size=pval,colour=(log2(prob))))+
  scale_size(limits = c(1,2))+
  scale_color_gradientn(colors = c("#F43B3F","#F3E544","#3EBEC3","#1E3164")[4:1])
ggsave("/media/user/sdd/cellchat/cellchat_dot.pdf",
       dotp,width = 8,height = 4)

# CellChat: compare two group: G7 vs G8
## merge G5, G6, G7 and G8

xxxx <- list(g5 = cc.g5, g6 = cc.g6, g7 = cc.g7, g8 = cc.g8)
mg.cc <- xxxx

mg.cc
mg.cc <- mergeCellChat(mg.cc,add.names = names(mg.cc),cell.prefix = T)

mg.cc

## general visualization

gg1 <- compareInteractions(mg.cc,show.legend = F, group = c(1,2,3,4), measure = "count")
gg1
gg2 <- compareInteractions(mg.cc,show.legend = F, group = c(1,2,3,4), measure = "weight")
gg2

par(mfrow = c(1,2))
netVisual_diffInteraction(mg.cc, weight.scale = T)
netVisual_diffInteraction(mg.cc, weight.scale = T, measure = "weight")

par(mfrow = c(1,2))
h1 <- netVisual_heatmap(mg.cc)
h1
h2 <- netVisual_heatmap(mg.cc, measure = "weight")
h2
h1 + h2

# general number of interaction
par(mfrow = c(2,2))
cnt1 <- mg.cc@net[[1]]$count
cnt2 <- mg.cc@net[[2]]$count
cnt3 <- mg.cc@net[[3]]$count
cnt4 <- mg.cc@net[[4]]$count

weight.max <- max(max(cnt1),max(cnt2),max(cnt3),max(cnt4))

netVisual_circle(cnt1, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Number of interaction: ",names(mg.cc@net)[1]))
netVisual_circle(cnt2, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Number of interaction: ",names(mg.cc@net)[2]))
netVisual_circle(cnt3, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Number of interaction: ",names(mg.cc@net)[3]))
netVisual_circle(cnt4, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Number of interaction: ",names(mg.cc@net)[4]))

# general weight of interaction
par(mfrow = c(2,2))
cnt1 <- mg.cc@net[[1]]$weight
cnt2 <- mg.cc@net[[2]]$weight
cnt3 <- mg.cc@net[[3]]$weight
cnt4 <- mg.cc@net[[4]]$weight

weight.max <- max(max(cnt1),max(cnt2),max(cnt3),max(cnt4))

netVisual_circle(cnt1, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Weight of interaction: ",names(mg.cc@net)[1]))
netVisual_circle(cnt2, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Weight of interaction: ",names(mg.cc@net)[2]))
netVisual_circle(cnt3, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Weight of interaction: ",names(mg.cc@net)[3]))
netVisual_circle(cnt4, weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.05,
                 edge.width.max = 10, title.name = paste0("Weight of interaction: ",names(mg.cc@net)[4]))

# split
cnt1 <- mg.cc@net[[1]]$weight
cnt2 <- mg.cc@net[[2]]$weight
cnt3 <- mg.cc@net[[3]]$weight
cnt4 <- mg.cc@net[[4]]$weight

cnt1 <- mg.cc@net[[1]]$count
cnt2 <- mg.cc@net[[2]]$count
cnt3 <- mg.cc@net[[3]]$count
cnt4 <- mg.cc@net[[4]]$count

split.circle(i = 7)

split.circle <- function(i){
  mat1 <- matrix(0,nrow = nrow(cnt1),ncol = ncol(cnt1),dimnames = dimnames(cnt1))
  mat2 <- matrix(0,nrow = nrow(cnt2),ncol = ncol(cnt2),dimnames = dimnames(cnt2))
  mat3 <- matrix(0,nrow = nrow(cnt3),ncol = ncol(cnt3),dimnames = dimnames(cnt3))
  mat4 <- matrix(0,nrow = nrow(cnt4),ncol = ncol(cnt4),dimnames = dimnames(cnt4))
  
  #mat1[i,] <- cnt1[i,]
  #mat2[i,] <- cnt2[i,]
  #mat3[i,] <- cnt3[i,]
  #mat4[i,] <- cnt4[i,]
  
  mat1[,i] <- cnt1[,i]
  mat2[,i] <- cnt2[,i]
  mat3[,i] <- cnt3[,i]
  mat4[,i] <- cnt4[,i]
  
  weight.max <- max(mat1,mat2,mat3,mat4)
  par(mfrow = c(2,2))
  netVisual_circle(mat1,weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.1,
                   title.name = paste0("Weight of interaction: ",names(mg.cc@net)[1]))
  netVisual_circle(mat2,weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.1,
                   title.name = paste0("Weight of interaction: ",names(mg.cc@net)[2]))
  netVisual_circle(mat3,weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.1,
                   title.name = paste0("Weight of interaction: ",names(mg.cc@net)[3]))
  netVisual_circle(mat4,weight.scale = T, label.edge = T, edge.weight.max = weight.max, arrow.size = 0.1,
                   title.name = paste0("Weight of interaction: ",names(mg.cc@net)[4]))
}



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

## conserved and specific pathway
pp1 <- rankNet(mg.cc, mode = "comparison", stacked = T, do.stat = T)
pp1
pp2 <- rankNet(mg.cc, mode = "comparison", comparison = c(2,4),stacked = F, do.stat = T)
pp2



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
