object.list<-xxxx
outgoing_lis<-list()
incoming_lis<-list()

pathway.union <- unique(c(object.list[[i]]@netP$pathways,
                          object.list[[i+1]]@netP$pathways,
                          object.list[[i+2]]@netP$pathways,
                          object.list[[i+3]]@netP$pathways))
for (i in c(1:4)){
  
  object<-object.list[[i]]
  pattern = "outgoing"
  signaling = pathway.union
  title = names(object.list)[i]
  width = 5;height = 7
  color.heatmap = "BuGn"
  slot.name = "netP"
  font.size = 8;font.size.title = 10
  cluster.rows = FALSE;cluster.cols = FALSE

centr <- slot(object, slot.name)$centr
outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(object@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)
for (i in 1:length(centr)) {
  outgoing[,i] <- centr[[i]]$outdeg
  incoming[,i] <- centr[[i]]$indeg
}

mat <- t(outgoing)
if (!is.null(signaling)) {
  mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
  mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
  idx <- match(rownames(mat1), signaling)
  mat[idx[!is.na(idx)], ] <- mat1
  dimnames(mat) <- list(signaling, colnames(mat1))
}
mat.ori <- mat
mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
mat[mat == 0] <- NA

inmat <- t(incoming)
if (!is.null(signaling)) {
  inmat1 <- inmat[rownames(inmat) %in% signaling, , drop = FALSE]
  inmat <- matrix(0, nrow = length(signaling), ncol = ncol(inmat))
  inidx <- match(rownames(inmat1), signaling)
  inmat[inidx[!is.na(inidx)], ] <- inmat1
  dimnames(inmat) <- list(signaling, colnames(inmat1))
}
inmat.ori <- inmat
inmat <- sweep(inmat, 1L, apply(inmat, 1, max), '/', check.margin = FALSE)
inmat[inmat == 0] <- NA

colnames(mat.ori)<-paste(colnames(mat.ori),title,sep = '_')
colnames(inmat.ori)<-paste(colnames(inmat.ori),title,sep = '_')

mat.ori<-list(mat.ori);names(mat.ori)<-title
inmat.ori<-list(inmat.ori);names(inmat.ori)<-title
outgoing_lis<-c(outgoing_lis,mat.ori)
incoming_lis<-c(incoming_lis,inmat.ori)
}

outgoing_all<-cbind(as.data.frame(outgoing_lis[[1]]),
                    as.data.frame(outgoing_lis[[2]]),
                    as.data.frame(outgoing_lis[[3]]),
                    as.data.frame(outgoing_lis[[4]]))

incoming_all<-cbind(as.data.frame(incoming_lis[[1]]),
                    as.data.frame(incoming_lis[[2]]),
                    as.data.frame(incoming_lis[[3]]),
                    as.data.frame(incoming_lis[[4]]))

overall_all<-outgoing_all+incoming_all

outgoing_all_ <- sweep(outgoing_all, 1L, 
                       apply(outgoing_all, 1, max), '/', 
                       check.margin = FALSE)
incoming_all_ <- sweep(incoming_all, 1L,
                       apply(incoming_all, 1, max), '/', 
                       check.margin = FALSE)
outgoing_all <- sweep(outgoing_all, 1L, 
                       apply(outgoing_all, 1, max), '/', 
                       check.margin = FALSE)
incoming_all <- sweep(incoming_all, 1L,
                       apply(incoming_all, 1, max), '/', 
                       check.margin = FALSE)
overall_all <- sweep(overall_all, 1L, 
                       apply(overall_all, 1, max), '/', 
                       check.margin = FALSE)

colnames(outgoing_all)
for( yu in c(1:31)){
outgoing_all$G7_G8[yu]<-sum(outgoing_all[yu,c(19:27)]) - sum(outgoing_all[yu,c(28:36)])
outgoing_all$G6_G8[yu]<-sum(outgoing_all[yu,c(10:18)]) - sum(outgoing_all[yu,c(28:36)])
outgoing_all$G5_G6[yu]<-sum(outgoing_all[yu,c(1:9)]) - sum(outgoing_all[yu,c(10:18)])
outgoing_all$G5_G7[yu]<-sum(outgoing_all[yu,c(1:9)]) - sum(outgoing_all[yu,c(19:27)])

incoming_all$G7_G8[yu]<-sum(incoming_all[yu,c(19:27)]) - sum(incoming_all[yu,c(28:36)])
incoming_all$G6_G8[yu]<-sum(incoming_all[yu,c(10:18)]) - sum(incoming_all[yu,c(28:36)])
incoming_all$G5_G6[yu]<-sum(incoming_all[yu,c(1:9)]) - sum(incoming_all[yu,c(10:18)])
incoming_all$G5_G7[yu]<-sum(incoming_all[yu,c(1:9)]) - sum(incoming_all[yu,c(19:27)])

overall_all$G7_G8[yu]<-sum(overall_all[yu,c(19:27)]) - sum(overall_all[yu,c(28:36)])
overall_all$G6_G8[yu]<-sum(overall_all[yu,c(10:18)]) - sum(overall_all[yu,c(28:36)])
overall_all$G5_G6[yu]<-sum(overall_all[yu,c(1:9)]) - sum(overall_all[yu,c(10:18)])
overall_all$G5_G7[yu]<-sum(overall_all[yu,c(1:9)]) - sum(overall_all[yu,c(19:27)])
}
write.csv(incoming_all,'/media/user/sdh/cmq_wj/cellChat/incoming_all_heatmap.csv')
write.csv(outgoing_all,'/media/user/sdh/cmq_wj/cellChat/outgoing_all_heatmap.csv')
write.csv(overall_all,'/media/user/sdh/cmq_wj/cellChat/overall_all_heatmap.csv')
pheatmap::pheatmap(outgoing_all[,37:40],scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='outgoing',
                   gaps_row = c(8, 16, 24),
                    border=FALSE)

color.heatmap = "OrRd"
pdf('/media/user/sdh/cmq_wj/cellChat/outgoing_all_heatmap.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(outgoing_all_,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='outgoing',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(9, 18, 27), border=FALSE)
dev.off()
pdf('/media/user/sdh/cmq_wj/cellChat/incoming_all_heatmap.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(incoming_all_,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='incoming',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(9, 18, 27), border=FALSE)
dev.off()

pdf('/media/user/sdh/cmq_wj/cellChat/overall_all_heatmap.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(overall_all,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='overall',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(9, 18, 27), border=FALSE)
dev.off()
grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
outgoing_all_<-outgoing_all_[,order(colnames(outgoing_all_))]
incoming_all_<-incoming_all_[,order(colnames(incoming_all_))]
overall_all_<-overall_all[,order(colnames(overall_all))]

pdf('/media/user/sdh/cmq_wj/cellChat/outgoing_all_heatmap_celltype.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(outgoing_all_,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='outgoing',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(4, 8, 12,16,20,24,28,32), border=FALSE)
dev.off()
pdf('/media/user/sdh/cmq_wj/cellChat/incoming_all_heatmap_celltype.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(incoming_all_,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='incoming',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(4, 8, 12,16,20,24,28,32), border=FALSE)
dev.off()

pdf('/media/user/sdh/cmq_wj/cellChat/overall_all_heatmap_celltype.pdf',
    width = 10,height = 6)
pheatmap::pheatmap(overall_all_,scale = 'none',
                   cluster_rows = T,cluster_cols = F,
                   color = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100),
                   main='overall',
                   gaps_row = c(8, 16, 24),
                   gaps_col = c(4, 8, 12,16,20,24,28,32), border=FALSE)
dev.off()

# color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
pheatmap::pheatmap(outgoing_all,scale = 'row',
                   cluster_rows = F,cluster_cols = F)
pheatmap::pheatmap(incoming_all,scale = 'row',
                   cluster_rows = F,cluster_cols = F)

if (!is.null(signaling)) {
  mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
  mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
  idx <- match(rownames(mat1), signaling)
  mat[idx[!is.na(idx)], ] <- mat1
  dimnames(mat) <- list(signaling, colnames(mat1))
}
mat.ori <- mat
mat <- sweep(mat, 1L, apply(mat, 1, max), '/', check.margin = FALSE)
mat[mat == 0] <- NA