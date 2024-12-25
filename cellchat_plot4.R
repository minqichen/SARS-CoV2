pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G5-8_',pathways.show,'_netVisual_heatmap.pdf'))
netVisual_heatmap(cc.g5, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(cc.g6, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(cc.g7, signaling = pathways.show, color.heatmap = "Reds")
netVisual_heatmap(cc.g8, signaling = pathways.show, color.heatmap = "Reds")
dev.off()



pairLR.CXCL <- extractEnrichedLR(cc.g5, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cc.g5, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cc.g5, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cc.g5, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G5-8_',LR.show,'_circle.pdf'),width=5,height=5)
netVisual_individual(cc.g5, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")       
netVisual_individual(cc.g6, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")   
netVisual_individual(cc.g7, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")   
netVisual_individual(cc.g8, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")   
dev.off()
pdf(paste0('/media/user/sdh/cmq_wj/cellChat/G5-8_',LR.show,'_chord.pdf'),width=8,height=8)
netVisual_individual(cc.g5, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")       
netVisual_individual(cc.g6, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")   
netVisual_individual(cc.g7, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")   
netVisual_individual(cc.g8, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")   
dev.off()

cellchat<-cc.g5
groupSize <- as.numeric(table(cellchat@idents))
pdf('/media/user/sdh/cmq_wj/cellChat/G5-8_chord.pdf')
# par(mfrow = c(1,2), xpd=TRUE)
AA<-netVisual_circle(cc.g5@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.g5@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.g6@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.g6@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.g7@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.g7@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.g8@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.g8@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


pathways.show <- c("CXCL") 
