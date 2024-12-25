setwd('/media/user/sdh/cmq_wj/G_allmarker')
library(org.Mm.eg.db)
library(clusterProfiler)
path<-getwd()
dir(path)
aa<-dir(path)[grep(".CSV", dir(path))]
aa<-unlist( strsplit(aa, ".CSV", fixed= T))

for (i in aa) {
  gene<-read.csv(paste0(i,'.CSV'))
  gene$thres<-  as.factor(ifelse(gene[,2] < 0.05 & abs(gene[,3]) > 0.58, 
                                 ifelse(gene[,3]> 0.58,'Up','Down'),'NoSignifi'))
  
  
  for (yy in c('Up','Down')){
    gene1<-subset(gene,gene$thres == yy )
    aa1<-gene1$X
    
    
    
    gene.df.trans<-bitr(aa1, fromType = "SYMBOL",
                        toType = c("ENTREZID",'ENSEMBL'),
                        OrgDb = org.Mm.eg.db)
    
    ego_bp <- enrichGO(
      gene  = gene.df.trans$ENTREZID,
      keyType = "ENTREZID",
      OrgDb   = org.Mm.eg.db,
      ont     = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.01,
      qvalueCutoff  = 0.05,
      readable      = TRUE)
    
    n = ego_bp@result
    
    write.csv(n,file = paste0('/media/user/sdh/cmq_wj/G_marker_go/',i,'_',yy,'_GO_thres0.58_0.05.csv'))
    
  }
  
}

# volcano
for (i in aa) {
  gene<-read.csv(paste0(i,'.CSV'))
  gene$thres<-  as.factor(ifelse(gene[,2] < 0.05 & abs(gene[,3]) > 0.58, 
                                 ifelse(gene[,3]> 0.58,'Up','Down'),'NoSignifi'))
  
  ggplot(data = gene, aes(x=avg_logFC, y=-(log10 (gene[,6])),colour= thres )) +
    geom_point() +
    scale_color_manual(values=c("blue", "grey","red"))+
    labs(x="avg_logFC",y="-log10 P adj",title=i)+ 
    # annotate("text", x = 10, y = 5, label =length(gene$thres[gene$thres=='Up']) ,
    #          colour= "red", hjust = -0.08, size = 5 )+
    # annotate("text", x = 10, y = -5, label =length(gene$thres[gene$thres=='Down']) ,
    #          colour= "blue", hjust = -0.08, size = 5 )+
    theme_bw()+#去除背景
    theme(panel.grid.minor = element_blank(),panel.grid.major =element_blank())#去除网格线
  ggsave(paste0('/media/user/sdh/cmq_wj/G_marker_go/PLOT/',i,'_volcano.pdf'),width = 6,height = 5)
}


#######  PLOT GO 
setwd('/media/user/sdh/cmq_wj/G_marker_go')
.libPaths("/media/user/sdf/R.lib/library/wj/reg")
library(ggplot2)
library(ggrepel)
.libPaths("/media/user/sdf/R.lib/library/cmq")
library(egg)

for (path in dir(getwd())[grep(".csv", dir(getwd()))]){
  go_term<-read.csv(path)
  go_term<-go_term[1:10,]
  go_term$GeneRatio1<-go_term$Count / as.numeric(strsplit(as.character(go_term$GeneRatio[1]) ,split = '/')[[1]][2]) 
  
  mytheme <- theme(axis.title=element_text(face="italic", size=14,colour = 'black'), #坐标轴标题
                   axis.text=element_text(face="italic", size=14,colour = 'black'), #坐标轴标签
                   axis.line = element_line(size=0.5, colour = 'black'), #轴线
                   panel.background = element_rect(color='black'), #绘图区边框
                   legend.key = element_blank() #关闭图例边框
  )
  
  
  p<-ggplot(go_term,aes(x=-1*log10(pvalue),y=reorder(go_term[,3],-go_term[,6]),
                        colour=GeneRatio1,size=Count))+
    geom_point()+
    scale_size(range=c(2, 8))+
    scale_colour_gradient(low = "blue",high = "red")+
    theme_bw()+
    ylab("GO Pathway Terms")+
    xlab("-log10 Pvalue")+
    
    labs(color=expression(GeneRatio))+
    theme(legend.title=element_text(size=14), legend.text = element_text(size=14))+
    theme(axis.title.y = element_text(margin = margin(r = 50)),
          axis.title.x = element_text(margin = margin(t = 20)))+
    theme(axis.text.x = element_text(face ="italic",color="black",angle=0,vjust=1))+
    ggtitle(strsplit(path,split = '_GO_thres0.58_0.05.csv')) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot <- p+mytheme  
  plot
  
  ggsave(paste0('/media/user/sdh/cmq_wj/G_marker_go/PLOT/',path,'.pdf'), 
         egg::set_panel_size(plot, width=unit(3.5, "in"), height=unit(5, "in")), 
         width = 20, height = 7, units = 'in', dpi = 300)
}
