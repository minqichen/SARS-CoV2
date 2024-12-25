setwd('/media/user/sdh/cmq_wj/markergene')
library(org.Mm.eg.db)
library(clusterProfiler)
path<-getwd()
dir(path)
aa<-dir(path)[grep(".txt", dir(path))]
aa<-unlist( strsplit(aa, ".txt", fixed= T))

for (i in aa) {
  gene<-read.table(paste0(i,'.txt'))
  gene$thres<-  as.factor(ifelse(gene[,1] < 0.05 & abs(gene[,2]) > 0.58, 
                                 ifelse(gene[,2]> 0.58,'Up','Down'),'NoSignifi'))

   
   for (yy in c('Up','Down')){
     gene1<-subset(gene,gene$thres == yy )
     aa1<-rownames(gene1)

     
     
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
     
     write.csv(n,file = paste0('../markergene_GO/',i,'_',yy,'_GO_thres0.58_0.05.csv'))

   }
 
}

# volcano
for (i in aa) {
   gene<-read.table(paste0(i,'.txt'))
   gene$thres<-  as.factor(ifelse(gene[,1] < 0.05 & abs(gene[,2]) > 0.58, 
                                  ifelse(gene[,2]> 0.58,'Up','Down'),'NoSignifi'))
   
   ggplot(data = gene, aes(x=avg_logFC, y=-(log10 (gene[,5])),colour= thres )) +
      geom_point() +
      scale_color_manual(values=c("blue", "grey","red"))+
      labs(x="avg_logFC",y="-log10 P adj",title=i)+ 
      # annotate("text", x = 10, y = 5, label =length(gene$thres[gene$thres=='Up']) ,
      #          colour= "red", hjust = -0.08, size = 5 )+
      # annotate("text", x = 10, y = -5, label =length(gene$thres[gene$thres=='Down']) ,
      #          colour= "blue", hjust = -0.08, size = 5 )+
      theme_bw()+#去除背景
      theme(panel.grid.minor = element_blank(),panel.grid.major =element_blank())#去除网格线
   ggsave(paste0('../volcano/',i,'_volcano.pdf'),width = 6,height = 5)
}
