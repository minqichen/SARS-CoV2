i=i+1
{
  path<-dir('/media/user/sdh/cmq_wj/markergene_GO')[i]
go_term<-read.csv(path)
go_term$Description[1:20]
}

go_term<-go_term[1:10,]



{
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

ggsave(paste0('../go_plot/',path,'.pdf'), 
       egg::set_panel_size(plot, width=unit(3.5, "in"), height=unit(5, "in")), 
       width = 20, height = 7, units = 'in', dpi = 300)


}