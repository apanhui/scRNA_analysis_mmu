library(ggplot2)

Data<-read.table("Fig2a.dbm-vs-dbdb.diff.stat.xls", head=T,sep = '\t')
colors <- c('#f0b13f','#ed3d47','#ffffff')
p <- ggplot(Data,aes(px,count,fill=color,group=type)) + geom_bar(stat = 'identity',position = 'identity',width=0.6,colour ='White') + 
scale_fill_manual(values=colors) + 
theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + theme(panel.border = element_blank()) +
theme(axis.ticks = element_blank(), axis.text = element_blank(),axis.title=element_blank()) + 
scale_y_continuous(limits=c(-2600,2600),expand = c(0, 0)) + scale_x_continuous(limits=c(1,10.5),expand = c(0, 0)) +
annotate('segment',x= 2,xend = 9, y= 500, yend=500) + annotate('segment',x=2,xend = 9, y= -500, yend=-500) + 
annotate('segment',x= 9, xend = 9, y= 500, yend=2500) + annotate('segment',x= 9, xend = 9, y= -2500, yend=-500) +
annotate('text',x = 9.8, y = 0, label = 'dbm-vs-dbdb',size=6) + 
annotate('text',x = 9.8, y = 1400, label = 'Induced genes',size = 5) + annotate('text',x = 9.8, y = -1400, label = 'Reduced genes', size = 5) +
annotate('text',x = 2.5, y = 0, label = 'Neurons',size = 4) + annotate('text',x = 3.5, y = 0, label = 'Oligodendrocytes',size = 4) + annotate('text',x = 4.5, y = 0, label = 'OPCs',size = 4) + annotate('text',x = 5.5, y = 0, label = 'Microglia',size = 4) + annotate('text',x = 6.5, y = 0, label = 'Fibroblasts',size = 4) + annotate('text',x = 7.5, y = 0, label = 'Endothelial_cells',size = 4) + annotate('text',x = 8.5, y = 0, label = 'Astrocytes',size = 4) +  
#### 
annotate('segment',x=9,xend=9.15,y=500,yend=500) + annotate('segment',x=9,xend=9.15,y=1000,yend=1000) + annotate('segment',x=9,xend=9.15,y=1500,yend=1500) + annotate('segment',x=9,xend=9.15,y=2000,yend=2000) + annotate('segment',x=9,xend=9.15,y=2500,yend=2500) + 
annotate('segment',x=9,xend=9.15,y=-500,yend=-500) + annotate('segment',x=9,xend=9.15,y=-1000,yend=-1000) + annotate('segment',x=9,xend=9.15,y=-1500,yend=-1500) + annotate('segment',x=9,xend=9.15,y=-2000,yend=-2000) + annotate('segment',x=9,xend=9.15,y=-2500,yend=-2500) + 
annotate('text',x=9.3, y = 500, label = '0',size = 3) + annotate('text',x=9.3, y = 1000, label = '500',size = 3) + annotate('text',x=9.3, y = 1500, label = '1000',size = 3) + annotate('text',x=9.3, y = 2000, label = '1500',size = 3) + annotate('text',x=9.3, y = 2500, label = '2000',size = 3) + 
annotate('text',x=9.3, y = -500, label = '0',size = 3) + annotate('text',x=9.3, y = -1000, label = '500',size = 3) + annotate('text',x=9.3, y = -1500, label = '1000',size = 3) + annotate('text',x=9.3, y = -2000, label = '1500',size = 3) + annotate('text',x=9.3, y = -2500, label = '2000',size = 3) + 
#### 
theme(legend.position = 'none') + coord_flip()
ggsave(paste(args[2],'.pdf',sep=''),width=11,height=8)

