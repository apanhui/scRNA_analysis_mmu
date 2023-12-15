library(ggplot2)
tg <- read.table('Fig2d.Mbp.draw4.txt',header=T ,sep="\t",check.names=FALSE, stringsAsFactors=F)

tg$sample <- factor(tg$sample,levels=c("dbm","dbdb"))
tg$type <- factor(tg$type,levels=unique(tg$type))

p0  <- ggplot(tg,aes(x=type,y=exp))+geom_violin(aes(fill = sample)) + 
         scale_fill_manual(name='', values = c("blue", "red")) +
         xlab("") + ylab("Dlg1") + scale_y_continuous(limits=c(0,2.3138)) + 
         geom_text(data = subset(tg, lable!="--"),aes(y = line*1.11, label = lable, group = type))  + 
         geom_text(data = subset(tg, lable!="--"),aes(y = line*1.1, label = "_______", group = type)) + 
         theme_bw() + 
         theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
         #theme(panel.border = element_blank()) + 
         theme(legend.position = 'top', #axis.ticks = element_blank(),
           axis.title.y = element_text(size=18),axis.text=element_text(size=10))

ggsave("Mbp.violin.pdf",p0, width=8,height=5)

