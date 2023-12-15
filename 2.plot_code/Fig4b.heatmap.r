library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)


file = "Fig4b.dbm.heatmap.xls"

data<-read.table(file, sep="\t", header=T, check=F, comment="", na.strings="", stringsAsFactors=F)

# set cluster order
cluster.order <- c("Neurons","Oligodendrocytes","OPCs","Microglia","Fibroblasts","Endothelial_cells","Astrocytes")

# set color
use.colors = c("#1d2f84","#69a3c8","#FFFFFF","#ffc77e","#c51c21")

# plot
data.new <- melt(data)
print(head(data.new))
colnames(data.new) <- c("target","source","count")
data.new$target <- factor(data.new$target, levels=rev(cluster.order))
data.new$source <- factor(data.new$source, levels=cluster.order)

p <- ggplot(data.new, aes(x=source, y=target, fill=count)) +
     geom_raster() + 
	 scale_fill_gradientn(colors=use.colors, limits=c(0,45),breaks = seq(0,40,10)) +
	 labs(x="",y="",fill="Number of pairs",title="") + 
	 theme_classic() + 
	 theme(axis.ticks = element_blank(),
		   axis.text.x = element_text(angle=45,hjust=1, vjust=1,size=10),
	       axis.text.y = element_text(size=10))

ggsave(p, file="dbm.heatmap.ggplot.pdf",width=8, height=6)
