library(pheatmap)
library(RColorBrewer)

data<-read.table("Fig2b.heatmap.xls", sep="\t", header=T, check=F, comment="", na.strings="", stringsAsFactors=F, row=1)

# col annotaion
cells = colnames(data)
colanno1 = rep("Oligodendrocytes",length(cells))
colanno2 = rep("dbm",length(cells))
colanno2[grep("dbdb_",cells)] = "dbdb"
annotation = data.frame(cluster=colanno1,contrast=colanno2)
rownames(annotation) = cells

# set color
cluster_color = c("#925E9F")
names(cluster_color) = c("Oligodendrocytes")
group_color = c("#3352a3","#e92124")
names(group_color) = c("dbm","dbdb")
annotation_colors = list(cluster=cluster_color,contrast=group_color)

color = c("purple", "black", "yellow")
color <- colorRampPalette(color)(100)

mat <- as.matrix(data)
MAX <- max(abs(mat))
breaks <- seq(-1 * MAX, MAX, length.out = length(color))
pheatmap(mat, cluster_cols = F, show_colnames = F,
		annotation_col = annotation, annotation_colors = annotation_colors,
		scale = "none", color = color, breaks = breaks,
		treeheight_row = 10,
		cellheight = 10, #cellwidth = 10,
		filename = "Heatmap.dbm-vs-dbdb.Oligodendrocytes.pdf")
