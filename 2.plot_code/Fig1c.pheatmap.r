library(pheatmap)

data.plot<-read.table("Fig1c.heatmap.xls",sep="\t",header=T,check=F,row=1)

mycolors = colorRampPalette(colors = c("#043062","#2167ab","white","#af152a","#680322"))(100)

p <- pheatmap(data.plot, scale="column", cluster_rows = F,cluster_cols = F,
		 color = mycolors,
		 cellwidth = 20,cellheight = 20,display_numbers = F,fontsize=10,
		 show_rownames = T,show_colnames = T,filename = "heatmap.pdf")

p <- pheatmap(data.plot, scale="column", cluster_rows = F,cluster_cols = F,
		 color = mycolors,
		 cellwidth = 20,cellheight = 20,display_numbers = F,fontsize=10,
		 show_rownames = T,show_colnames = T,filename = "heatmap.png")

