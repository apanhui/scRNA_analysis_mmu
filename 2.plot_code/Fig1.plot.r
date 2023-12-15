library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

load("obj.Rda")

# b. plot tSNE
p <- DimPlot(object, reduction = 'tsne', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = FALSE)
gsave(p, file = "tSNE.pdf")

# d-e. marker tSNE
p <- FeaturePlot(object, features = features, order = TRUE, reduction = "tsne", combine = TRUE, label = FALSE, cols = c("lightgrey", "#ff0000", "#00ff00"))
ggsave(p, file = "marker.tSNE.pdf")

# f. plot tSNE of each sample
for ( i in unique(object@meta.data$orig.ident ){
	cells.use <- rownames(object@meta.data)[object@meta.data$orig.ident == i]
	p <- DimPlot(object, reduction = 'tsne', group.by = "seurat_clusters", cols=object@misc$color.cluster, label = FALSE, cells=cells.use)
	gsave(p, file = paste0("tSNE.",i,".pdf"))
}

