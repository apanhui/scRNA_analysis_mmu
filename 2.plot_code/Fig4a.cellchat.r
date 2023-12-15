reticulate::use_python("/usr/bin/python3", required = T)
options(stringsAsFactors = FALSE)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(CellChat)
library(ggalluvial)
library(svglite)
library(SeuratData)
library(NMF)
library(ComplexHeatmap)
library(yaml)
library(showtext)
showtext_auto(enable = TRUE)

## load obj
obj <- load("obj.Rda")
# get sample
obj <- subset(obj, orig.ident=="dbm")

## Creat cellchat object
cellchat.list <- CreatCellChatList(obj, group.by = "seurat_clusters")
# db
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB 

cellchat <- identifyOverExpressedInteractions(cellchat, features = rownames(cellchat@data.signaling))
cellchat <- computeCommunProb(cellchat, population.size = TRUE)
cellchat <- aggregateNet(cellchat, thresh = 0.05)
saveRDS(cellchat, file = 'cellchat.dbm.Rds')

# plot net
groupSize = as.numeric(table(cellchat@idents))
pdf("LR.interaction.count.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "db/m")
dev.off()
pdf("LR.interaction.weight.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "db/m")
dev.off()

