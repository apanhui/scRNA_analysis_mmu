library(CellChat)
library(reshape)

file1 <- "cellchat.dbm.rds"
file2 <- "cellchat.dbdb.rds"
sample1 <- "dbm"
sample2 <- "dbdb"

samples <- c(sample1, sample2)

object.list <- list()
object.list[[sample1]] <- readRDS(file1)
object.list[[sample2]] <- readRDS(file2)

ident1 <- levels(object.list[[sample1]]@idents)
ident2 <- levels(object.list[[sample2]]@idents)

idents <- union(ident1, ident2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

count1 = cellchat@net[[sample1]]$count
count2 = cellchat@net[[sample2]]$count
dt <- melt(count1, varnames = c('Source', 'Target'), value.name = "Count_dbm") %>% full_join(y = melt(count2, varnames = c('Source', 'Target'), value.name = "Count_dbdb"))
dt[["Count_dbm-vs-db"]] <- dt[[4]] - dt[[3]]

weight1 = cellchat@net[[sample1]]$weight
weight2 = cellchat@net[[sample2]]$weight
dt <- melt(weight1, varnames = c('Source', 'Target'), value.name = "Prob_dbm") %>% full_join(y = melt(weight2, varnames = c('Source', 'Target'), value.name = "Prob_dbdb"))
dt[["Prob_dbm-vs-db"]] <- dt[[7]] - dt[[6]]

diff_number <- dt

# get top5
Source <- as.character(head((diff_number %>% group_by(Source) %>% summarise(count = sum(abs(diff_count))) %>% arrange(-count))$Source, 5))
Target <- as.character(head((diff_number %>% group_by(Target) %>% summarise(count = sum(abs(diff_count))) %>% arrange(-count))$Target, 5))
Signaling <- head((subset(df, source %in% Source & target %in% Target) %>% group_by(pathway_name) %>% summarise(prob = sum(prob.original)) %>% arrange(-prob))$pathway_name, 5)

# plot bubble
color.text = c("#377EB8","#FF7F00")
dtp <- netVisual_bubble(cellchat, comparison = c(1,2), angle.x = 45, remove.isolate = T, color.text = color.text, return.data = TRUE)
p <- dtp[['gg.obj']]
ggsave(p, file = "Diff.CommunProb.top5.bubble.pdf", width = 10, height = 6, limitsize = FALSE)


