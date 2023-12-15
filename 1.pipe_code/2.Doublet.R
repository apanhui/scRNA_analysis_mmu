library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(DoubletFinder)

dir("dbdb")
# [1] "barcodes.tsv" "genes.tsv"
# [2] "matrix.mtx"
### Creat Seurat Object
object.list <- list()
data_name <- c("dbdb")
for (i in seq(data_name)) {
		mat <- Read10X(data.dir = data_name[i], gene.column = 1)
		object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = "RNA")
}

## merge Seurat object
object <- merge(x = object.list[[1]], y = unlist(object.list[-1]),add.cell.ids = data_name)
sample <- "dbdb"
### Normalization Data
object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)

### Reduce dimension
object <- RunPCA(object,assay = DefaultAssay(object), npcs = 50,features = VariableFeatures(object), verbose = FALSE)
object <- RunUMAP(object, dims = 1:50, umap.method = "uwot", n.neighbors = 30, n.components = 2, reduction.name = 'umap', reduction = "pca")
perplexity <- min(30, floor((ncol(object) - 1)/3))
object <- RunTSNE(object, dims = seq(50),dim.embed = 2,reduction.name = 'tsne',perplexity = perplexity,reduction ="pca",check_duplicates = FALSE)

dims <- seq(object@reductions$pca)

### get pN
pN <- 0.25

### get pK
sweep.res.list <- paramSweep_v3(object, PCs = dims, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list)
bcmvn <- find.pK(sweep.stats)
pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)[1]]))

pdf(paste0("pK.", sample, ".pdf"))
plot(x = as.vector(bcmvn$pK), y = bcmvn$BCmetric, type = "b", xlab = "pK", ylab = "BCmetric", col = "blue", pch = 19)
abline(v = pK, lty = 2, col = "red")
dev.off()

### get nExp 
if ( is.na(rate) ) 
		rate <- 7.6 * 10^-6 * ncol(object) + 5.27 * 10^-4
nExp_poi <- round(as.numeric(rate) * ncol(object)) 

if( exists("seurat_clusters", object@meta.data) ){
		annotations <- object@meta.data$seurat_clusters
		homotypic.prop <- modelHomotypic(annotations)
		nExp_poi <- round(nExp_poi * (1 - homotypic.prop))
}

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
object <- doubletFinder_v3(object, PCs = dims, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(object@meta.data)[grep('pANN', colnames(object@meta.data))] <- "pANN"
colnames(object@meta.data)[grep('DF.classifications', colnames(object@meta.data))] <- "classifications"

embeddings <- cbind(object@reductions[["tsne"]]@cell.embeddings, object@reductions[["umap"]]@cell.embeddings)
data <- object@meta.data[,c("pANN", "classifications")]
write.table(cbind(Cells = rownames(data), data), file = paste0("DF.classify.", sample, ".xls"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

p1 <- DimPlot(object, reduction = "umap", group.by = "classifications", cols = c("Singlet" = "black", "Doublet" = "red")) + ggtitle(NULL)
ggsave(p1, file = paste0("DF.classify.UMAP.", sample, ".pdf"), width = 6, height = 5)
p2 <- DimPlot(object, reduction = "tsne", group.by = "classifications", cols = c("Singlet" = "black", "Doublet" = "red"))  + ggtitle(NULL)
ggsave(p2, file = paste0("DF.classify.tSNE.", sample, ".pdf"), width = 6, height = 5)


sessionInfo()
#R version 4.1.3 (2022-03-10)
#Platform: x86_64-conda-linux-gnu (64-bit)
#Running under: CentOS release 6.9 (Final)

#Matrix products: default

#locale:
#[1] LC_CTYPE=en_US.iso885915       LC_NUMERIC=C
#[3] LC_TIME=en_US.iso885915        LC_COLLATE=en_US.iso885915
#[5] LC_MONETARY=en_US.iso885915    LC_MESSAGES=en_US.iso885915
#[7] LC_PAPER=en_US.iso885915       LC_NAME=C
#[9] LC_ADDRESS=C                   LC_TELEPHONE=C
#[11] LC_MEASUREMENT=en_US.iso885915 LC_IDENTIFICATION=C

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base

#other attached packages:
#[1] DoubletFinder_2.0.3 patchwork_1.1.1     ggplot2_3.3.5
#[4] dplyr_1.0.9         SeuratObject_4.0.4  Seurat_4.1.0

#loaded via a namespace (and not attached):
#[1] nlme_3.1-157          matrixStats_0.61.0    spatstat.sparse_3.0-2
#[4] RcppAnnoy_0.0.19      RColorBrewer_1.1-3    httr_1.4.2
#[7] sctransform_0.3.3     tools_4.1.3           utf8_1.2.2
#[10] R6_2.5.1              irlba_2.3.5           rpart_4.1.16
#[13] KernSmooth_2.23-20    uwot_0.1.11           mgcv_1.8-40
#[16] DBI_1.1.3             lazyeval_0.2.2        colorspace_2.0-3
#[19] withr_2.5.0           tidyselect_1.1.2      gridExtra_2.3
#[22] compiler_4.1.3        cli_3.5.0             plotly_4.10.0
#[25] scales_1.2.0          lmtest_0.9-40         spatstat.data_3.0-1
#[28] ggridges_0.5.3        pbapply_1.5-0         goftest_1.2-3
#[31] stringr_1.4.0         digest_0.6.29         spatstat.utils_3.0-3
#[34] pkgconfig_2.0.3       htmltools_0.5.2       parallelly_1.33.0
#[37] fastmap_1.1.0         htmlwidgets_1.5.4     rlang_1.0.6
#[40] shiny_1.7.1           generics_0.1.2        zoo_1.8-10
#[43] jsonlite_1.8.0        ica_1.0-2             spatstat.random_3.1-5
#[46] magrittr_2.0.3        Matrix_1.4-1          Rcpp_1.0.8.3
#[49] munsell_0.5.0         fansi_1.0.3           abind_1.4-5
#[52] reticulate_1.26-9000  lifecycle_1.0.3       yaml_2.3.5
#[55] stringi_1.7.6         MASS_7.3-56           Rtsne_0.16
#[58] plyr_1.8.7            grid_4.1.3            parallel_4.1.3
#[61] listenv_0.8.0         promises_1.2.0.1      ggrepel_0.9.1
#[64] crayon_1.5.1          miniUI_0.1.1.1        deldir_1.0-6
#[67] lattice_0.20-45       cowplot_1.1.1         splines_4.1.3
#[70] tensor_1.5            pillar_1.7.0          igraph_1.3.0
#[73] spatstat.geom_3.2-1   future.apply_1.8.1    reshape2_1.4.4
#[76] codetools_0.2-18      leiden_0.3.9          glue_1.6.2
#[79] data.table_1.14.2     png_0.1-7             vctrs_0.5.1
#[82] httpuv_1.6.5          gtable_0.3.0          RANN_2.6.1
#[85] purrr_0.3.4           spatstat.core_2.4-2   polyclip_1.10-0
#[88] tidyr_1.2.0           scattermore_0.8       future_1.30.0
#[91] assertthat_0.2.1      mime_0.12             xtable_1.8-4
#[94] RSpectra_0.16-0       later_1.2.0           survival_3.3-1
#[97] viridisLite_0.4.0     tibble_3.1.6          cluster_2.1.3
#[100] globals_0.16.2        fitdistrplus_1.1-8    ellipsis_0.3.2
#[103] ROCR_1.0-11
