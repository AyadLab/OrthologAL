library(ggpubr)
library(reshape2)
library(Seurat)
library(msigdbr)
library(tidyverse)
library(cowplot)
library(data.table)
library(dplyr)
library(readr)

#obj_sub is pre conversion mouse data
obj_sub <- subset(obj,status == "mouse")
obj_sub <- RunPCA(
  obj_sub,
  assay = "Spatial",
  rev.pca = TRUE,
  reduction.name = "pca.rev",
  reduction.key = "PCr_", npcs = 50,
  verbose = FALSE
)
#not_scaled_tmp_sub is post converted data
not_scaled_tmp_sub <- subset(not_scaled_tmp,status == "mouse")
not_scaled_tmp_sub <- NormalizeData(not_scaled_tmp_sub,assay = "Spatial_ortho")
not_scaled_tmp_sub <- ScaleData(not_scaled_tmp_sub,assay = "Spatial_ortho")
not_scaled_tmp_sub <- RunPCA(
  not_scaled_tmp_sub,
  assay = "Spatial_ortho",
  rev.pca = TRUE,
  reduction.name = "pca.rev",
  reduction.key = "PCr_", npcs = 50,
  verbose = FALSE
)

E <- obj_sub@reductions$pca.rev@feature.loadings
E_post <- not_scaled_tmp_sub@reductions$pca.rev@feature.loadings
gesecaRes <- geseca(pathways_mouse, E, minSize = 1, maxSize = 50000, center = FALSE, eps=1e-100)
head(gesecaRes)
topPathways <- gesecaRes[, pathway]
gesecaRes$pathway

#length(topPathways) #6010
common_pathways <- intersect(topPathways,topPathways_post)
#length(topPathways_post) #6066
gesecaRes_post <- geseca(pathways,E_post, minSize = 1, maxSize = 50000, center = FALSE, eps=1e-100)
head(gesecaRes_post)
topPathways_post <- gesecaRes_post[, pathway]
gesecaRes_post$pathway

score.geseca <- function(object = patient_tumor, pathways, assay = DefaultAssay(object), prefix = "", scale = FALSE){
  dat <- GetAssay(object, assay)
  expr <- dat$scale.data
  res <- object
  for (i in seq_along(pathways)) {
    genes <- pathways[[i]]
    overlap <- intersect (unique(genes), rownames(expr))
    score <- colSums(expr[overlap, , drop = FALSE])/sqrt(length(overlap))
    scaled <- scale (score, center = TRUE, scale = scale)
    res <- AddMetaData(res, metadata = score, col.name = paste0("raw_", names(pathways)[i]))
    res <- AddMetaData(res, metadata = scaled, col.name = paste0(prefix, names(pathways)[i]))
    #res@meta.data[[paste0((prefix, names(pathways)[i]))]] <- score #BAD
  }
  return(res)
}

obj_sub <- score.geseca(pathways = pathways_mouse[common_pathways], object = obj_sub, assay = "Spatial")
not_scaled_tmp_sub <- score.geseca(pathways = pathways[common_pathways], object = not_scaled_tmp_sub, assay = "Spatial_ortho")

# pulling out only raw values
metacoldf <- data.frame(cols =colnames(obj_sub@meta.data))
for (row in 1:length(rownames(metacoldf))){
  metacoldf$first[row] <- unlist(strsplit(metacoldf$cols[row], split = "_", fixed = T))[1]
}
unique(metacoldf$first)
metacoldf <- subset(metacoldf, metacoldf$first == "raw")

# coercing raw value to matrix
raw_matrix <- as.matrix(t(obj_sub@meta.data[,which(colnames(obj_sub@meta.data) %in% metacoldf$cols)]))
raw_geseca <- CreateAssayObject(data = raw_matrix)
rownames(raw_geseca)
obj_sub[["raw_GESECA"]] <- raw_geseca
DefaultAssay(obj_sub) <- "raw_GESECA"
obj_sub <- ScaleData(obj_sub, assay = "raw_GESECA", features = rownames(raw_geseca))
variablePathways <- gesecaRes[which(gesecaRes$padj < 0.05),]$pathway
obj_sub <- FindVariableFeatures(obj_sub,assay = "raw_GESECA",features = rownames(raw_geseca))
obj_sub <- RunPCA(obj_sub, assay = "raw_GESECA", reduction.name = "pca_geseca",features = variablePathways)

###post converted data

metacoldf1 <- data.frame(cols =colnames(not_scaled_tmp_sub@meta.data))
for (row in 1:length(rownames(metacoldf1))){
  metacoldf1$first[row] <- unlist(strsplit(metacoldf1$cols[row], split = "_", fixed = T))[1]
}
unique(metacoldf1$first)
metacoldf1 <- subset(metacoldf1, metacoldf1$first == "raw")

# coercing raw value to matrix
raw_matrix1 <- as.matrix(t(not_scaled_tmp_sub@meta.data[,which(colnames(not_scaled_tmp_sub@meta.data) %in% metacoldf1$cols)]))
raw_geseca1 <- CreateAssayObject(data = raw_matrix1)
rownames(raw_geseca1)
not_scaled_tmp_sub[["raw_GESECA"]] <- raw_geseca1
DefaultAssay(not_scaled_tmp_sub) <- "raw_GESECA"
not_scaled_tmp_sub <- ScaleData(not_scaled_tmp_sub, assay = "raw_GESECA",features = rownames(raw_geseca1))

variablePathways_post <- gesecaRes_post[which(gesecaRes_post$padj < 0.05),]$pathway
# length(gesecaRes[which(gesecaRes$padj < 0.05),]$pathway)
not_scaled_tmp_sub <- RunPCA(not_scaled_tmp_sub, assay = "raw_GESECA", reduction.name = "pca_geseca", features = variablePathways_post)

r1 <- obj_sub@assays$raw_GESECA$scale.data

r2 <- not_scaled_tmp_sub@assays$raw_GESECA$scale.data

# length(rownames(r2))
# length(rownames(r1))
common_cells <- intersect(rownames(r1), rownames(r2))
r1 <- r1[common_cells, , drop = FALSE]
r2 <- r2[common_cells, , drop = FALSE]

fisher_z <- function(r) {
  r <- pmax(pmin(r, 0.999), -0.999)
  return(0.5 * log((1 + r) / (1 - r)))
}

r1_zcor <- fisher_z(r1)
r1_zcor <- as.data.frame(t(r1_zcor))

r2_zcor <- fisher_z(r2)
r2_zcor <- as.data.frame(t(r2_zcor))
#How to get z_corr value per spot ? subset the healthy tissue aka cerebellum ?
meta <- obj_sub@meta.data
meta_post <- not_scaled_tmp_sub@meta.data
meta$z_cor <- r1_zcor
meta_post$z_cor <- r2_zcor
meta$GOBP_SYNAPTIC_SIGNALING <- meta$GOBP_SYNAPTIC_SIGNALING
common_cells <- intersect(rownames(meta), rownames(meta_post))
df_combined <- data.frame(
  z_cor_pre = meta[common_cells, "z_cor"],
  z_cor_post = meta_post[common_cells, "z_cor"]
)
meta$z_cor <- as.numeric(r1_zcor)
GOBP_SYNAPTIC_SIGNALING_pre <- meta$GOBP_SYNAPTIC_SIGNALING
GOBP_SYNAPTIC_SIGNALING_post <- meta_post$GOBP_SYNAPTIC_SIGNALING
df <- data.frame(GOBP_SYNAPTIC_SIGNALING_pre,GOBP_SYNAPTIC_SIGNALING_post)
df$GOBP_SYNAPTIC_SIGNALING_pre <- as.numeric(df$GOBP_SYNAPTIC_SIGNALING_pre)
df$GOBP_SYNAPTIC_SIGNALING_post <- as.numeric(df$GOBP_SYNAPTIC_SIGNALING_post)

pathway1 <- ggscatter(df, x = "GOBP_SYNAPTIC_SIGNALING_pre", y = "GOBP_SYNAPTIC_SIGNALING_post",
          cor.method = "spearman",
          add = "reg.line", add.params = list(color = "red")) +
  stat_cor(method = "spearman", label.x = min(df$GOBP_SYNAPTIC_SIGNALING_pre),label.y = max(df$GOBP_SYNAPTIC_SIGNALING_post)) +
  geom_smooth(method = "lm", color = "red", fill = "orange", alpha = 0.8)


GOBP_GENERATION_OF_NEURONS_pre <- meta$GOBP_GENERATION_OF_NEURONS
GOBP_GENERATION_OF_NEURONS_post <- meta_post$GOBP_GENERATION_OF_NEURONS

df3 <- data.frame(GOBP_GENERATION_OF_NEURONS_pre,GOBP_GENERATION_OF_NEURONS_post)
df3$GOBP_GENERATION_OF_NEURONS_pre <- as.numeric(df3$GOBP_GENERATION_OF_NEURONS_pre)
df3$GOBP_GENERATION_OF_NEURONS_post <- as.numeric(df3$GOBP_GENERATION_OF_NEURONS_post )

pathway3 <- ggscatter(df3, x = "GOBP_GENERATION_OF_NEURONS_pre", y = "GOBP_GENERATION_OF_NEURONS_post",
                      cor.method = "spearman",
                      add = "reg.line", add.params = list(color = "red")) +
  stat_cor(method = "spearman", label.x = min(df3$GOBP_GENERATION_OF_NEURONS_pre),label.y = max(df3$GOBP_GENERATION_OF_NEURONS_post)) +
  geom_smooth(method = "lm", color = "red", fill = "orange", alpha = 0.8)

#GOBP_COLLAGEN_FIBRIL_ORGANIZATION

GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre <- meta$GOBP_COLLAGEN_FIBRIL_ORGANIZATION
GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post <- meta_post$GOBP_COLLAGEN_FIBRIL_ORGANIZATION

df4 <- data.frame(GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre,GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post)
df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre <- as.numeric(df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre)
df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post <- as.numeric(df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post)

pathway4 <- ggscatter(df4, x = "GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre", y = "GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post",
                      cor.method = "spearman",
                      add = "reg.line", add.params = list(color = "red")) +
  stat_cor(method = "spearman", label.x = min(df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_pre),label.y = max(df4$GOBP_COLLAGEN_FIBRIL_ORGANIZATION_post)) +
  geom_smooth(method = "lm", color = "red", fill = "orange", alpha = 0.8)

#### selected pathway1,pathway3 and pathway4

SpatialFeaturePlot(obj_sub,features = c("GOBP_SYNAPTIC_SIGNALING","GOBP_GENERATION_OF_NEURONS","GOBP_COLLAGEN_FIBRIL_ORGANIZATION"))

p1 <- SpatialFeaturePlot(obj_sub,features = c("GOBP_SYNAPTIC_SIGNALING"))

image_names <- names(obj_sub@images)
p_list_pre <- list()
p_list_post <- list()
# Loop through each image and generate a plot
for (image_name in image_names) {
  p1 <- SpatialFeaturePlot(not_scaled_tmp_sub,feature = c("GOBP_COLLAGEN_FIBRIL_ORGANIZATION"),images = image_name,alpha =c(0.8,1))
  p_list_post[[image_name]] <- p1
}

pdf("GOBP_COLLAGEN_FIBRIL_ORGANIZATION_spatial_plot_post_data.pdf",width =25,height =30)
cowplot::plot_grid(plotlist = p_list_post)
dev.off()
####

meta_common <- meta[, common_pathways, drop = FALSE]
meta_post_common <- meta_post[, common_pathways, drop = FALSE]

# Initialize dataframe to store correlations
cor_results <- data.frame(Pathway = common_pathways, cor_results = NA)

# Calculate Spearman correlation for each pathway
for (i in seq_along(common_pathways)) {
  path <- common_pathways[i]
  cor_results$cor_results[i] <- cor(meta_common[[path]], meta_post_common[[path]], method = "spearman", use = "pairwise.complete.obs")
}

write.csv(cor_results, "common_pathways_correlation_before_and_after_OrthologAL_app.csv", row.names = T)



# length(pathways$GOBP_SYNAPTIC_SIGNALING)
# [1] 796
# > length(pathways_mouse$GOBP_SYNAPTIC_SIGNALING)
# [1] 770
# length(pathways_mouse$GOBP_GENERATION_OF_NEURONS)
# [1] 1486
# > length(pathways$GOBP_GENERATION_OF_NEURONS)
# [1] 1510
# length(pathways_mouse$GOBP_COLLAGEN_FIBRIL_ORGANIZATION)
# [1] 66
# > length(pathways$GOBP_COLLAGEN_FIBRIL_ORGANIZATION)
# [1] 67
#
