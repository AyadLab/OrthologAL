#to process data fromn GEO Database,extract thebfiles and create a seurat object

library(GEOquery)
library(limma)
library(R.utils)
library(dittoSeq)
library(MAST)
library(GEOquery)
# Download the supplementary files for GSE129730
getGEOSuppFiles("GSE213240_RAW")
gunzip("GSE213240/GSE129730_series_matrix.txt.gz", remove = FALSE)
untar("GSE213240_RAW.tar")

# Load the matrix into R
matrix_data <- read.delim("GSE129730/GSE129730_series_matrix.txt", header = TRUE, row.names = 1)
#track all .gz files in the directory
gz_files <- list.files("/data/ayad_lab/app_analysis/GSE129730_RAW", full.names = TRUE)
lapply(gz_files, gunzip, remove = FALSE) 

gz_files <- list.files("/data/ayad_lab/app_analysis/GSE129730_RAW", full.names = TRUE)
lapply(gz_files, gunzip, remove = FALSE) 


dge_files <- list.files("/data/ayad_lab/Rat_SCI_data/GSE213240", full.names = TRUE)

# Create a list to store each Seurat object
seurat_list <- list()

# Loop over each file and create a Seurat object
for (i in seq_along(dge_files)) {
  dge_data <- read.delim(dge_files[i], header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  dge_matrix <- as.matrix(dge_data)
  
  # Create Seurat object
  medullo <- CreateSeuratObject(counts = dge_matrix, project = "medullo_scRNA")
  
  # Store the Seurat object in the list
  seurat_list[[i]] <- medullo
}

medullo <- NormalizeData(medullo)
medullo <- ScaleData(medullo)
medullo <- FindVariableFeatures(medullo)
medullo <- RunPCA(medullo,dims =1:30)
medullo <- RunUMAP(medullo,dims =1:30)
medullo <- FindNeighbors(medullo, dims = 1:30,features = rownames(medullo))
medullo <- FindClusters(medullo,algorithm = 1)


med_cluster <-as.data.frame(Idents(medullo))
omed_cluster <- as.data.frame(Idents(o_med))
o_med$med_cluster <- med_cluster
o_med$omed_cluster <- omed_cluster


m1 <- DimPlot(med) + theme_void()
m2 <- DimPlot(o_med) + theme_void()
pdf("ortho_medulllosc_UMAP.pdf")
cowplot::plot_grid(m1, m2, labels = c("before", "after"))
dev.off()

Idents(medullo) <- "seurat_clusters"
Idents(o_med) <- "seurat_clusters"
m_list <- o_med@meta.data %>%
  group_by(med_cluster,omed_cluster) %>%
  tally() %>%
  ungroup() %>%
  gather_set_data(1:2) %>%
  dplyr::mutate(
    x = factor(med_cluster, levels = unique(med_cluster)),
    y = factor(omed_cluster, levels = unique(omed_cluster))
  ) %>% as.data.frame()


ps1 <- ggplot(m_list, 
              aes(y = n, axis1 = med_cluster, axis2 = omed_cluster)) + 
  geom_alluvium(aes(fill = med_cluster), width = 1/10) + 
  geom_stratum(width = 1/10, fill = "black", color = "grey") + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) + 
  scale_x_discrete(limits = c( "before", "after"), expand = c(0.05, 0.05)) + 
  scale_fill_manual(values = colors) +
  
  #scale_fill_brewer(type = "qual", palette = colors) + 
  theme_bw()
pdf("cluster_mapping_mousemedullo_afterconversion.pdf",height =35,width =15)
ps1
dev.off()



o_med <- FindVariableFeatures(o_med)
o_med <- RunPCA(o_med,dims =1:30)
o_med <- RunUMAP(o_med,dims =1:30)
o_med <- FindNeighbors(o_med, dims = 1:30,features = rownames(medullo))
o_med <- FindClusters(o_med,algorithm = 1)



cluster_map <- c("0" = "A",
                "1" = "C",
                "2" = "B" ,
                "3" = "D",
                "4" = "E",
                "5" = "F",
                "6" = "G",
                "7" = "H",
                "8" = "I",
                "9" = "J"
                    
)

omed <- RenameIdents(omed, cluster_map)
clusters_1 <- as.numeric(Idents(medullo))
clusters_2 <- as.numeric(Idents(o_med))

rand_index_mouse <- rand.index(clusters_1,clusters_2)
# rand_index_mouse
# [1] 0.9306893
c <- med@assays$RNA$counts
gpc <- Matrix::colSums(c > 0)
gpc <- as.data.frame(gpc)

co <- omed@assays$RNA_ortho$counts
gpc_o <- Matrix::colSums(co > 0)
gpc_o <- as.data.frame(gpc_o)

df_retained <- cbind(gpc,gpc_o)
df_retained <- as.data.frame(df_retained)
df_retained$percent_retained <- (df_retained$gpc_o / df_retained$gpc) * 100

omed$percent_retained <- df_retained$percent_retained


vlnplot_percentretained <- VlnPlot(omed,features = "percent_retained",pt.size = 0)

pdf("vln_ratsci_percent_retained.pdf")
vlnplot_percentretained
dev.off()



