#pre-processing figures
library(Seurat)
library(biomaRt)
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

species_lookup <- data.frame(
Mouse = c(ensembl_id = "mmusculus_gene_ensembl", attributes = 'mgi_symbol', filters = 'mgi_symbol'),
Human = c(ensembl_id = "hsapiens_gene_ensembl", attributes = 'ensembl_gene_id', filters = 'ensembl_gene_id'),
Zebrafish = c(ensembl_id = "drerio_gene_ensembl", attributes = 'zfin_id_symbol', filters = 'zfin_id_symbol'),
Rat = c(ensembl_id = "rnorvegicus_gene_ensembl", attributes = 'rgd_symbol', filters = 'rgd_symbol'),
stringsAsFactors = FALSE
)
species_info <- species_lookup[["Mouse"]]

mart.species <- useEnsembl("ensembl", species_info[[1]], mirror = 'useast', host = "https://dec2021.archive.ensembl.org")
mart.human <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = 'useast', host = "https://dec2021.archive.ensembl.org")
mart.rat <- useEnsembl("ensembl", rat_info[[1]], mirror = 'useast', host = "https://dec2021.archive.ensembl.org")

#species-genes gives us all the genes available in the BioMart DB
species_genes <- getBM(
  attributes = c("ensembl_gene_id","gene_biotype","mgi_symbol"),
  uniqueRows = TRUE,
  mart = mart.species
) #21884 pco transcripts

rat_genes <- getBM(
  attributes = c("ensembl_gene_id","rgd_symbol","gene_biotype"),
  uniqueRows = TRUE,
  mart = mart.rat
) 
#rat
species_unique_genes_MOUSE <- species_genes[!duplicated(species_genes$mgi_symbol), ] # 21789 genes
m_class <- as.data.frame(table(species_unique_genes_MOUSE$gene_biotype))
#mouse
species_unique_genes_RAT <- rat_genes[!duplicated(rat_genes$ensembl_gene_id), ] # 21789 genes

human_unique_genes <- human_genes[!duplicated(human_genes$ensembl_gene_id), ] 
##### mouse_converted__hg is mouse species orthologous to human 
mouse_converted__hg <- biomaRt::getLDS(
  attributes =  c(species_info[[2]],"gene_biotype","ensembl_gene_id"),
  filters = species_info[[3]],
  values = as.character(species_genes$mgi_symbol),
  mart = mart.species,
  attributesL = c('hgnc_symbol'),
  martL = mart.human,
  uniqueRows = T
)
##### rat_converted__hg is mouse species orthologous to human 

rat_converted__hg <- biomaRt::getLDS(
  attributes =  c(rat_info[[2]],"gene_biotype","ensembl_gene_id"),
  filters = rat_info[[3]],
  values = as.character(rat_genes$rgd_symbol),
  mart = mart.rat,
  attributesL = c('hgnc_symbol'),
  martL = mart.human,
  uniqueRows = T
)

converted_unique_m_h <- mouse_converted__hg[!duplicated(mouse_converted__hg$HGNC.symbol), ]
converted_unique_r_h <- rat_converted__hg[!duplicated(rat_converted__hg$HGNC.symbol), ]

m_h_class <- as.data.frame(table(converted_unique_m_h$Gene.type))
colnames(m_h_class) <- c("Gene_Type", "Freq")
graph_r_h <- ggplot(r_h_class, aes(x = "", y = Freq, fill = Gene_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

pdf("graph1_BIOMART_rat_human_ortholgs.pdf")
graph_r_h
dev.off()

r_h_class <- as.data.frame(table(converted_unique_r_h$Gene.type))
colnames(r_h_class) <- c("Gene_Type", "Freq")
species_unique_genes <- species_genes[!duplicated(species_genes$ensembl_gene_id), ] # 21789 genes
Mouse_Data <- as.data.frame(table(species_unique_genes_MOUSE$gene_biotype))
Rat_Data <- as.data.frame(table(species_unique_genes$gene_biotype))
colnames(Rat_Data) <- c("Gene_Type", "Freq")

Human_Data <- as.data.frame(table(human_unique_genes$gene_biotype))
colnames(Human_Data) <- c("Gene_Type", "Freq")
# Assuming Rat_Data is your dataframe
total_freq <- sum(Rat_Data$Freq)
# total_freq
# [1] 30560
colnames(Mouse_Data) <- c("Gene_Type", "Freq")


human_graph <- ggplot(Human_Data, aes(x = "", y = Freq, fill = Gene_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

unique_hg <- human_genes[!duplicated(human_genes$hgnc_symbol), ]

pdf("human_graph_biomart_genes.pdf")
human_graph
dev.off()

gen_hg <- getBM(
  attributes = c('gene_biotype', 'hgnc_symbol','ensembl_gene_id'),
  filters = 'hgnc_symbol',
  values = gen,
  mart = mart.human
)
#human genes detected after conversion and rat is the the dataset fro reference and orat is the seurat object after conversion using the app
c <- rat@assays$RNA$counts
genes_matrix <- as.data.frame(Matrix::rowSums(c > 0))
c_o <- orat@assays$RNA_ortho$counts
ortho_genes_matrix <- as.data.frame(Matrix::rowSums(c_o > 0))

c_o <- orat@assays$RNA_ortho$counts
ortho_genes_matrix <- as.data.frame(Matrix::rowSums(c_o > 0))

# 1 Human Detected   64.65737 MB scrna seq data
# 2   Not Detected   35.34263
human_detected <- nrow(ortho_genes_matrix) / nrow(genes_matrix) * 100
#human_detected <- sum(nrow(converted_unique)+ nrow(gen_hg_unique)) / sum(nrow(genobj_unique)+ nrow(gen_hg_unique)) * 100
not_detected <- 100 - human_detected
data_detected <- data.frame(
  category = c("Human Detected", "Not Detected"),
  percentage = c(human_detected, not_detected)
)

p1 <- ggplot(data_detected, aes(x = "", y = percentage, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Category") +
  theme_void() +
  theme(legend.title = element_blank())


# data_detected for MB scRNAseq
# category percentage
# 1 Human Detected   71.81613
# 2   Not Detected   28.18387

# data_detected PDOX data
# category percentage
# 1 Human Detected   76.63479
# 2   Not Detected   23.36521
pdf("humangenes_detect_from_med.pdf")
p1
dev.off()
# data_detected for tabula dataset
# category percentage
# 1 Human Detected     57.465
# 2   Not Detected     42.535

gene_classification <- as.data.frame(table(converted_unique$Gene.type))
colnames(gene_classification) <- c("Gene_Type", "Freq")
#17,485 pco genes were converted

graph2 <- ggplot(gene_classification, aes(x = "", y = Freq, fill = Gene_Type)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", margin = margin(b = 10)),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12))

pdf("graph2-BIOMART_rat_SCI_gene.pdf")
graph2
dev.off()

c <- tab@assays$RNA$counts
gpc <- Matrix::colSums(c > 0)
gpc <- as.data.frame(gpc) 
co <- orat@assays$RNA_ortho$counts
gpc_o <- Matrix::colSums(co > 0)
gpc_o <- as.data.frame(gpc_o)
df_retained <- cbind(gpc,gpc_o)
df_retained <- as.data.frame(df_retained)
df_retained$percent_retained <- (df_retained$gpc_o / df_retained$gpc) * 100
orat$percent_retained <- df_retained$percent_retained

med$percent_retained <- gpc

ps1 <- DimPlot(med, group.by = "percent_retained",label = FALSE) + 
  scale_color_viridis_d()


ps1 <- DimPlot(orat, group.by = "seurat_clusters",split.by = "percent_retained",raster = FALSE) + 
  scale_color_viridis_d() + 
  theme(legend.text = element_blank(), legend.title = element_blank())

pdf("figure1_data.pdf",width = 80)
ps1
dev.off()

# #gen extracts gene names from seurat object uploaded 
# gen <- rownames(tmp)
# gen <- sub("^hg38-|^mm10-", "", gen)
# 
# gen_obj <- getBM(
#   attributes = c('gene_biotype', 'mgi_symbol','ensembl_gene_id'),
#   filters = species_info[[2]],
#   values = gen,
#   mart = mart.species
# )
# gen_pco <- gen_obj[gen_obj$gene_biotype == "protein_coding",] 
# 
# genobj_unique_tab <- gen_obj[!duplicated(gen_obj$ensembl_gene_id), ]
# 
# gen_unique_hg <- gen_hg[!duplicated(gen_hg$ensembl_gene_id), ]
# unqiue_ensembl_human <-  as.data.frame(table(genobj_unique_hg$gene_biotype))  #protein_coding 21236
# gclass_hg_data <- as.data.frame(table(genobj_unique_tab$gene_biotype))
#protein_coding 18837
# length(gen_pco_hg$ensembl_gene_id)
# [1] 38487
# length(gen_pco$ensembl_gene_id) (PDOX mgi data)
# [1] 36905

# length(gen_pco$ensembl_gene_id)
# [1] 21446 (tab data)
# 
# length(gen_pco$ensembl_gene_id)
# [1] 19305 (rat SCI data)
# genobj_unique <- gen_obj[!duplicated(gen_obj$mgi_symbol), ] #29,751
# gen_hg_unique <- gen_hg[!duplicated(gen_hg$hgnc_symbol), ] #22,840
# 
# gene_classification_Data <- as.data.frame(table(genobj_unique$gene_biotype))
# colnames(gene_classification_Data) <- c("Gene_Type", "Freq")



converted_all <- biomaRt::getLDS(
  attributes =  c(species_info[[2]],"gene_biotype","ensembl_gene_id"),
  filters = species_info[[3]],
  values = as.character(species_genes$mgi_symbol),
  mart = mart.species,
  attributesL = c('hgnc_symbol'),
  martL = mart.human,
  uniqueRows = T
)
genobj_unique <- gen_obj[!duplicated(gen_obj$mgi_symbol), ]
uconverted_all <- converted_all[!duplicated(converted_all$HGNC.symbol), ]

genobj_unique <- gen_obj[!duplicated(gen_obj$mgi_symbol), ] #29,751
gen_hg_unique <- gen_hg[!duplicated(gen_hg$hgnc_symbol), ] #22,840

