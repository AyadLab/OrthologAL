## Saving the Ortho table of genes of Rat,Mouse and Zebrafish to avoid using BioMart server
### as this server can be not a reliable option as we want to make the app locally available and independant of BioMart
library(biomaRt)
library(dplyr)

mart_species_mouse <- useEnsembl(biomart = "genes",
                                 dataset = "mmusculus_gene_ensembl")
mart_human <- useEnsembl(biomart = "genes",
                         dataset = "hsapiens_gene_ensembl")
df_mouse <- getBM(attributes = c("mgi_symbol","gene_biotype"),
                  mart = mart_species_mouse)
df_human <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                  mart = mart_human)

#  Prepare Mouse data: Create a "JOIN_KEY" in ALL CAPS
df_mouse_clean <- df_mouse %>%
  mutate(JOIN_KEY = toupper(mgi_symbol))

#  Prepare Human data: Create a "JOIN_KEY" in ALL CAPS
df_human_clean <- df_human %>%
  mutate(JOIN_KEY = toupper(hgnc_symbol))

#  Merge based on the name (JOIN_KEY) instead of the IDs
ortho_final_df <- inner_join(df_mouse_clean, df_human_clean, by = "JOIN_KEY") %>%
  dplyr::select(
    mgi_symbol,
    hgnc_symbol,
    gene_biotype,
    ensembl_gene_id = ensembl_gene_id.x # Keeping the Mouse ID for Seurat
  ) %>%
  distinct()

# Check the results
head(ortho_final_df)

ortho_final_df$ensembl_gene_id <- NULL

colnames(ortho_final_df) <- c("MGI.symbol", "HGNC.symbol", "Gene.type")
# Save this for the OrthologAL app

saveRDS(ortho_final_df, "data/ortho_df_Mouse_Human.rds")

# 1. Load/Connect (Note: Use archive if main is down, or use the DFs you already have)
# Already have the df_human from your previous step, skip re-fetching it
mart_rat <- useEnsembl(biomart = "genes", dataset = "rnorvegicus_gene_ensembl")

df_rat <- getBM(attributes = c("rgd_symbol", "gene_biotype"), mart = mart_rat)

#Clean and Match
df_rat_clean <- df_rat %>%
  mutate(JOIN_KEY = toupper(rgd_symbol))

df_human_clean <- df_human %>%
  mutate(JOIN_KEY = toupper(hgnc_symbol))

# Join by Name
ortho_rat_human <- inner_join(df_rat_clean, df_human_clean, by = "JOIN_KEY") %>%
  dplyr::select(
    RGD.symbol = rgd_symbol,
    HGNC.symbol = hgnc_symbol,
    Gene.type = gene_biotype
  ) %>%
  distinct()

ortho_rat_human$ensembl_gene_id <- NULL

saveRDS(ortho_rat_human, "data/ortho_df_Rat_Human.rds")
### Zebrafish
#  Fetch Zebrafish data
mart_zebrafish <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")

df_zebrafish <- getBM(attributes = c("zfin_id_symbol", "gene_biotype"), mart = mart_zebrafish)

#  Clean and Match
df_zebrafish_clean <- df_zebrafish %>%
  mutate(JOIN_KEY = toupper(zfin_id_symbol))

# Join by Name
ortho_zebrafish_human <- inner_join(df_zebrafish_clean, df_human_clean, by = "JOIN_KEY") %>%
  dplyr::select(
    ZFIN.symbol = zfin_id_symbol,
    HGNC.symbol = hgnc_symbol,
    Gene.type = gene_biotype
  ) %>%
  distinct()
ortho_zebrafish_human$ensembl_gene_id <- NULL

#  Save the rds file
saveRDS(ortho_zebrafish_human, "data/ortho_df_Zebrafish_Human.rds")
