### Import required libraries
library(tidyverse)
library(DreamAI)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(MSnSet.utils)
library(PCAtools)
library(limma)
library(viridis)
library(gprofiler2)
library(rstudioapi)
#######################


########## Required functions
run_limma_from_seurat <- function(
    seurat_obj,
    assay        = "Prot",
    slot         = "data",
    group_col    = "celltype",
    group1,
    group2
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(limma)
    library(Matrix)
  })
  
  # --- checks ---
  stopifnot(group_col %in% colnames(seurat_obj@meta.data))
  stopifnot(all(c(group1, group2) %in% seurat_obj@meta.data[[group_col]]))
  stopifnot(assay %in% Assays(seurat_obj))
  
  # Drop unused levels
  seurat_obj[[group_col]] <- droplevels(factor(seurat_obj[[group_col]][,1]))
  
  # Subset by metadata safely
  md <- seurat_obj@meta.data
  cells_use <- rownames(md)[md[[group_col]] %in% c(group1, group2)]
  
  seurat_sub <- subset(seurat_obj, cells = cells_use)
  
  # --- expression matrix (proteins × cells) ---
  expr <- GetAssayData(
    seurat_sub,
    assay = assay,
    slot  = slot
  )
  
  expr <- as.matrix(expr)  # limma needs dense
  
  # --- design matrix ---
  group <- factor(seurat_sub[[group_col]][,1], levels = c(group1, group2))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c(group1, group2)
  
  # --- contrast ---
  contrast <- makeContrasts(
    contrasts = paste0(group2, " - ", group1),
    levels = design
  )
  
  # --- limma ---
  fit <- lmFit(expr, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)
  
  res <-  topTable(
    fit,
    coef = paste0(group2, " - ", group1),
    number = Inf
  )
  
  return(res)
}

#####################
compute_hypergeo_pvals_at_threshold <- function(df, threshold) {
  # Filter by prediction score
  df_filt <- df %>% 
    dplyr::filter(prediction.score.max >= threshold)
  
  # Total number of high-confidence cells (N)
  N <- nrow(df_filt)
  if (N == 0) {
    return(
      tibble::tibble(
        cluster = unique(df$cluster),
        predicted.id = NA_character_,
        k = NA_integer_,
        K = NA_integer_,
        n = NA_integer_,
        N = 0,
        threshold = threshold,
        p_value = NA_real_,
        p_adj_BH = NA_real_
      )
    )
  }
  
  # Counts from high-confidence cells
  cluster_counts <- df_filt %>% 
    dplyr::count(cluster, name = "n")
  
  celltype_counts <- df_filt %>% 
    dplyr::count(predicted.id, name = "K")
  
  # Hypergeometric test
  results <- df_filt %>%
    dplyr::count(cluster, predicted.id, name = "k") %>%
    dplyr::left_join(celltype_counts, by = "predicted.id") %>%
    dplyr::left_join(cluster_counts, by = "cluster") %>%
    dplyr::mutate(
      N = N,
      threshold = threshold,
      p_value = phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    )
  
  # Ensure all clusters appear
  all_clusters <- tibble::tibble(cluster = unique(df$cluster))
  results <- all_clusters %>%
    dplyr::left_join(results, by = "cluster")
  
  # BH correction (across all tests at this threshold)
  results <- results %>%
    dplyr::mutate(
      p_adj_BH = p.adjust(p_value, method = "bonferroni")
    )
  
  return(results)
}
################################ Set wd and Import Data

# Get the directory of the current source file
script_dir <- dirname(getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)


########### Load subset of 20,000 PBMCs
integrated_data_new.1 <- readRDS(file = "./Data/scRNAseq_RefData_20k_PBMC.rds")


######## FragPipe output for all TMT channels

df <- read_tsv("./Data/abundance_protein_MD.tsv")

#### Dual column system metadata

LC_Sample_Meta <- read.csv("./Metadata/LC_Column_Meta.csv") %>%
  mutate(File_Name = gsub("_rerun", "",File_Name))


########### Estimated pg data, for later use
pg_ref <- readRDS("./Data/PBMC_pg_peptide_estimate.rds")

########################
prot_meta <- df %>%
  distinct(Index, Gene) %>%
  arrange(Index)

######### Pivot longer, Create Channel Types, filter for single-cell samples, create Batches
######### Create target(CHIP) groups
######### Join with meta data

df_long <- df %>%
  pivot_longer(cols = -1:-11, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(SampleID = gsub("_rerun", "", SampleID)) %>%
  mutate(Intensity = as.numeric(log2(Intensity))) %>%
  mutate(
    Channel_Type = case_when(
      str_detect(SampleID, "(^|[^A-Z0-9])127C([^A-Z0-9]|$)") | str_detect(SampleID, "(^|[^A-Z0-9])134ND([^A-Z0-9]|$)") ~ "Isotope",
      str_detect(SampleID, "(^|[^A-Z0-9])135CD([^A-Z0-9]|$)") |
        str_detect(SampleID, "(^|[^A-Z0-9])134C([^A-Z0-9]|$)") |
        str_detect(SampleID, "(^|[^A-Z0-9])135N([^A-Z0-9]|$)") ~ "Blank",
      str_detect(SampleID, "RefInt|RefDInt|ReferenceIntensity") ~ "Bridge",
      TRUE ~ "Single-Cell"
    )) %>%
  filter(Channel_Type %in% c("Single-Cell", "Isotope")) %>%
  mutate(File_Name = sub("_[^_]+$", "", SampleID),) %>%
  mutate(Channel = str_extract(SampleID, "[^_]+$")) %>%
  mutate(Chip = case_when(grepl("Target3", SampleID) ~ "Target3",
                            grepl("Target4_A1_", SampleID) ~ "Target3",
                            grepl("Target4_A2_", SampleID) ~ "Target3",
                            TRUE ~ "Target4")) %>%
  mutate(Well = sub(".*_(.*)$", "\\1", File_Name)) %>%
  left_join(.,LC_Sample_Meta) %>%
  filter(!is.na(Intensity)) %>%
  filter(!grepl("KRT", Gene)) %>%
  group_by(SampleID,Channel_Type) %>%
  add_count(name = "Protein_Obs") %>%
  mutate(LC_Column = paste("LC", LC,sep = "")) %>%
  ungroup() %>%
  select(-ReferenceIntensity)


###############

protein_wide <- df_long %>%
  select(Gene, SampleID, Intensity) %>%
  ungroup() %>%
  spread(SampleID, Intensity) %>% 
  column_to_rownames(var= "Gene") 

############


###### Full meta data for all TMT channels
meta_all <- df_long %>%
  distinct(File_Name, SampleID, Channel_Type, LC_Column, Channel, Well, Chip) %>%
  arrange(match(SampleID, colnames(protein_wide))) %>%
  mutate(File_Name = paste(File_Name, ".raw", sep = "" )) %>%
  group_by(File_Name) %>%
  mutate(TMT_Plex_ID = cur_group_id()) %>%
  column_to_rownames(var = "SampleID") %>%
  mutate(across(where(is.character), as.factor))  




#### Correct batches in order of largest to smallest impact
protein_BC <- protein_wide %>%
  MSnSet.utils::ComBat.NA(., meta_all$LC_Column) 

protein_BC <- protein_BC$`corrected data` %>%
  MSnSet.utils::ComBat.NA(., meta_all$Chip) 

protein_BC <- protein_BC$`corrected data` %>%
  MSnSet.utils::ComBat.NA(., meta_all$Channel_Type) 

protein_BC <- protein_BC$`corrected data` %>%
  MSnSet.utils::ComBat.NA(., meta_all$Channel) 

protein_BC <- protein_BC$`corrected data` %>%
  MSnSet.utils::ComBat.NA(., meta_all$Well) 

protein_BC <- protein_BC$`corrected data` 

protoP_impute <-  DreamAI(protein_BC,
                          k = 10, method = "KNN", 
                          out = c("KNN")) 

protoP_impute <- as.data.frame(protoP_impute$Ensemble)

######

###########################
#Transferring cell identities to proteomics dataset

seurat_protein <- CreateSeuratObject(counts = as.sparse(protoP_impute),
                                     min.features = 0)


seurat_protein <-  SetAssayData(seurat_protein,
                           layer = "data",
                           seurat_protein@assays$RNA$counts)

protoP__scale <- protoP_impute %>%
  t() %>%
  scale() %>%
  t() 

seurat_protein <-  SetAssayData(seurat_protein,
                                layer = "scale.data",
                                protoP__scale)


seurat_protein <- FindVariableFeatures(seurat_protein, 
                                       slot = "scale.data",
                                       selection.method = "vst",
                                       mean.cutoff = c(0,30),
                                       nfeatures = 600)

overlapping <- intersect(VariableFeatures(seurat_protein),
                         integrated_data_new.1@assays[["integrated"]]@var.features)

seurat_protein <- RunPCA(seurat_protein)

############################
pca_j<- PCAtools::pca(protoP_impute[VariableFeatures(seurat_protein), ],
                      metadata = meta_all,scale = TRUE, center = T) 


# ### 
#  sup1a <- eigencorplot(pca_j,
#              components = getComponents(pca_j, 1:12), 
#              metavars = colnames(meta_all),
#              scale = F)  
# 
# png("Supplementary_Fig1a.png", width = 750, height = 400)
# sup1a
# dev.off()

############################
###  18 PCs are selected from ref dataset based on elbowplot


anchors_transfer <- FindTransferAnchors(
  reference = integrated_data_new.1, 
  features = overlapping,
  query = seurat_protein, 
  max.features = 500,
  dims = 1:18,
  npcs = 18,
  reduction = "cca",
)

dim1 <- c(1:14)

predictions <- TransferData(
  anchorset = anchors_transfer,
  refdata = integrated_data_new.1$cell_type,
  dims = dim1, 
  weight.reduction = "cca"
)


seurat_protein <- AddMetaData(seurat_protein, metadata = predictions)
seurat_protein$cell_type <- seurat_protein$predicted.id


################# Mapping IDs back to proteomics UMAP


neighbors_k <- 20
seurat_protein <- FindNeighbors(seurat_protein, dims = dim1,
                                k.param = neighbors_k)
seurat_protein <- FindClusters(seurat_protein, 
                               resolution = 0.5, 
                               algorithm =1)

seurat_protein <- RunUMAP(seurat_protein, dims = dim1,
                          n.neighbors = neighbors_k,
                          metric = "euclidean")

DimPlot(seurat_protein, reduction = "umap", 
        pt.size = 1,
       group.by = "seurat_clusters"
)

DimPlot(seurat_protein, reduction = "umap", 
        pt.size = 1,
        group.by = "cell_type"
)


 predictions$cluster <- seurat_protein@meta.data$seurat_clusters

#######################


celltypes <- predictions %>%
  mutate(SampleID = rownames(.)) %>%
  compute_hypergeo_pvals_at_threshold(., threshold = 
                                        mean(predictions$prediction.score.max)-
                                        sd(predictions$prediction.score.max)) %>%  ### Main parameter to set
  group_by(cluster) %>%
  slice_min(order_by = p_value, n =1, with_ties = F) %>%
  mutate(predicted.id = case_when(p_value > 0.01 ~ "Unknown",
                                  is.na(p_value) ~ "Unknown",
                                  TRUE ~ predicted.id)) %>%
  ungroup() %>%
  mutate(Cell_Type_Final = factor(paste(predicted.id, cluster, sep = "_"))) %>%
  mutate(Cell_Type_Final= fct_reorder(Cell_Type_Final, as.integer(cluster))) %>%
  select(cluster, Cell_Type_Final)

predictions2 <- full_join(predictions, celltypes) %>%
  mutate(SampleID = rownames(predictions)) %>%
  mutate(Cell_Type_Final = case_when(grepl("CD4T_2", Cell_Type_Final) ~ "CD4T_1",
                                     grepl("CD4T_0", Cell_Type_Final) ~ "CD4T_1_2",
                                     grepl("monocyte_3", Cell_Type_Final) ~ "monocyte_0_2",
                                     grepl("monocyte_1", Cell_Type_Final) ~ "monocyte_0",
                                    # TRUE ~ Cell_Type_Final
                                     TRUE ~ NA
                                     ))


seurat_protein  <- AddMetaData(seurat_protein, predictions2$Cell_Type_Final, col.name = "Cell_Type_Final")

 #### UMAP Plots 

# Step 4: Visualize clusters
p1 <- DimPlot(seurat_protein, reduction = "umap", 
       # group.by = "seurat_clusters",
        group.by = "Cell_Type_Final",
       pt.size = 1,
       cols = c("#B79F00","#FFB000","#F8766D","#DC267F")
        #group.by = "cell_type"
)

#pdf(file = "Figure_S3.pdf", height =7.5, width = 9)
p1
#dev.off()

#######################################

protein_compare <- df_long %>%
  distinct(SampleID, Channel_Type, LC_Column, Channel, Well, Chip, Protein_Obs) %>%
  arrange(match(SampleID, colnames(protein_wide)))


pg_ref <- pg_ref %>%
  mutate(SampleID = gsub("_rerun","", SampleID))

protein_compare <- full_join(protein_compare, pg_ref)
protein_compare$Cell_Type <- predictions2$Cell_Type_Final
protein_compare$prediction_score <- predictions2$prediction.score.max


protein_compare %>%
  filter(grepl("T", Cell_Type)) %>%
  ggplot()+
  aes(x = Estimated_pgcell,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 1)+
  scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Estimated peptide (pg)",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#B79F00","#FFB000"))+

protein_compare %>%
  filter(grepl("mono", Cell_Type)) %>%
  ggplot()+
  aes(x = Estimated_pgcell,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 1)+
  scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Estimated peptide (pg)",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#F8766D","#DC267F")) #+
 # ggsave("Figure_S4a.png",
   #      width = 8, height = 3)




protein_compare %>%
  filter(grepl("T", Cell_Type)) %>%
  ggplot()+
  aes(x = Protein_Obs,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 10)+
 # scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Protein Identifications (n)",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#B79F00","#FFB000"))+

  protein_compare %>%
  filter(grepl("mono", Cell_Type)) %>%
  ggplot()+
  aes(x = Protein_Obs,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 10)+
 # scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Protein Identifications (n)",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#F8766D","#DC267F"))# +
 # ggsave("Figure_S4b.png",
  #       width = 8, height = 3)



protein_compare %>%
  filter(grepl("T", Cell_Type)) %>%
  ggplot()+
  aes(x = prediction_score,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 0.02)+
  # scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Seurat Prediction Score",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#B79F00","#FFB000"))+

  protein_compare %>%
  filter(grepl("mono", Cell_Type)) %>%
  ggplot()+
  aes(x = prediction_score,
      fill = Cell_Type)+
  geom_histogram( alpha = 0.8, color  = "black", position = "identity",
                  binwidth = 0.02)+
  # scale_x_continuous(limits = c(0, 45), oob = scales::squish)+
  theme_minimal()+
  labs(
    x = "Seurat Prediction Scores",
    y = "Cells (n)"
  )+
  scale_fill_manual(values = c("#F8766D","#DC267F")) #+
  #ggsave("Figure_S4c.png",
   #      width = 8, height = 3)

########################

predictions3 <- predictions2 %>%
  filter(!is.na(predictions2$Cell_Type_Final)) %>%
  select(SampleID, Cell_Type_Final) %>%
  mutate(Cell_Type_Final = as.factor(Cell_Type_Final))


protein_BC2 <- protein_BC[,colnames(protein_BC) %in% predictions3$SampleID]


seurat_protein2 <- protein_BC2 %>%
  as.sparse() %>%
  CreateSeuratObject() %>%
  RenameAssays(.,
               assay.name = "RNA",
               new.assay.name = "Prot",
               verbose = TRUE) %>%
  AddMetaData(
    object = .,
    metadata = predictions3$Cell_Type_Final,
    col.name = 'celltype') %>%
  SetAssayData(., layer = "data",.@assays$Prot$counts)

TCELL <- run_limma_from_seurat(seurat_protein2, group1 = "CD4T_1", group2 = "CD4T_1_2")
TCELL_Down <- TCELL %>%
  filter(adj.P.Val <= 0.01 & logFC < 0
         )%>%
  mutate(Gene = rownames(.)) %>%
  distinct(Gene) %>%
  as.list()
monocyte <- run_limma_from_seurat(seurat_protein2, group1 = "monocyte_0", group2 = "monocyte_0_2")

monocyte_Down <- monocyte %>%
  filter(adj.P.Val <= 0.01 & logFC < 0
         )%>%
  mutate(Gene = rownames(.)) %>%
  distinct(Gene) %>%
  as.list()

bg <- rownames(protein_BC2)


TCELL_DOWN_GO<- gost(query = TCELL_Down,
                   organism = "hsapiens", ordered_query = F, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                   measure_underrepresentation = FALSE, evcodes = TRUE, 
                   user_threshold = 0.05, correction_method = "g_SCS", 
                   domain_scope = "annotated", custom_bg = bg, 
                   sources = "GO:CC",
                   numeric_ns = "", as_short_link = FALSE)


TCELL_1 <- TCELL_DOWN_GO$result %>%
  select(term_id, p_value) %>%
  distinct(p_value,.keep_all = T) %>%
  dplyr::slice_min(order_by = p_value, n = 3, with_ties = FALSE)  %>%
  pull(term_id)

monocyte_DOWN_GO<- gost(query = monocyte_Down,
                        organism = "hsapiens", ordered_query = F, 
                        multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.05, correction_method = "g_SCS", 
                        domain_scope = "annotated", custom_bg = bg, 
                        sources = "GO:CC",
                        numeric_ns = "", as_short_link = FALSE) 

monocyte_1 <- monocyte_DOWN_GO$result %>%
  select(term_id, p_value) %>%
  distinct(p_value,.keep_all = T) %>%
  dplyr::slice_min(order_by = p_value, n = 3, with_ties = FALSE)  %>%
  pull(term_id)

highlight_terms <- c(TCELL_1, monocyte_1)

p1 <- gostplot(TCELL_DOWN_GO, capped = F, interactive = F)

######## Figure S5a
#publish_gostplot(p1,highlight_terms = highlight_terms,
#                  height = 6, width = 7, filename = "Figure_S5a.png" )


p2 <- gostplot(monocyte_DOWN_GO, capped = F, interactive = F)

######## Figure S5b
#publish_gostplot(p2,highlight_terms = highlight_terms,
#                 height = 6, width = 7, filename = "Figure_S5b.png" )




################# Identify single cells to retain for further analysis

keep_datasets <- predictions2 %>%
  mutate(Cell_Type_Final = case_when(is.na(Cell_Type_Final) ~ "KEEP",
                                     TRUE ~ Cell_Type_Final)) %>%
  filter(Cell_Type_Final != "CD4T_1_2") %>%
  filter(Cell_Type_Final != "monocyte_0_2")  %>%
  distinct(SampleID)

 # write.csv(keep_datasets, file= "First_Pass_PBMCs_KEEP.csv",
  #          row.names = F)
  
############# Supplementary Table S5
  
supp_table <- meta_all %>%
  rownames_to_column(var = "SampleID")

supp_table <- predictions2 %>%
  select(SampleID, Cell_Type_Final) %>%
  full_join(., supp_table) %>%
  full_join(., pg_ref) %>%
  ungroup() %>%
  mutate(Estimated_pgcell = round(Estimated_pgcell, digits = 1) )

#write.csv(supp_table, file = "Supplementary_Table_S5.csv",
 #         row.names = F)  
  