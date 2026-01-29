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
library(gt)
library(scales)
library(rstudioapi)
#########################


########## Required functions

MyFindAllMarkers <- function(
    object,
    assay = "Prot",
    slot = "data",
    logfc.threshold = 0,
    verbose = TRUE,
    min_obs = 0,
    ...
) {
  # Ensure the assay is available in the Seurat object
  DefaultAssay(object) <- assay
  
  # Extract data from the specified slot
  data <- GetAssayData(object, slot = slot)
  check1 <- rowSums(is.na(data))
  check1 <- check1[check1 < ncol(data)*min_obs]
  data <- data[names(check1),]
  # Cluster information
  clusters <- Idents(object)
  
  # Initialize a data frame to store marker information
  all_markers <- data.frame()
  
  # Loop through each cluster and find markers
  for (cluster in levels(clusters)) {
    if (verbose) {
      message("Processing cluster: ", cluster)
    }
    
    # Cells in the current cluster
    cluster_cells <- WhichCells(object, idents = cluster)
    
    # Cells in other clusters
    other_cells <- setdiff(Cells(object), cluster_cells)
    
    # Create design matrix for limma
    design <- model.matrix(~ factor(c(rep(1, length(cluster_cells)), rep(0, length(other_cells)))))
    colnames(design) <- c("Other", "Cluster")
    
    # Combine expressions into a matrix
    expr_data <- cbind(data[, cluster_cells], data[, other_cells])
    
    cntrsts <- "Cluster-Other"
    
    # Fit linear model using limma
    fit <- lmFit(expr_data, design)
    made_contrasts <- makeContrasts(contrasts = cntrsts, levels = design)
    contrast_fit <- contrasts.fit(fit, made_contrasts)
    fit <- eBayes(fit)
    
    # Extract results
    results <- topTable(fit, coef = "Cluster", number = Inf, adjust.method = "bonferroni")
    results <- results[!is.na(results$t), ]
    # Calculate log2 fold-change
    log2FC <- results$logFC
    
    # Create a data frame of results
    markers_df <- data.frame(
      gene = rownames(results),
      cluster = cluster,
      log2FC_cluster = log2FC,
      p_value = results$P.Value,
      p_value_adj_cluster = results$adj.P.Val
    )
    
    # Filter based on logfc.threshold and min.pct
    markers_df <- markers_df %>%
      filter(log2FC > logfc.threshold)
    
    # Append to the all_markers data frame
    all_markers <- rbind(all_markers, markers_df)
  }
  
  return(all_markers)
}
##############


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
      p_adj = p.adjust(p_value, method = "bonferroni")
    )
  
  return(results)
}

################################ Set wd and Import Data

# Get the directory of the current source file
script_dir <- dirname(getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)


########### Load subset of 20,000 PBMCs
load(file = "./Data/scRNAseq_RefData_20k_PBMC.RData")

############## Intact PBMC cell identifiers to keep
df_filter <- read.csv(file = "./Metadata/First_Pass_PBMCs_KEEP.csv")

######## FragPipe output for all TMT channels

df <- read_tsv("./Data/abundance_protein_MD.tsv") 

#### Dual column system metadata

LC_Sample_Meta <- read.csv("./Metadata/LC_Column_Meta.csv") %>%
  mutate(File_Name = gsub("_rerun", "",File_Name))


########################
prot_meta <- df %>%
  distinct(Index, Gene) %>%
  arrange(Index)

######### Pivot longer, Create Channel Types, filter for single-cell samples, create Batches
######### Create target(CHIP) groups
######### Join with meta data, filter for minimum 600 protein observations

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
  left_join(., LC_Sample_Meta) %>%
  filter(!is.na(Intensity)) %>%
  filter(!grepl("KRT", Gene)) %>%
  filter(SampleID %in% df_filter$SampleID) %>%
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
  distinct(SampleID, Channel_Type, LC_Column, Channel, Well, Chip) %>%
  arrange(match(SampleID, colnames(protein_wide))) %>%
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
var_plot <- VariableFeaturePlot(seurat_protein,
                                pt.size = 2)
sup7 <- LabelPoints(plot = var_plot, points = overlapping, repel = T) 
sup7[["data"]] <- sup7[["data"]] %>%
  mutate(colors = case_when(row.names(.) %in% overlapping ~ "yes",
                            TRUE ~ "No"))

#pdf(file = "Sup_FigS7.pdf", height =8, width = 9)
sup7
#dev.off()
############################
pca_j<- PCAtools::pca(protoP_impute[VariableFeatures(seurat_protein), ],
                      metadata = meta_all,scale = TRUE, center = T) 


############################

seurat_protein <- RunPCA(seurat_protein)

###18 PCs are selected from ref dataset based on elbowplot

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

neighbors_k <- 13

seurat_protein <- FindNeighbors(seurat_protein, dims = dim1,
                                k.param = neighbors_k)
seurat_protein <- FindClusters(seurat_protein,
                               resolution = 0.5,
                               algorithm =1,
                               group.singletons = F
                               )
seurat_protein <- RunUMAP(seurat_protein, dims = dim1,
                          n.neighbors = neighbors_k,
                          metric = "euclidean")

#pdf(file = "Figure_S8a.pdf", height = 5, width = 5)
#DimPlot(seurat_protein, reduction = "umap", 
#        pt.size = 1,
#        group.by = "seurat_clusters"
       ### group.by = "cell_type"
#)
#dev.off()

umap_seurat <- seurat_protein@reductions[["umap"]]@cell.embeddings %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.),
         predicted.id = ifelse(predictions$prediction.score.max <= mean(predictions$prediction.score.max)-sd(predictions$prediction.score.max),
                                       NA,predictions$predicted.id)) %>%
  arrange(!is.na(predicted.id), predicted.id)



colrs <- c("#00BFC4",
           "#B79F00",
           "#53B400",
           "#619CFF",
           "#F8766D",
           "#00C094"
           )
#pdf(file = "Figure_S8b.pdf", height = 4, width = 6)
# ggplot(umap_seurat, aes(x = umap_1, y = umap_2, color = predicted.id)) +
#   geom_point(size = 1.5 ) +
#  # theme_minimal(base_size = 10) +
#   theme(
#     axis.line = element_line(color='black'),
#     panel.border = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.background = element_rect(fill = 'white'),
#     plot.background = element_rect(fill = "white")
#                        )+
#   xlab("umap_1") +
#   ylab("umap_2") +
#   theme() +
#   labs(color = "Seurat Predicted \n Cell Type \n (prediction score > 0.56)")+
#   scale_color_manual(values = colrs)
# dev.off()

# seurat_protein@meta.data %>%
#   select(cell_type, prediction.score.max) %>%
#   arrange(cell_type) %>%
#   ggplot()+
#   aes(x = prediction.score.max, fill = cell_type)+
#   geom_histogram( alpha = 0.6, color  = "black",
#                   binwidth = 0.025)+
#   geom_vline(xintercept = (mean(predictions$prediction.score.max)-sd(predictions$prediction.score.max)), linetype  = "dashed",
#              size = 0.7)+
#   labs(
#     x = "Prediction Score",
#     y = "Cells (n)"
#   )+
#   theme_minimal()+
#   scale_fill_manual(values = colrs)+
#   ggsave(file = "Figure_S6.png", height = 5, width = 8)

predictions$cluster <- seurat_protein@meta.data$seurat_clusters

#######################


fig_s8c <- predictions %>%
  mutate(SampleID = rownames(.)) %>%
  compute_hypergeo_pvals_at_threshold(mean(predictions$prediction.score.max)-sd(predictions$prediction.score.max)) %>%
  group_by(cluster)  %>%
 slice_min(order_by = p_adj, n =1, with_ties = F) %>%
  mutate(predicted.id = case_when(p_adj > 0.01 ~ "Unknown",
                                  is.na(p_adj) ~ "Unknown",
                                  TRUE ~ predicted.id)) %>%
  mutate(Cell_Type_Final = factor(paste(predicted.id, cluster, sep = "_"))) %>%
  mutate(Cell_Type_Final= fct_reorder(Cell_Type_Final, as.integer(cluster))) %>%
  select(cluster, predicted.id, p_adj, Cell_Type_Final) %>%
  mutate(cluster = as.character(cluster)) %>%
  ungroup()

 
 #colrs <- hue_pal()(7)
 # fig_s8c %>%
 #   mutate(p_adj = round(-log10(p_adj),2))%>%
 #  dplyr::rename(`-log10(p_adj)` = p_adj) %>%
 #  gt() %>%
 #    data_color(columns =  `-log10(p_adj)`,
 #               palette = "viridis") %>%
 #    data_color(columns = `cluster`,
 #               palette = colrs) %>%
 #    tab_style(
 #      style = cell_borders(sides = "all",
 #                           color = "#000000",
 #                           style = "solid",
 #                           weight = px(1)),
 #      locations = cells_body()) %>%
 #    cols_align(
 #      align = "center",
 #      columns = everything()
 #    ) %>%
 #    tab_options(table.font.size = 16) %>%
 #    gtsave("Figure_S8c.png",
 #           expand = 10)

celltypes <- predictions %>%
  mutate(SampleID = rownames(.)) %>%
  compute_hypergeo_pvals_at_threshold(.,
                                      threshold = mean(predictions$prediction.score.max)-sd(predictions$prediction.score.max)
                                          ) %>%
  group_by(cluster) %>%
  slice_min(order_by = p_adj, n =1, with_ties = F) %>%
  mutate(predicted.id = case_when(p_adj >= 0.01 ~ "Unknown",
                                  is.na(p_adj) ~ "Unknown",
                                  TRUE ~ predicted.id)) %>%
  ungroup() %>%
  mutate(Cell_Type_Final = factor(paste(predicted.id, cluster, sep = "_"))) %>%
  mutate(Cell_Type_Final= fct_reorder(Cell_Type_Final, as.integer(cluster))) %>%
  select(cluster, Cell_Type_Final)

predictions2 <- full_join(predictions, celltypes)

seurat_protein  <- AddMetaData(seurat_protein, predictions2$Cell_Type_Final, col.name = "Cell_Type_Final")

#### UMAP Plots 

# Step 4: Visualize clusters
#pdf(file = "Figure_4a.pdf", height = 5, width = 6)
#DimPlot(seurat_protein, reduction = "umap",
   #    group.by = "Cell_Type_Final",
   #    pt.size = 1
#)
#dev.off()
#######################################

Seurat_scale <- protein_BC %>%
  t() %>%
  scale() %>%
  t() %>%
  as.sparse() 

seurat_protein2 <- protein_BC %>%
  as.sparse() %>%
  CreateSeuratObject() %>%
  RenameAssays(.,
               assay.name = "RNA",
               new.assay.name = "Prot",
               verbose = TRUE) %>%
  AddMetaData(
    object = .,
    metadata = predictions2$Cell_Type_Final,
    col.name = 'Cluster') %>%
  SetAssayData(., layer = "data",.@assays$Prot$counts) %>%
  SetAssayData(layer = "scale.data",Seurat_scale) 

seurat_protein[["Prot"]] <- seurat_protein2@assays

seurat_protein <- SetIdent(seurat_protein, value = "Cell_Type_Final")

Prot.markers <- MyFindAllMarkers(seurat_protein, assay = "Prot",
                                 slot = "data",
                                 min_obs = 0.9) %>%
  filter(p_value_adj_cluster <= 0.001) %>%
  filter(log2FC_cluster >= 0.5) %>%
  ungroup()

#write.csv(Prot.markers, file = "Supplementary_Table_S3.csv",
#          row.names = F)

################ Established/canonical immune markers
manual_markers <- c("CD3D", "CD14", "CD36",
                    "IFI30",
                    "RHOC",
                    "DGKA",
                    "CD74",
                    "HLA_DRB1",
                    "WARS1",
                    "FCER1G",
                    "LTA4H",
                    "ASAH1",
                    "GIMAP4",
                    "CAMP", "ELANE", "CTSG", "MMP9", "LTF"
                    )
  
manual_markers2 <- Prot.markers %>%
  filter(gene %in% manual_markers) 
top5 <- Prot.markers %>% 
  group_by(cluster) %>% 
  slice_min(., order_by = p_value_adj_cluster, n = 3, 
            with_ties = F ) %>%
  full_join(., manual_markers2) %>%
  group_by(gene) %>%
  slice_max(., order_by = log2FC_cluster, n = 1, 
            with_ties = F ) %>%
  ungroup() %>%
  mutate(cluster= factor(cluster,
                         levels = levels(seurat_protein@active.ident))) %>%
  arrange(cluster,
          -log2FC_cluster)


heat2<- DoHeatmap(seurat_protein, assay  ="Prot", features = top5$gene,
                  slot = "scale.data",
                  group.bar.height = 0.05,
                  draw.lines = F,
                  size = 4,
                  label = T
) + scale_fill_viridis( na.value = "black")

#pdf(file = "Figure_4b.pdf", width = 12, height  = 7)
heat2
#dev.off()


############ For supplementary table S5


#second_cell_types <- as.data.frame(seurat_protein@active.ident) %>%
#  rownames_to_column(var = "SampleID")%>%
 # dplyr::rename(Cell_Type2 = "seurat_protein@active.ident")
