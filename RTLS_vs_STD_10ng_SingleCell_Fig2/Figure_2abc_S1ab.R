library(tidyverse)
library(SingleCellExperiment)
library(scp)
library(rstudioapi)
#######################

# Get the directory of the current source file
script_dir <- dirname(getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)

### Import data/metadata

df <- read_tsv("./Data/abundance_peptide_None.tsv") %>%
  pivot_longer(-Index,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  filter(!grepl("_SC_", SampleID)) %>%
  filter(!grepl("pg_", SampleID)) %>%
  mutate(SampleID = gsub("RefDInt_", "",SampleID)) %>%
  mutate(SampleID = gsub("RefInt_", "",SampleID)) %>%
  mutate(SampleID = gsub("_135ND", "",SampleID)) %>%
  mutate(SampleID = gsub("_126", "",SampleID)) %>%
  filter(!is.na(Intensity)) %>%
  group_by(SampleID, Index) %>%
  summarize(Intensity = mean(Intensity)) %>%
  ungroup() %>%
  mutate(Method = case_when(grepl("RTLS", SampleID) ~ "RTLS",
                            grepl("STD_MS2", SampleID) ~ "STD_MS2"))

metadata <- df %>%
  distinct(SampleID, Method)

################ Generate Figures

df2 <- df %>%
  dplyr::select(SampleID, Index, Intensity) %>%
  pivot_wider(values_from = "Intensity",
              names_from = "SampleID") %>%
  column_to_rownames(var = "Index") %>%
  as.matrix()
  

sce <- SingleCellExperiment(assays = List(df2),)
sim <- QFeatures(experiments = List(Assay1 = sce))

metadata <- metadata %>%
  arrange(SampleID,sim@colData@rownames)

sim$SampleType <- metadata$Method


data_completedness <- reportMissingValues(sim, "Assay1", by = sim$SampleType) %>%
  mutate(Completeness = Completeness*100) %>%
  mutate(Type = rownames(.))

dt1 <- data_completedness %>%
  filter(Type == "RTLS")
dt2 <- data_completedness %>%
  filter(Type == "STD_MS2")

csc <- cumulativeSensitivityCurve(
  sim, "Assay1", by = sim$SampleType,
) %>%
  mutate(Type = by)
csc %>%
  ggplot()+
  aes(x = SampleSize, y = Sensitivity, colour = Type) +
  geom_point(size =0.5)+
  geom_hline(yintercept = dt1$TotalSensitivity,
             color = "#F8766D")+
  geom_hline(yintercept = dt2$TotalSensitivity,
             color = "#00BFC4")+
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom",
    panel.background = element_rect(fill= 'white'),
    axis.text.y=element_text(color = 'black'),
    axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line())+
  ylab("Sensitivity")+
  xlab("Sample Size")+
  scale_y_continuous(limits = c(0,9500))+
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10))+
  labs(colour = "Method")
 #+
 # ggsave("Figure_S1a.png", width = 3.3, height = 2.4)


data_completedness %>%
  ggplot()+
  aes(x = Completeness, y = LocalSensitivityMean, color = Type,
      group = Type)+
  geom_point(size =6,
             shape = 5)+
  theme_bw()+
  scale_x_continuous(limits = c(10,55))+
  scale_y_continuous(limits = c(2000,5500))+
  geom_pointrange(aes(ymin=LocalSensitivityMean-LocalSensitivitySd,
                      ymax=LocalSensitivityMean+LocalSensitivitySd),
                  position=position_dodge(0.05),
                  size = 0.1)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
    panel.background = element_rect(fill= 'white'),
    axis.text.y=element_text(color = 'black'),
    axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black'),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line())+
  ylab("Local Sensitivty")+
  xlab("Data completeness (%)") #+
 # ggsave("Figure_S1b.png",
    #     height = 2,
    #     width = 3.3)

###########################################


df2 <- read_tsv("./Data/abundance_peptide_None.tsv") %>%
  pivot_longer(-Index,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  filter(grepl("_SC_", SampleID)) %>%
  filter(!is.na(Intensity)) %>%
  mutate(Channel_Type = case_when(
    str_detect(SampleID, "(^|[^A-Z0-9])135CD([^A-Z0-9]|$)") |
      str_detect(SampleID, "(^|[^A-Z0-9])134C([^A-Z0-9]|$)") |
      str_detect(SampleID, "(^|[^A-Z0-9])135N([^A-Z0-9]|$)") ~ "Blank",
    str_detect(SampleID, "RefInt|RefDInt") ~ "Bridge",
    TRUE ~ "Single-Cell")) %>%
  filter(Channel_Type == "Single-Cell") %>%
  dplyr::select(SampleID, Index, Intensity) %>%
  pivot_wider(values_from = "Intensity",
              names_from = "SampleID") %>%
  column_to_rownames(var = "Index") %>%
  as.matrix()

metadata <- data.frame(SampleID = colnames(df2)) %>%
  mutate(Method = case_when(grepl("RTLS", SampleID) ~ "RTLS",
                                     grepl("STD_MS2", SampleID) ~ "STD_MS2"))


sce <- SingleCellExperiment(assays = List(df2),)
sim <- QFeatures(experiments = List(Assay1 = sce))

metadata <- metadata %>%
  arrange(SampleID,sim@colData@rownames)

sim$SampleType <- metadata$Method


data_completedness <- reportMissingValues(sim, "Assay1", by = sim$SampleType) %>%
  mutate(Completeness = Completeness*100) %>%
  mutate(Type = rownames(.))
dt1 <- data_completedness %>%
  filter(Type == "RTLS")
dt2 <- data_completedness %>%
  filter(Type == "STD_MS2")

csc <- cumulativeSensitivityCurve(
  sim, "Assay1", by = sim$SampleType,
) %>%
  mutate(Type = by)

csc %>%
  ggplot()+
  aes(x = SampleSize, y = Sensitivity, colour = Type) +
  geom_point(size =0.5)+
  # geom_hline(yintercept = dt1$LocalSensitivityMean,
  #            linetype = "dashed",
  #           color = "#F8766D")+
  geom_hline(yintercept = dt1$TotalSensitivity,
             color = "#F8766D")+
  # geom_hline(yintercept = dt2$LocalSensitivityMean,
  #            linetype = "dashed",
  #           color = "#00BFC4")+
  geom_hline(yintercept = dt2$TotalSensitivity,
             color = "#00BFC4")+
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black'),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  ylab("Sensitivity")+
  xlab("Sample Size")+
  scale_y_continuous(limits = c(0,4500))
 # ggsave("Figure_2a.png", width = 3.3, height = 2.4)


data_completedness %>%
  ggplot()+
  aes(x = Completeness, y = LocalSensitivityMean, color = Type,
      group = Type)+
  geom_point(size =6,
             shape = 5)+
  theme_bw()+
  scale_x_continuous(limits = c(10,55))+
  scale_y_continuous(limits = c(0,2300))+
  geom_pointrange(aes(ymin=LocalSensitivityMean-LocalSensitivitySd,
                      ymax=LocalSensitivityMean+LocalSensitivitySd),
                  position=position_dodge(0.05),
                  size = 0.1)+
  theme_bw(base_size = 8) +
  theme(legend.position = "none",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black'),
        axis.text.x=element_text(angle = 0, vjust = 0.5, hjust= 0.5, color='black'),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  ylab("Local Sensitivty")+
  xlab("Data completeness (%)")
#  ggsave("Figure_2b.png",
   #      height = 2,
   #      width = 3.3)


################################

df_long <- read_tsv("./Data/abundance_peptide_None.tsv") %>%
  pivot_longer(-Index,
               names_to = "SampleID",
               values_to = "Intensity") %>%
  filter(grepl("pg_", SampleID)) 

# Summarize per sample
sample_summary <- df_long %>%
  group_by(SampleID) %>%
  summarize(
    Mean = mean(Intensity, na.rm = TRUE),
    Count = sum(!is.na(Intensity))
  ) %>%
  ungroup()

# Extract Method (RTLS or Standard)
sample_summary <- sample_summary %>%
  mutate(
    Group_Input_Level = str_extract(SampleID, "1000pg|250pg|100pg|50pg"),
    Channel_Type = case_when(
    #  str_detect(SampleID, "(^|[^A-Z0-9])127C([^A-Z0-9]|$)") | str_detect(SampleID, "(^|[^A-Z0-9])134ND([^A-Z0-9]|$)") ~ "Isotope",
      str_detect(SampleID, "(^|[^A-Z0-9])135CD([^A-Z0-9]|$)") |
        str_detect(SampleID, "(^|[^A-Z0-9])134C([^A-Z0-9]|$)") |
        str_detect(SampleID, "(^|[^A-Z0-9])135N([^A-Z0-9]|$)") ~ "Blank",
      str_detect(SampleID, "RefInt|RefDInt") ~ "Bridge",
      TRUE ~ "Single-Cell"
    )) 
  


sample_summary$Group_Input_Level <- factor(
  sample_summary$Group_Input_Level,
  levels = c("1000pg", "250pg", "100pg", "50pg")
)


# --- PLOT 2: Count of Non-NA Values per Sample ---
p1 <- ggplot(sample_summary, aes(x = Group_Input_Level, y = Count, fill = Channel_Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_jitter(
    aes(color = Channel_Type),
    #width = 0.2,
    size = 0.5,
    alpha = 0.8,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)
  ) +
  scale_color_manual(values = c("#DC267F",
                               "#648FFF",
                               "#FFB000"))+
  scale_fill_manual(values = c("#DC267F",
                               "#648FFF",
                                "#FFB000"))+
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom",
        panel.background = element_rect(fill= 'white'),
        axis.text.y=element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  labs(
    # title = "Number of Peptides by Channel Type",
    x = "Bridge Input Level",
    y = "Peptide Count (n)",
    fill = "Channel Type",
    color = "Channel Type"
  ) +
  theme(
    legend.position = "bottom"
  )

#p1+ggsave("Figure_2c.png",width = 3.3, height = 2.4)
       
