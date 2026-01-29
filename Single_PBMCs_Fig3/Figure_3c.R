library(tidyverse)
library(rstudioapi)
#######################

# Get the directory of the current source file
script_dir <- dirname(getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)

### Import data/metadata

load("./Data/MASIC_TMT_ReportIon_10ppm_Tolerance.RData")


meta1 <-read.csv(file = "./Metadata/PBMC_ReporterIon_Metadata.csv") %>%
  mutate(Channel = sprintf("%.4f",Channel))

############################


masic_data2 <- masic_data %>%
  filter(InterferenceScore >= 0.5) %>%
  select(-contains("SignalToNoise")) %>%
  pivot_longer(-c(Dataset,ScanNumber,InterferenceScore),
               values_to = "Intensity",
               names_to = "Channel") %>%
  mutate(Channel = str_sub(Channel,-8)) %>%
  group_by(Dataset,ScanNumber) %>%
  mutate(Sum1 = sum(Intensity)) %>%
  ungroup() %>%
  filter(Sum1 > 0) %>%
  left_join(., meta1) %>%
  filter(!(ScanNumber %in% ScanNumber[Type == "Bridge" & Intensity == 0])) %>%
  select(-Sum1, -InterferenceScore) 


###############

summed_intensities <- masic_data2 %>%
  filter(Intensity != 0) %>%
  group_by(Dataset, Channel) %>%
  mutate(
    Q1 = quantile(Intensity, 0.25, na.rm = TRUE),
    Q3 = quantile(Intensity, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1
  ) %>%
  filter(
    Intensity >= Q1 - 1.5 * IQR,
    Intensity <= Q3 + 1.5 * IQR
  ) %>%
  group_by(Dataset, Channel, Type) %>%
  summarize(SumInt = log10(sum(Intensity, na.rm = T)))  %>%
  ungroup()# %>%
# group_by(Type) %>%
 #summarize(median(SumInt))


summed_intensities %>%
  ggplot()+
  aes(x = SumInt,  fill = factor(Type,
                                 levels = c("Single-cell",
                                            "Bridge",
                                            "Blank"
                                 )))+
  geom_histogram( alpha = 0.6, color  = "black", position = "identity",
                  bins = 100)+
  #  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5))+
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(title = "Channel Type"), color = "none")+
  geom_vline(xintercept = 7.08, linetype  = "dashed",
             linewidth = 0.7)+
  geom_vline(xintercept = 6.06, linetype  = "dashed",
             linewidth = 0.7)+
  geom_vline(xintercept = 3.48, linetype  = "dashed",
             linewidth = 0.7)+
  scale_fill_manual(values = c("#FFB000",
                               "#648FFF",
                               "#DC267F"))+
  xlab("log10(Sum(Reporter Ion Intensity))")+
  ylab("Cells (n)") +
  ggsave("Figure_3c.png", width = 6, height = 4)
