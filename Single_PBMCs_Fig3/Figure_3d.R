###################
library(data.table)
library(tidyverse)
library(rstudioapi)
#######################

# Get the directory of the current source file
script_dir <- dirname(getActiveDocumentContext()$path)

# Set the working directory to the script's directory
setwd(script_dir)

### Import data/metadata

load("./Data/MASIC_TMT_ReportIon_10ppm_Tolerance.RData")

reference <- read.csv(file = "./Metadata/Channels_Meta.csv",
                      colClasses = "character") %>%
  mutate(Channel = as.character(Channel))

meta1 <-read.csv(file = "./Metadata/PBMC_ReporterIon_Metadata.csv") %>%
  mutate(Channel = sprintf("%.4f",Channel))


#### Processing steps

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


### -------------------------------------------------------
### 0. Convert df to data.table
### -------------------------------------------------------

dt <- as.data.table(masic_data2)
dt <- dt[Type %in% c("Single-cell", "Blank", "Bridge")]

### -------------------------------------------------------
### 1. Per-scan sample-positive indicator
### -------------------------------------------------------

dt[, SamplePositive := (Type %in% c("Single-cell","Blank") & Intensity > 0)]

### -------------------------------------------------------
### 2. Per-scan total bridge intensity by Dataset × Scan × Deuterated
### -------------------------------------------------------

bridge_per_scan <- dt[Type == "Bridge",
                      .(BridgeScanSum = sum(Intensity)),
                      by = .(Dataset, ScanNumber, Deuterated)]

dt <- merge(dt, bridge_per_scan,
            by = c("Dataset", "ScanNumber", "Deuterated"),
            all.x = TRUE)

### -------------------------------------------------------
### 3. Assign bridge intensity to sample rows for each scan
### -------------------------------------------------------

dt[, BridgeAssigned :=
     fifelse(Type %in% c("Single-cell","Blank"),
             fifelse(SamplePositive, BridgeScanSum, NA_real_),   # 0 → NA
             Intensity)                                           # bridge row keeps raw intensity
]

### -------------------------------------------------------
### 4. Assign sample intensity per-scan; convert 0 → NA
### -------------------------------------------------------

dt[, SampleIntensityUsed :=
     fifelse(Type %in% c("Single-cell","Blank"),
             fifelse(Intensity > 0, Intensity, NA_real_),         # 0 → NA
             NA_real_)
]

### -------------------------------------------------------
### 5. Robust median summarization (ignoring NA values)
### -------------------------------------------------------

result <- dt[, {
  
  sample_vals <- SampleIntensityUsed[Type %in% c("Single-cell","Blank")]
  bridge_vals <- BridgeAssigned[Type != "Bridge"]   # per-scan bridge intensity assigned to this channel
  
  .(
    SampleSummedIntensity = median(sample_vals, na.rm = TRUE),
    BridgeSummedIntensity = median(bridge_vals, na.rm = TRUE)
  )
  
}, by = .(Dataset, Channel, Deuterated)]

setorder(result, Dataset, Channel)

### Return final table

result <- result %>%
  left_join(., reference) %>%
  filter(grepl("SC", SingleCell_Channel)) %>%
  mutate(Ratio = BridgeSummedIntensity / SampleSummedIntensity) %>%
  mutate(Estimated_pgcell = 150/(Ratio)) 
  
result %>%
  ggplot()+
  aes(x = Estimated_pgcell)+
  geom_histogram( alpha = 0.6, color  = "black", position = "identity",
                  bins = 100)+
  theme_bw(base_size = 12) +
  geom_vline(xintercept = median(result$Estimated_pgcell), linetype  = "dashed",
             color = "black",
             size = 1)+
  annotate("text", x = 27, y = 100, label = "Median = 14.2 pg", color = "black", size = 5)+
  scale_x_continuous(limits = c(0, 75), oob = scales::squish)+
  labs(
    x = "Estimated peptide (pg)",
    y = "Cells (n)"
  )+
  ggsave("Figure_3d.png", width = 6, height = 4)

  
pg_ref <- result %>%
  mutate(SampleID = paste(Dataset, TMT, sep = "_")) %>%
  select(SampleID, Estimated_pgcell)

save(pg_ref, file = "PBMC_pg_peptide_estimate.RData")
