

#SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/sigProfilerAssignment_ouput/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
## filtered maf
SigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"


df <- read_delim(SigPath, delim = "\t")
data_long <- df %>%
  pivot_longer(cols = -Samples, names_to = "signature", values_to = "count") %>%
  filter(count > 0)
data_long %>% print(n=Inf)

sample <- "DNA03974_S11.1"
df1 <- read_delim(paste0("/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_allSamples_filtered/",sample,".filtered.maf"), delim = "\t", skip=1)
df1
df2 <- maf_masked_coding %>% filter(Samp == sample)
df2




signature_colors <- c(
  ## original colors used in sigProfiler plotting scheme
  "SBS1" = "#acf2d0", "SBS5" = "#63d69e", "SBS2" = "#f8b6b3", "SBS13" = "#f17fb2",
  "SBS3" = "#c4abc4", "SBS4" = "#bcf2f5", "SBS7a" = "#b5d7f5", "SBS7b" = "#9ecef7",
  "SBS7c" = "#84bdf0", "SBS7d" = "#6cb2f0", "SBS8" = "#dfc4f5", "SBS9" = "#ebf5bc",
  "SBS10a" = "#f2aeae", "SBS10b" = "#f08080", "SBS17a" = "#d9f7b0", "SBS17b" = "#8cc63f",
  "SBS40" = "#c4c4f5", 
  
  ## added colors with known aetiology that are present in our dataset
  ## associated with mismatch repair (sigprofiler color)
  "SBS6" = "#faf1dc", "SBS14" = "#faecca", "SBS15" = "#fcebc2",
  "SBS20" = "#fae4af", "SBS21" = "#fae1a5", "SBS26" = "#fcde97", "SBS44" = "#fad682",
  
  ## associated with sequencing artifact (sigprofiler color)
  "SBS27" = "#C8C8C8", "SBS43" = "#C0C0C0", "SBS45" = "#BEBEBE", "SBS46" = "#B8B8B8",
  "SBS47" = "#B0B0B0", "SBS48" = "#A9A9A9", "SBS49" = "#A8A8A8", "SBS50" = "#A0A0A0",
  "SBS51" = "#989898", "SBS52" = "#909090", "SBS53" = "#888888", "SBS54" = "#808080",
  "SBS55" = "#787878", "SBS56" = "#707070", "SBS57" = "#696969", "SBS58" = "#686868",
  "SBS59" = "#606060", "SBS60" = "#585858", 
  
  ## unknown aetiology (or not present in our dataset)
  "SBS10c" = "#505050", "SBS10d" = "#4a4a4a", "SBS11" = "#444444", "SBS12" = "#262626", 
  "SBS16" = "#383838", "SBS18" = "#323232", "SBS19" = "#2c2c2c", "SBS22" = "#262626",
  "SBS23" = "#202020", "SBS24" = "#1a1a1a", "SBS25" = "#141414", "SBS28" = "#3cc63f",
  "SBS29" = "#999933", "SBS30" = "#ff7f00", "SBS31" = "#aa4499", "SBS32" = "#282828", 
  "SBS33" = "#222222",  "SBS34" = "#1c1c1c", "SBS35" = "#884455", "SBS36" = "#101010", 
  "SBS37" = "#4c4c4c", "SBS38" = "#bcf2f5","SBS39" = "#404040", "SBS40c" = "#3a3a3a", 
  "SBS41" = "#3a3a3a", "SBS42" = "#343434","SBS84" = "#282828", "SBS85" = "#282828", 
  "SBS86" = "red", "SBS87" = "pink", "SBS88" = "#161616", "SBS89" = "#101010",
  "SBS90" = "#0a0a0a", "SBS91" = "#040404", "SBS92" = "#303030", "SBS93" = "#2a2a2a",
  "SBS94" = "#242424", "SBS95" = "#1e1e1e", "SBS96" = "#181818", "SBS97" = "#121212",
  "SBS98" = "#0c0c0c", "SBS99" = "#060606", "SBS100" = "#020202"
  )

## add known etiologies
data_long <- data_long %>%
  mutate(sig_etiology = case_when(
    signature %in% paste0("SBS", c("6","14","15","20","21","26","44")) ~ "MMR_deficiency",
    signature %in% paste0("SBS", c("10a","10b","10c","10d","28"))     ~ "POL_deficiency",
    signature %in% "SBS3"                                              ~ "HR_deficiency",
    signature %in% paste0("SBS", c("30","36"))                         ~ "BER_deficiency",
    #signature %in% paste0("SBS", c("11","25","31","35","86","87","90","99")) ~ "Chemotherapy",
    signature %in% "SBS32"                                             ~ "Immunosuppressants",
    signature %in% paste0("SBS", c("11","25","31","32","35","86","87","90","99")) ~ "Treatment",
    signature %in% paste0("SBS", c("2","13"))                          ~ "APOBEC",
    signature %in% paste0("SBS", c("4","29","92"))                      ~ "Tobacco",
    signature %in% paste0("SBS", c("7a","7b","7c","7d","38"))           ~ "UV",
    signature %in% paste0("SBS", c("22a","22b"))                       ~ "Aristolochic Acid",
    signature %in% "SBS88"                                             ~ "Colibactin",
    signature %in% paste0("SBS", c("27","43","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","95")) 
    ~ "Artifact",
    signature %in% paste0("SBS", c("9","84","85"))                      ~ "Lymphoid",
    signature %in% "SBS8"                                              ~ "ROS",
    signature %in% paste0("SBS", c("1","5","40"))                      ~ "Clock-like",
    TRUE                                                               ~ "Unknown"
  ))

etiology_order <- c("Unknown", "Artifact","MMR_deficiency","POL_deficiency","ROS","HR_deficiency",
                    "BER_deficiency","Immunosuppressants","Treatment",
                    "APOBEC","Tobacco","UV","Aristolochic Acid","Colibactin",
                    "Lymphoid","Clock-like")
data_long <- data_long %>% 
  mutate(sig_etiology = factor(sig_etiology, levels = etiology_order))
data_long %>% filter(sig_etiology == "Unknown") %>% print(n=Inf)
etiology_colors <- c(
  "Clock-like"            = "#63d69e",
  "MMR_deficiency"        = "#1f78b4",
  "POL_deficiency"        = "#33a02c",
  "HR_deficiency"         = "#c4abc4",
  "BER_deficiency"        = "#ff7f00",
  "Immunosuppressants"    = "#6a3d9a",
  "Treatment"             = "#b15928",
  "APOBEC"                = "#f8b6b3",
  "Tobacco"               = "#882255",
  "UV"                    = "#84bdf0",
  "Aristolochic Acid"     = "#cab2d6",
  "Colibactin"            = "#ffff99",
  "ROS"                   = "#dfc4f5",
  "Lymphoid"              = "#1b9e77",
  "Artifact"              = "#DDDDDD",
  "Unknown"               = "#111111"
)

## add etiology to SBS signatures
sig_etio_map <- data_long %>% 
  distinct(signature, sig_etiology) %>%
  { setNames(.$sig_etiology, .$signature) }
sig_etio_map
legend_labels <- paste0(
  names(signature_colors),          # e.g. "SBS6"
  " (",
  sig_etio_map[names(signature_colors)],  # e.g. "MMR_deficiency"
  ")"
)
names(legend_labels) <- names(signature_colors)


# Assign sample information
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)
sampleID <- sampleID %>%
  mutate(Samples = sample) %>%
  dplyr::select(-sample)
data_long <- left_join(data_long, sampleID, by = "Samples")
data_long <- data_long %>%
  mutate(x_axis_label = if_else(subject != "Patient:", subject, tissue))

# Reverse order for plotting while keeping legend in original order
tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Liver", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")
data_long$signature <- factor(data_long$signature, levels = rev(names(signature_colors)))
data_long <- data_long %>%
  mutate(x_axis_label = factor(x_axis_label, levels = tissue_order))

#### Cosine Similarity stats
sigStats_all_panel <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt"

all_panel_stats <- read_delim(sigStats_all_panel, delim = "\t")
data_long <- data_long %>%
  left_join(all_panel_stats %>% select("Sample Names", "Cosine Similarity"), by = c("Samples" = "Sample Names"))
data_long

custom_label <- function(x) {
  sapply(x, function(t) {
    if (t %in% c("Esophageal cancer 1",
                 "Esophageal cancer 2",
                 "Mediastinal metastasis",
                 "Lung metastasis",
                 "Liver metastasis 1",
                 "Liver metastasis 2")){
      paste0("<span style='color:red;'>", t, "</span>")
    } else {
      t
    }
  })
}
data_long
p1 <- ggplot(data_long, aes(x = x_axis_label, y = count, fill = signature)) +
#p1 <- ggplot(data_long, aes(x = x_axis_label, y = count, fill = sig_etiology)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_bar(stat = "identity") +
  labs(y = "Count") +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  geom_text(
    data = data_long,
    aes(x = x_axis_label, y = 1.02, label = round(`Cosine Similarity`, 3)),
    inherit.aes = FALSE
  ) +
  theme(axis.text.x = element_markdown(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  #scale_fill_manual(values = etiology_colors, breaks = names(etiology_colors))
  #scale_fill_manual(values = signature_colors, breaks = names(signature_colors))
  scale_fill_manual(values = signature_colors,breaks = names(signature_colors),labels = legend_labels)

show(p1)
path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/out/plots/SBS96_signatures_patient_sigProfiler_3.2.png"
#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_sigProfiler_3.4.png"
ggsave(path, p1, width = 12, height = 6)

path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_patient_tissues_filtered/out/plots/SBS96_signatures_patient_sigProfiler_3.2_grouped.png"
ggsave(path, p1, width = 12, height = 6)
######### MuSiCaL
musicalSigPath_0015 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_0015.csv"
musicalSigPath_002 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_002.csv"
musicalSigPath_005 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_005.csv"
musicalSigPath_01 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_01.csv"
musicalSigPath_015 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_015.csv"

musicalSigPath_02 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_02.csv"
musicalSigPath_025 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_025.csv"
musicalSigPath_03 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_03.csv"
musicalSigPath_035 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_035.csv"
musicalSigPath_04 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_04.csv"

musicalSigPath_05 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_05.csv"
musicalSigPath_1 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_1.csv"
musicalSigPath_2 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_2.csv"
musicalSigPath_3 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_3.csv"
musicalSigPath_5 <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-03-31-musical_signatures/SBS_musical_5.csv"


import_musical_file <- function(file_path, threshold_value) {
  read_csv(file_path) %>%
    rename(signature = 1) %>%
    pivot_longer(
      cols = -signature,
      names_to = "Samples",
      values_to = "count"
    ) %>%
    filter(count >0) %>%
    mutate(threshold = threshold_value, 
    Samples = factor(Samples, levels = tissue_order), 
    signature = factor(signature, levels = rev(names(signature_colors))))
}

df_musical_0015 <- import_musical_file(musicalSigPath_0015, "0.0015")
df_musical_002 <- import_musical_file(musicalSigPath_002, "0.002")
df_musical_005  <- import_musical_file(musicalSigPath_005, "0.005")
df_musical_01  <- import_musical_file(musicalSigPath_01, "0.01")
df_musical_015   <- import_musical_file(musicalSigPath_015, "0.015")

df_musical_02 <- import_musical_file(musicalSigPath_02, "0.02")
df_musical_025 <- import_musical_file(musicalSigPath_025, "0.025")
df_musical_03  <- import_musical_file(musicalSigPath_03, "0.03")
df_musical_035  <- import_musical_file(musicalSigPath_035, "0.035")
df_musical_04   <- import_musical_file(musicalSigPath_04, "0.04")

df_musical_05 <- import_musical_file(musicalSigPath_05, "0.05")
df_musical_1 <- import_musical_file(musicalSigPath_1, "0.1")
df_musical_2  <- import_musical_file(musicalSigPath_2, "0.2")
df_musical_3  <- import_musical_file(musicalSigPath_3, "0.3")
df_musical_5   <- import_musical_file(musicalSigPath_5, "0.5")


data_long_comparison <- data_long %>% mutate(threshold = "sigProfiler")

df_musical_normal_range <- dplyr::bind_rows(
  data_long_comparison,
  df_musical_0015,
  df_musical_002,
  df_musical_005,
  df_musical_01,
  df_musical_015 
)

df_musical_high_range <- dplyr::bind_rows(
  data_long_comparison,
  df_musical_02,
  df_musical_025,
  df_musical_03,
  df_musical_035,
  df_musical_04 
)

df_musical_wide_range <- dplyr::bind_rows(
  data_long_comparison,
  df_musical_05,
  df_musical_1,
  df_musical_2,
  df_musical_3,
  df_musical_5 
)


p1_musical <- ggplot(df_musical_high_range, aes(x = Samples, y = count, fill = signature)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_bar(stat = "identity") +
  labs(y = "Count") +
  #theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  theme(axis.text.x = element_markdown(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors)) +
  #facet_wrap(~factor(threshold, c("sigProfiler", "0.0015", "0.002", "0.005", "0.01", "0.015")), nrow = 1)
  facet_wrap(~factor(threshold, c("sigProfiler", "0.02", "0.025", "0.03", "0.035", "0.04")), nrow = 1)
  #facet_wrap(~factor(threshold, c("sigProfiler", "0.05", "0.1", "0.2", "0.3", "0.5")), nrow = 1)

show(p1_musical)



#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_normal_range_COSMICv3.2.png"
#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_high_range_COSMICv3.2.png"
#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_wide_range_COSMICv3.2.png"

#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_normal_range_MuSiCalv3.2.png"
path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_high_range_MuSiCalv3.2.png"
#path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_wide_range_MuSiCalv3.2.png"

ggsave(path, p1_musical, width = 12, height = 6)

#########################
### plot by sample
#########################


threshold_order <- c("sigProfiler", "0.0015", "0.002", "0.005", "0.01", "0.015")
df_musical_normal_range$threshold <- factor(df_musical_normal_range$threshold, levels = threshold_order)
threshold_order <- c("sigProfiler", "0.02", "0.025", "0.03", "0.035", "0.04")
df_musical_high_range$threshold <- factor(df_musical_high_range$threshold, levels = threshold_order)
threshold_order <- c("sigProfiler", "0.05", "0.1", "0.2", "0.3", "0.5")
df_musical_wide_range$threshold <- factor(df_musical_wide_range$threshold, levels = threshold_order)

p1_musical_sample <- ggplot(df_musical_high_range, aes(x = factor(threshold), y = count, fill = signature)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_bar(stat = "identity") +
  labs(y = "Count") +
  #theme_minimal() +
  #scale_x_discrete(labels = threshold_order) +
  theme(axis.text.x = element_markdown(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors)) +
  facet_wrap(~Samples)

show(p1_musical_sample)

path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_musical_comparison_high_range_COSMICv3.2_threshold.png"

ggsave(path, p1_musical_sample, width = 12, height = 6)



########################

data_long_cancer <- data_long %>%
  filter(Samples %in% c("Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis"))

p3 <- ggplot(data_long_cancer, aes(x = Samples, y = count, fill = signature)) +
  geom_bar(stat = "identity") +
  labs(y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors))

show(p3)
path3 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_patient_cancer.png"
ggsave(path3, p3, width = 6, height = 6)





##### SBS1/SBS5 ratio
sbs_ratio <- df %>% 
  mutate(ratio = SBS1/SBS5) %>%
  dplyr::select(Samples, ratio)

sbs_1_5_only <- data_long %>%
  filter(signature == "SBS1" | signature == "SBS5")
sbs_1_5_only
sbs_ratio <- sbs_1_5_only %>%
  left_join(sbs_ratio, by = "Samples")
sbs_ratio
ratio_labels <- sbs_ratio %>%
  #filter(signature == "SBS1" | signature == "SBS5") %>%
  group_by(Samples) %>%
  mutate(count = sum(count)) %>%
  filter(signature == "SBS1") %>%
  mutate(label = sprintf("%.2f", ratio))
ratio_labels

p2 <- ggplot(sbs_1_5_only, aes(x = Samples, y = count, fill = signature)) +
  geom_bar(stat = "identity") +
  labs(y = "Count", title="SBS1/SBS5") +
  geom_text(data = ratio_labels, aes(label=label, y = count + 5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors)

show(p2)

path2 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS_1_5_sigs.png"
ggsave(path2, p2, width = 10, height = 6)


data_long

#### heatmap of common mutational signatures in patient normal and cancer

signature_order <- data_long %>%
  group_by(signature) %>%
  summarise(Sample_Count = dplyr::n_distinct(Samples)) %>%
  arrange(Sample_Count) %>%
  pull(signature)
signature_order

data_long$signature <- factor(data_long$signature, levels = signature_order)

manual_tissue_order <- c("Whole blood", "Buffy coat", "Plasma", "Bone marrow", "Colon", "Skin", "Skin, non-sun-exposed", "Liver", "Esophagus 1", "Esophageal cancer 1", "Esophageal cancer 2", "Liver metastasis 1", "Liver metastasis 2", "Lung metastasis", "Mediastinal metastasis")
data_long$Samples <- factor(data_long$Samples, levels = manual_tissue_order)

signature_heatmap <- ggplot(data_long, aes(x = Samples, y = signature, fill = count)) +
  geom_tile() +  # Create the heatmap
  scale_fill_gradient(low = "white", high = "red") +  # Color scale from white to red
  geom_text(aes(label = count), color = "black", size = 4) +  # Add counts in the middle
  theme_minimal() +
  labs(x = "Sample", y = "Signature") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability

show(signature_heatmap)
path3 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signature_patient_heatmap.png"
ggsave(path3, signature_heatmap, width = 12, height = 6)




