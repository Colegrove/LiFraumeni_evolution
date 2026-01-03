
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_chemo_clock/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
# filtered matrix
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"


df <- read_delim(sigPath, delim = "\t")
data_long <- df %>%
  pivot_longer(cols = -Samples, names_to = "signature", values_to = "count") %>%
  filter(count > 0)
data_long %>% print(n=Inf)

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
  
  ## unknown aetiology (or not grouped within github)
  "SBS10c" = "#505050", "SBS10d" = "#4a4a4a", "SBS11" = "#aa4499", "SBS12" = "#262626", 
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



# Reverse order for plotting while keeping legend in original order
data_long$signature <- factor(data_long$signature, levels = rev(names(signature_colors)))

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


data_long %>% print(n=Inf)
#manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")
subjects <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Patient:", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
subjects <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Buffy coat", "Whole blood", "Plasma" , "Bone marrow" ,"Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")

data_long$x_axis_label <- factor(data_long$x_axis_label, levels = subjects)
data_long <- data_long %>% filter(tissue == "PBMC" | tissue == "Buffy coat")

#### Cosine Similarity stats
sigStats_all_panel <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt"
#sigStats_all_panel <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_chemo_clock/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt"
sigStats_all_panel_filtered <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt"

all_panel_stats <- read_delim(sigStats_all_panel_filtered, delim = "\t")
all_panel_stats
data_long
data_long <- data_long %>%
  left_join(all_panel_stats %>% select("Sample Names", "Cosine Similarity"), by = c("Samples" = "Sample Names"))
data_long


p1 <- ggplot(data_long, aes(x = x_axis_label, y = count, fill = signature)) +
  geom_bar(stat= "identity", position = "fill") +
  #geom_bar(stat = "identity") +
  labs(y = "Count") +
  theme_minimal() +
  #ylim(0,800) +
  geom_text(
    data = data_long,
    aes(x = x_axis_label, y = 1.02, label = round(`Cosine Similarity`, 3)),
    inherit.aes = FALSE
  ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  #scale_fill_manual(values = signature_colors, breaks = names(signature_colors))
  scale_fill_manual(values = signature_colors,breaks = names(signature_colors),labels = legend_labels)

show(p1)

path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/out/plots/SBS96_signatures_blood_allPanels.png"
ggsave(path, p1, width = 12, height = 6)


##### Add age and LFS information
#subjects <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Buffy coat", "Whole blood", "Plasma" , "Bone marrow" ,"Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")

subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
ages <- c("20-30", "20-30", "20-30", "20-30", "30s", "39", "61", "60s", "69", "70s", "70s")
LFS <- c("-", "-", "-", "-", "+", "+", "-","-", "+","-", "-")
CTx <- c("-", "-", "-", "-", "+", "-", "-","-", "-","-", "+")

age_lf <- data.frame(
  x_axis_var = subjects, 
  CTx = CTx,
  LFS = LFS,
  age = ages
)

age_lf_long <- age_lf %>%
  pivot_longer(cols = c(age, LFS, CTx), names_to = "variable", values_to = "value") %>%
  mutate(y = case_when(
    variable == "age" ~ 2,
    variable == "LFS" ~ 1.5,
    variable == "CTx" ~ 1, 
  )) %>% 
  mutate(x_axis_var = factor(x_axis_var, levels = unique(x_axis_var)),
         x_axis_var = factor(x_axis_var, levels = unique(x_axis_var))  # preserve original order
  )


legend <- cowplot::get_legend(
  p1 + theme(legend.position = "right")
)
legend_aligned <- cowplot::plot_grid(
  legend,
  NULL,  # spacer to push legend upward
  ncol = 1,
  rel_heights = c(1, .1)  # stretch legend space
)

p1_nolegend <- p1 +
  theme(legend.position = "none")

combined_plot <- cowplot::plot_grid(
  p1_nolegend + theme(plot.margin = margin(0, 0, 2, 0, "pt")),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE, 
  align = "v", 
  axis = "lr"
)

combined_plot_legend <- cowplot::plot_grid(
  combined_plot, 
  legend_aligned,
  ncol = 2, 
  rel_widths = c(1, 0.25), 
  align = "h",
  axis = "tb"
  
)

show(combined_plot_legend)



path1 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_blood_allPanels_filtered/out/plots/SBS96_signatures_blood_allPanels_clinical_info.png"
ggsave(path1, combined_plot_legend, width = 12, height = 6, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)


#################################################################
### compare mut panel only to mut + chip + tp53
#################################################################
sigPath_mut_panel_only <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
#sigPath_mut_panel_only <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/sigProfilerAssignment_ouput/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"

df_mut_panel <- read_delim(sigPath_mut_panel_only, delim = "\t")
data_long_mut_panel <- df_mut_panel %>%
  pivot_longer(cols = -Samples, names_to = "signature", values_to = "count") %>%
  filter(count > 0)
data_long_mut_panel %>% print(n=Inf)

data_long_mut_panel$signature <- factor(data_long_mut_panel$signature, levels = rev(names(signature_colors)))
# Assign sample information
sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)
sampleID <- sampleID %>%
  mutate(Samples = sample) %>%
  dplyr::select(-sample)
data_long_mut_panel <- left_join(data_long_mut_panel, sampleID, by = "Samples")
data_long_mut_panel <- data_long_mut_panel %>%
  mutate(x_axis_label = if_else(subject != "Patient:", subject, tissue))
subjects <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Buffy coat", "Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")
data_long_mut_panel$x_axis_label <- factor(data_long_mut_panel$x_axis_label, levels = subjects)
data_long_mut_panel <- data_long_mut_panel %>% filter(tissue == "PBMC" | tissue == "Buffy coat")

#### Cosine Similarity stats
sigStats_mut_panel <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/sigProfilerAssignment_ouput_3.2/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt"
mut_panel_stats <- read_delim(sigStats_mut_panel, delim = "\t")
mut_panel_stats
data_long_mut_panel
data_long_mut_panel <- data_long_mut_panel %>%
  left_join(mut_panel_stats %>% select("Sample Names", "Cosine Similarity"), by = c("Samples" = "Sample Names"))
data_long_mut_panel

p1 <- ggplot(data_long_mut_panel, aes(x = x_axis_label, y = count, fill = signature)) +
  geom_bar(stat = "identity", position = "fill") +
  #geom_bar(stat = "identity") +
  labs(y = "Count") +
  theme_minimal() +
  #ylim(0,800) +
  geom_text(
    data = data_long_mut_panel,
    aes(x = x_axis_label, y = 1.02, label = round(`Cosine Similarity`, 3)),
    inherit.aes = FALSE
  ) +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors))

show(p1)



#### compare mut panel to all panel
data_long <- data_long %>%
  mutate(panel = "ALL")
data_long_mut_panel <- data_long_mut_panel %>%
  mutate(panel = "MUT")
all_data_long <- data_long %>%
  rbind(data_long_mut_panel)

p2 <- ggplot(all_data_long, aes(x = panel, y = count, fill = signature)) +
  geom_bar(
    stat = "identity",
    position = "fill", 
    width = 1
  ) +
  labs(y = "Count") +
  theme_minimal() +
  geom_text(
    data = all_data_long,
    aes(x = panel, y = 1.06, label = round(`Cosine Similarity`, 3)),
    inherit.aes = FALSE
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) +
  facet_wrap(~x_axis_label, nrow=2) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors))

show(p2)
ggsave("results/mutational_signature_comparison.png", p2, width = 10, height = 5, units = "in", dpi = 300)
#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)


#################################################################
### SBS 1/SBS5 ratio
#################################################################
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





#### heatmap of common mutational signatures

signature_order <- data_long %>%
  group_by(signature) %>%
  summarise(Sample_Count = dplyr::n_distinct(subject)) %>%
  arrange(Sample_Count) %>%
  pull(signature)
signature_order


data_long$signature <- factor(data_long$signature, levels = signature_order)

manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")
manual_subject_order <- c("UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "Buffy coat", "Whole blood", "Plasma" , "Bone marrow" ,"Family member A:", "Family member B:", "UW volunteer 5:", "Family member C:", "UW volunteer 6:", "UW volunteer 7:")

data_long$subject <- factor(data_long$subject, levels = manual_subject_order)

signature_heatmap_family <- ggplot(data_long, aes(x = x_axis_label, y = signature, fill = count)) +
  geom_tile() +  # Create the heatmap
  scale_fill_gradient(low = "white", high = "red") +  # Color scale from white to red
  geom_text(aes(label = count), color = "black", size = 4) +  # Add counts in the middle
  theme_minimal() +
  labs(x = "Sample", y = "Signature") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
show(signature_heatmap_family)
path3 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signature_family_control_heatmap.png"
ggsave(path3, signature_heatmap_family, width = 12, height = 6)
