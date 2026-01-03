
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures_family_controls/sigProfilerAssignment_ouput/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"

df <- read_delim(sigPath, delim = "\t")
df
data_long <- df %>%
  pivot_longer(cols = -Samples, names_to = "signature", values_to = "count") %>%
  filter(count > 0)
data_long %>% print(n=Inf)


signature_colors <- c(
  "SBS1" = "#acf2d0", "SBS5" = "#63d69e", "SBS2" = "#f8b6b3", "SBS13" = "#f17fb2",
  "SBS3" = "#c4abc4", "SBS4" = "#bcf2f5", "SBS7a" = "#b5d7f5", "SBS7b" = "#9ecef7",
  "SBS7c" = "#84bdf0", "SBS7d" = "#6cb2f0", "SBS8" = "#dfc4f5", "SBS9" = "#ebf5bc",
  "SBS10a" = "#f2aeae", "SBS10b" = "#f08080", "SBS17a" = "#d9f7b0", "SBS17b" = "#8cc63f",
  "SBS40" = "#c4c4f5", "SBS6" = "#faf1dc", "SBS14" = "#faecca", "SBS15" = "#fcebc2",
  "SBS20" = "#fae4af", "SBS21" = "#fae1a5", "SBS26" = "#fcde97", "SBS44" = "#fad682",
  "SBS27" = "#C8C8C8", "SBS43" = "#C0C0C0", "SBS45" = "#BEBEBE", "SBS46" = "#B8B8B8",
  "SBS47" = "#B0B0B0", "SBS48" = "#A9A9A9", "SBS49" = "#A8A8A8", "SBS50" = "#A0A0A0",
  "SBS51" = "#989898", "SBS52" = "#909090", "SBS53" = "#888888", "SBS54" = "#808080",
  "SBS55" = "#787878", "SBS56" = "#707070", "SBS57" = "#696969", "SBS58" = "#686868",
  "SBS59" = "#606060", "SBS60" = "#585858"
)

extra_colors <- c(
  "pink", "orange", "purple", "brown", "red", 
  "green", "cyan", "deeppink", "orangered", "blueviolet", "chocolate", 
  "darkgreen", "dodgerblue", "mediumvioletred", "salmon", "magenta", "sandybrown", 
  "forestgreen", "royalblue", "orchid", "indigo", "darkseagreen", "blue", 
  "palevioletred", "darkslateblue", "olivedrab", "cyan", "hotpink"
)

# Identify missing signatures
missing_signatures <- setdiff(unique(data_long$signature), names(signature_colors))

# Assign extra colors to missing signatures
if (length(missing_signatures) > length(extra_colors)) {
  # Generate random colors if not enough defined ones
  random_colors <- sprintf("#%06X", sample(0x000000:0xFFFFFF, length(missing_signatures) - length(extra_colors), replace = FALSE))
  assigned_colors <- c(extra_colors, random_colors)
} else {
  assigned_colors <- extra_colors[1:length(missing_signatures)]
}

signature_colors <- c(signature_colors, setNames(assigned_colors, missing_signatures))
signature_colors
# Reverse order for plotting while keeping legend in original order
data_long$signature <- factor(data_long$signature, levels = rev(names(signature_colors)))

sampleID_map <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sampleID <- read_delim(sampleID_map, delim = "\t")
sampleID <- sampleID %>%
  dplyr::select(-sampleID, -tissue)
sampleID[] <- lapply(sampleID, function(x) if(is.character(x)) trimws(x) else x)
sampleID <- sampleID %>%
  mutate(Samples = sample) %>%
  dplyr::select(-sample)
sampleID
data_long
data_long <- left_join(data_long, sampleID, by = "Samples")
data_long

manual_subject_order <- c("Patient:", "Family member A:", "Family member C:", "Family member B:", "UW volunteer 1:", "UW volunteer 2:", "UW volunteer 3:", "UW volunteer 4:", "UW volunteer 5:", "UW volunteer 6:", "UW volunteer 7:")

data_long$subject <- factor(data_long$subject, levels = manual_subject_order)


p1 <- ggplot(data_long, aes(x = subject, y = count, fill = signature)) +
  geom_bar(stat = "identity") +
  labs(y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(), 
        legend.title=element_blank()) +
  scale_fill_manual(values = signature_colors, breaks = names(signature_colors))

show(p1)

path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signatures_family_controls.png"
ggsave(path, p1, width = 10, height = 6)


### SBS 1/SBS5 ratio

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
data_long$subject <- factor(data_long$subject, levels = manual_subject_order)

signature_heatmap_family <- ggplot(data_long, aes(x = subject, y = signature, fill = count)) +
  geom_tile() +  # Create the heatmap
  scale_fill_gradient(low = "white", high = "red") +  # Color scale from white to red
  geom_text(aes(label = count), color = "black", size = 4) +  # Add counts in the middle
  theme_minimal() +
  labs(x = "Sample", y = "Signature") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
show(signature_heatmap_family)
path3 <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/out/plots/SBS96_signature_family_control_heatmap.png"
ggsave(path3, signature_heatmap_family, width = 12, height = 6)
