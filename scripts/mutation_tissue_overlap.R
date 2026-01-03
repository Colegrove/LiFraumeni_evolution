
require(tidyverse)
require(stringr)


tissue_order <- c("Whole blood", 
                 "Buffy coat", 
                 "Plasma", 
                 "Bone marrow", 
                 "Buccal mucosa", 
                 "Thyroid", 
                 "Mainstem bronchus",
                 "Lung", 
                 "Esophagus 1", 
                 "Esophagus 2", 
                 "Gastric 1",
                 "Gastric 2",
                 "Cardiac muscle",
                 "Spleen",
                 "Liver",
                 "Colon",
                 "Omentum",
                 "Peritoneum",
                 "Renal",
                 "Testis",
                 "Skeletal muscle",
                 "Skin",
                 "Skin, non-sun-exposed",
                 "Mediastinal metastasis",
                 "Lung metastasis",
                 "Esophageal cancer 1",
                 "Esophageal cancer 2",
                 "Liver metastasis 1",
                 "Liver metastasis 2")

# MAF_table %>% filter(Tissue == "Skin, non-sun-exposed") %>% 
#   arrange(desc(VAF)) %>% 
#   filter(Hugo_Symbol == "TP53") %>% 
#   filter(Variant_Classification != "Intron") %>% 
#   #filter(prot.pos > 150 & prot.pos < 200) %>%
#   filter(Subject == "Patient") %>%
#   print(n=Inf, width = Inf)

# MAF_table %>% #filter(Tissue == "Skin") %>% 
#   arrange(desc(VAF)) %>% 
#   filter(Hugo_Symbol == "TP53") %>% 
#   filter(Variant_Classification != "Intron") %>% 
#   filter(prot.pos > 150 & prot.pos < 200) %>%
#   filter(Subject == "Patient") %>%
#   print(n=Inf, width = Inf)


df_filtered <- MAF_table %>%
  filter(nchar(as.character(mutPosition)) <= 50)

#MAF_table %>% filter(HGVSc == "c.993+390C>T") %>% filter(Subject == "Patient")
#df_filtered
# Aggregate the data to find the maximum t_alt_count for each HGVSc and Tissue combination
df_agg <- df_filtered %>%
  filter(Subject == "Patient") %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(InBed_StartPosition == TRUE | InBed_EndPosition == TRUE) %>%
  group_by(mutPosition, Tissue, VAF) %>%
  summarise(max_t_alt_count = max(t_alt_count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(plot_max = if_else(max_t_alt_count > 100, VAF, max_t_alt_count),
         plot_max = round(plot_max, 2)
  )

#df_agg %>% filter(is.na(HGVSc)) %>% print(width = Inf)

df_complete <- df_agg %>%
  #complete(HGVSc, Tissue, fill = list(max_t_alt_count = 0, VAF=0))
  complete(mutPosition, Tissue, fill = list(max_t_alt_count = 0, VAF=0))  # Fill missing combinations with 0


df_agg <- df_complete %>% filter(VAF <= 1) %>% arrange(desc(VAF))
df_agg %>% arrange(desc(VAF)) %>% print(n=100)


mutation_tissue_count <- df_agg %>%
  group_by(mutPosition) %>%
  #group_by(HGVSc) %>%
  summarise(tissue_count = sum(max_t_alt_count > 0), .groups = "drop")
mutation_tissue_count

# Merge the tissue count with df_agg
df_agg <- df_agg %>%
  left_join(mutation_tissue_count, by = "mutPosition")
  #left_join(mutation_tissue_count, by = "HGVSc")

# Order HGVSc based on descending number of tissue matches
df_agg <- df_agg %>%
  mutate(mutPosition = fct_reorder(mutPosition, tissue_count, .desc = FALSE))
  #mutate(HGVSc = fct_reorder(HGVSc, tissue_count, .desc = FALSE))


df_unique <- df_agg %>%
  filter(max_t_alt_count > 0) %>%
  group_by(mutPosition) %>%
  #group_by(HGVSc) %>%
  summarise(
    tissue_count = n_distinct(Tissue),  # Count number of unique tissues
    tissues = list(unique(Tissue)),  # Store unique tissues as a list
    .groups = "drop"
  ) %>%
  arrange(desc(tissue_count))  # Sort by tissue count descending
df_unique

# Convert HGVSc and Tissue to factors to ensure proper ordering in the heatmap
df_agg$mutPosition <- as.factor(df_agg$mutPosition)
#df_agg$HGVSc <- as.factor(df_agg$HGVSc)
df_agg$Tissue <- as.factor(df_agg$Tissue)

# Create a custom color scale
# Black for 0, white-to-red gradient for values >= 1
df_agg
custom_fill_scale <- scale_fill_gradientn(
  colours = c("black", "white", "red"),
  values = scales::rescale(c(0, 1, 30)),  # Rescale to ensure 0 is black and 1 starts white
  limits = c(0, 30),  # Set limits to include 0
  oob = scales::squish
)

df_agg <- df_agg %>% 
  mutate(genomic_pos = word(mutPosition, 2, sep = ":")) %>%
  arrange(desc(genomic_pos)) %>%
  mutate(mutPosition = factor(mutPosition, levels = unique(mutPosition)))

df_agg <- df_agg %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))
custom_label <- function(x) {
  sapply(x, function(t) {
    if (t %in% c("Mediastinal metastasis",
                "Lung metastasis",
                "Esophageal cancer 1",
                "Esophageal cancer 2",
                "Liver metastasis 1",
                "Liver metastasis 2")){
      paste0("<span style='color:red;'>", t, "</span>")
    } else {
      t
    }
  })
}

df_agg %>% filter(VAF > 0.2)
# Create the heatmap using ggplot2
heatmap_position <- ggplot(df_agg, aes(x = Tissue, y = mutPosition, fill = max_t_alt_count)) +
#heatmap <- ggplot(df_agg, aes(x = Tissue, y = HGVSc, fill = max_t_alt_count)) +
  geom_tile() +
  custom_fill_scale +
  geom_text(aes(label = plot_max), color = "black", size = 3) +  # Display values inside tiles
  labs(title = "Mutation Heatmap",
       x = "Tissue",
       y = "Mutation/genomic position (+ strand)",
       fill = "t_alt_count") +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 8))
show(heatmap_position)

ggsave("results/mutation_heatmap_position.png", heatmap_position, width = 12, height = 20, units = "in")
#ggsave("results/mutation_heatmap.png", heatmap, width = 75, height = 25, units = "in", limitsize = FALSE)


################################################################################
############# overlap using just the codon 
################################################################################


#MAF_table %>% filter(HGVSc == "c.993+390C>T") %>% filter(Subject == "Patient")
#df_filtered
# Aggregate the data to find the maximum t_alt_count for each HGVSc and Tissue combination
df_agg <- MAF_table %>%
  filter(Subject == "Patient") %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(InBed_StartPosition == TRUE | InBed_EndPosition == TRUE) %>%
  group_by(prot.pos, Tissue, VAF) %>%
  summarise(max_t_alt_count = max(t_alt_count, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(plot_max = if_else(max_t_alt_count > 100, VAF, max_t_alt_count),
         plot_max = round(plot_max, 2)
  )

#df_agg %>% filter(is.na(HGVSc)) %>% print(width = Inf)

df_complete <- df_agg %>%
  #complete(HGVSc, Tissue, fill = list(max_t_alt_count = 0, VAF=0))
  complete(prot.pos, Tissue, fill = list(max_t_alt_count = 0, VAF=0))  # Fill missing combinations with 0


df_agg <- df_complete %>% filter(VAF <= 1) %>% arrange(desc(VAF))
df_agg %>% arrange(desc(VAF)) %>% print(n=100)


mutation_tissue_count <- df_agg %>%
  group_by(prot.pos) %>%
  #group_by(HGVSc) %>%
  summarise(tissue_count = sum(max_t_alt_count > 0), .groups = "drop")
mutation_tissue_count

# Merge the tissue count with df_agg
df_agg <- df_agg %>%
  left_join(mutation_tissue_count, by = "prot.pos")
#left_join(mutation_tissue_count, by = "HGVSc")

df_agg
# Order HGVSc based on descending number of tissue matches
df_agg <- df_agg %>%
  mutate(prot.pos = fct_reorder(factor(prot.pos), tissue_count, .desc = FALSE))
#mutate(HGVSc = fct_reorder(HGVSc, tissue_count, .desc = FALSE))


df_unique <- df_agg %>%
  filter(max_t_alt_count > 0) %>%
  group_by(prot.pos) %>%
  #group_by(HGVSc) %>%
  summarise(
    tissue_count = n_distinct(Tissue),  # Count number of unique tissues
    tissues = list(unique(Tissue)),  # Store unique tissues as a list
    .groups = "drop"
  ) %>%
  arrange(desc(tissue_count))  # Sort by tissue count descending
df_unique

# Convert HGVSc and Tissue to factors to ensure proper ordering in the heatmap
df_agg$prot.pos <- as.factor(df_agg$prot.pos)
#df_agg$HGVSc <- as.factor(df_agg$HGVSc)
df_agg$Tissue <- as.factor(df_agg$Tissue)

# Create a custom color scale
# Black for 0, white-to-red gradient for values >= 1
df_agg
custom_fill_scale <- scale_fill_gradientn(
  colours = c("black", "white", "red"),
  values = scales::rescale(c(0, 1, 30)),  # Rescale to ensure 0 is black and 1 starts white
  limits = c(0, 30),  # Set limits to include 0
  oob = scales::squish
)

df_agg <- df_agg %>% 
  #mutate(genomic_pos = word(prot.pos, 2, sep = ":")) %>%
  arrange(desc(prot.pos)) %>%
  mutate(prot.pos = factor(prot.pos, levels = unique(prot.pos)))
df_agg
df_agg <- df_agg %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order))
custom_label <- function(x) {
  sapply(x, function(t) {
    if (t %in% c("Mediastinal metastasis",
                 "Lung metastasis",
                 "Esophageal cancer 1",
                 "Esophageal cancer 2",
                 "Liver metastasis 1",
                 "Liver metastasis 2")){
      paste0("<span style='color:red;'>", t, "</span>")
    } else {
      t
    }
  })
}
df_agg
df_agg %>% filter(VAF > 0.2)
# Create the heatmap using ggplot2
heatmap_codon <- ggplot(df_agg, aes(x = Tissue, y = prot.pos, fill = max_t_alt_count)) +
  #heatmap <- ggplot(df_agg, aes(x = Tissue, y = HGVSc, fill = max_t_alt_count)) +
  geom_tile() +
  custom_fill_scale +
  geom_text(aes(label = plot_max), color = "black", size = 3) +  # Display values inside tiles
  labs(title = "Mutation Heatmap",
       x = "Tissue",
       y = "Mutation/genomic position (+ strand)",
       fill = "t_alt_count") +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 8))
show(heatmap_codon)
ggsave("results/mutation_heatmap_codon_test.png", heatmap_codon, width = 12, height = 20, units = "in")

################################################################################
############# 
################################################################################
df_agg <- df_agg %>% 
  filter(tissue_count > 1) %>%
  arrange(tissue_count) %>%
  mutate(Tissue = factor(Tissue, levels = tissue_order),
         mutPosition = factor(mutPosition, levels = unique(mutPosition)),
         plot_max = round(plot_max, 1))

heatmap_overlap <- ggplot(df_agg, aes(x = Tissue, y = mutPosition, fill = max_t_alt_count)) +
  #heatmap <- ggplot(df_agg, aes(x = Tissue, y = HGVSc, fill = max_t_alt_count)) +
  geom_tile() +
  custom_fill_scale +
  geom_text(aes(label = plot_max), color = "black", size = 3) +  # Display values inside tiles
  labs(title = "Mutation Heatmap",
       x = "Tissue",
       y = "Mutation/genomic position (+ strand)",
       fill = "t_alt_count") +
  theme_minimal() +
  scale_x_discrete(labels = custom_label) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme(
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 8))
show(heatmap_overlap)

ggsave("results/mutation_heatmap_overlap.png", heatmap_overlap, width = 8, height = 5, units = "in")


#####################
manual_order <- tibble(
  manual_order = c(1,2,3,4,5),
  Start_Position = c(7674220, 7673587, 7674894, 7676387, 7673535)
)

MAF_table %>% filter(Hugo_Symbol == "TP53") %>% print(width=Inf)
top_common_muts <- MAF_table %>% 
  filter(mutPosition == "chr17:7674220:C:T" | mutPosition == "chr17:7673587:G:A" | mutPosition == "chr17:7674894:G:A" | mutPosition == "chr17:7676387:C:T" | mutPosition == "chr17:7673535:C:T") %>% 
  dplyr::select(Start_Position, Variant_Classification, Variant_Type, HGVSp_Short, am_pathogenicity, am_class, Exon_Number) %>% 
  left_join(manual_order, by = "Start_Position") %>%
  print()
top_common_muts %>% group_by(across(everything())) %>% summarise(n_tissues = n()) %>% arrange(manual_order) %>% dplyr::select(-manual_order)
