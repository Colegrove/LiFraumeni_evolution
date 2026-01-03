## Li-Fraumeni chip panel skyscraper plot

family <- c("Family member A", "Family member B", "Family member C")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

heatmap_prep <- 
  maf_masked_coding %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>% # only chip panel
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>% # exclude non-coding
  filter(!inRepeatMask| Variant_Classification == "Splice_Site") %>% #exclude repeat masking
  group_by(Subject, Hugo_Symbol) %>%
  summarise(total_mut_count = n(), 
            total_burden_count = sum(t_alt_count), .groups = "drop") %>%
  mutate(Subject = factor(
    Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
  ))


heatmap_prep <-
  maf_masked_coding %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>% #exclude repeat masking
  group_by(Subject, Tissue) %>%
  mutate(SampCodingOrder = rank(dplyr::desc(t_alt_count), ties.method = "first"))

heatmap_prep



################################################################################
############################# heatmaps #########################################
################################################################################



##############################################
##### heatmap by gene and total reads
##############################################
gene_order <- heatmap_prep %>% 
  group_by(Hugo_Symbol) %>% 
  summarise(count = n()) %>% 
  arrange(count) %>%
  pull(Hugo_Symbol)

heatmap_data <- heatmap_prep %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Subject, Hugo_Symbol) %>%
  summarise(Total_t_alt_count = sum(t_alt_count, na.rm = TRUE), .groups = "drop")   %>%
  group_by(Hugo_Symbol) %>%
  mutate(sum_count = sum(Total_t_alt_count)) %>%
  ungroup() %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = gene_order))
  
heatmap_data

# skyscraper_blood_ms_3A.R
sample_codes <- ann_wide %>% dplyr::select(Subject, SampleCode)
heatmap_data <- heatmap_data %>% left_join(sample_codes, by="Subject") %>% 
  mutate(SampleCode = factor(SampleCode, levels = c("CON01", "CON04", "CON03", "CON02", "LFS01", "CON05", "LFS02", "CON06", "REL01", "LFS03", "CON07")))
  
heatmap_data
heatmap_plot <- ggplot(heatmap_data, aes(x = SampleCode, y = Hugo_Symbol, fill = Total_t_alt_count)) +
  geom_tile(color = "#DDDDDD") +  # Add white borders between tiles
  geom_text(aes(label = Total_t_alt_count), color = "black", size = 8*25.4/72.27) +
  scale_fill_gradient(low = "#DDDDDD", high = "#88CCEE", na.value = "#DDDDDD", limits=c(10,150), oob = scales::squish) +  # Heatmap color scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8, margin=margin(-4,0,0,0)),
        axis.text.y = element_text(size = 8, margin=margin(0,0,0,0)),
        panel.spacing = unit(0, "pt"),
        strip.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 8, angle=45, vjust = 0.5),
        legend.title = element_text(size = 8, margin = margin(b=5, r=2), vjust = 1),
        #legend.margin = margin(l=-9, b=-9),
        legend.margin = margin(-8,0,0,0),
        legend.key.size = unit(8, "pt"),
        legend.key.spacing = unit(0.7, "pt"),
        plot.margin=margin(t=1, l=1, r=1)) +
  labs(fill = "Count")

# Show the heatmap
show(heatmap_plot)

ggsave("results/CHIP_gene_heatmap_total_count.png", heatmap_plot, width = 3.1, height = 3.3, units = "in", dpi = 300)


##############################################
##### heatmap by gene and unique mutation
##############################################

# heatmap_data <- skyscraper_prep %>%
#   group_by(Subject, Hugo_Symbol) %>%
#   filter(Hugo_Symbol != "Unknown") %>%
#   summarise(Unique_Mutations = n_distinct(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2), 
#             .groups = "drop") %>%
#   group_by(Hugo_Symbol) %>%
#   mutate(total_mut = sum(Unique_Mutations)) %>%
#   ungroup() %>%
#   mutate(Hugo_Symbol = fct_reorder(Hugo_Symbol, total_mut)) %>%
#   mutate(Subject = factor(
#     Subject, levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
#   ))
# 
# 
# heatmap_plot <- ggplot(heatmap_data, aes(x = Subject, y = Hugo_Symbol, fill = Unique_Mutations)) +
#   
#   geom_tile(color = "white") +  # Add white borders between tiles
#   geom_text(aes(label = Unique_Mutations), color = "black", size = 3, fontface = "bold") +
#   
#   scale_fill_gradient(low = "white", high = "red", na.value = "white", oob = scales::squish) +  # Heatmap color scale
#   
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
#         axis.text.y = element_text(size = 8),
#         panel.spacing = unit(0, "pt"),
#         strip.text.x = element_text(size = 10, face = "bold"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()) +
#   
#   labs(fill = "Unique Mutations")
# 
# # Show the heatmap
# show(heatmap_plot)
# 
# subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Patient", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
# ages <- c("20-30", "20-30", "20-30", "20-30", "30s", "39", "61", "60s", "69", "70s", "70s")
# LFS <- c("-", "-", "-", "-", "+", "+", "-","-", "+","-", "-")
# CTx <- c("-", "-", "-", "-", "+", "-", "-","-", "-","-", "+")
# 
# age_lf <- data.frame(
#   x_axis_var = subjects, 
#   CTx = CTx,
#   LFS = LFS,
#   age = ages
# )
# 
# age_lf_long <- age_lf %>%
#   pivot_longer(cols = c(age, LFS, CTx), names_to = "variable", values_to = "value") %>%
#   mutate(y = case_when(
#     variable == "age" ~ 2,
#     variable == "LFS" ~ 1.5,
#     variable == "CTx" ~ 1, 
#   )) %>% 
#   mutate(x_axis_var = factor(x_axis_var, levels = unique(x_axis_var)),
#          x_axis_var = factor(x_axis_var, levels = unique(x_axis_var))  # preserve original order
#   )
# 
# 
# annotation_plot <- ggplot(age_lf_long, aes(x = x_axis_var, y = y)) +
#   geom_text(aes(label = value), size = 4, vjust = 0) +
#   scale_y_continuous(
#     breaks = c(1.1, 1.6, 2.1),
#     labels = c("CTx", "LFS", "age"),
#     #expand = c(0, 0),
#     limits = c(1, 2.3)
#   ) +
#   scale_x_discrete(expand = c(0.07, 0)) +
#   theme_minimal(base_size = 8) +
#   theme(
#     axis.title = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_text(size = 8),
#     panel.grid = element_blank(),
#     plot.margin = margin(0, 0, 0, 0, "pt")
#   )
# show(annotation_plot)
# 
# 
# combined_plot <- cowplot::plot_grid(
#   heatmap_plot + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = "top"),
#   annotation_plot,
#   ncol = 1,
#   rel_heights = c(1, 0.15),
#   greedy = TRUE
# )
# show(combined_plot)
# ggsave("results/CHIP_gene_heatmap_unique_muts.png", combined_plot, width = 6, height = 6, units = "in", dpi = 300)
# 
# 





