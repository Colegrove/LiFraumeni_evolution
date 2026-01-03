# ## Li-Fraumeni mutation frequencies all chip genes
#install.packages(c("ggplot2", "scales", "rlang"))
#### First run 
library(tidyverse)
library(ggtext)

### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")
family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

CHIP_genes <- c("NRAS", "BRINP3", "DNMT3A", "IDH1", "GATA2", "KIT", "TET2", "NPM1", 
                "EZH2", "RAD21", "HNRNPK", "PTEN", "SMC3", "WT1", "KMT2A", "CBL", "KRAS", 
                "PTPN11", "FLT3", "IDH2", "MYH11", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "SMC1A", 
                "STAG2", "PHF6", "TP53")

sample_id_mapping_path <- "/Volumes/feder-vol1/project/li_fraumeni/dat/sampleID_mapping.txt"
sample_map <- read_delim(sample_id_mapping_path, delim = "\t", quote="\"") %>%
  mutate(tissue = str_trim(str_replace_all(tissue, '"', '')), 
         subject = str_trim(str_remove(subject, ":$")))



###################################################
########### Pull sequencing depths in coding regions
###################################################

GEOM_POINT_SIZE = 1.5

age_map <- c("UW volunteer 1" = 25,
             "UW volunteer 2" = 30,
             "UW volunteer 3" = 27,
             "UW volunteer 4" = 25,
             "Patient" = 35,
             "Family member A" = 39,
             "Family member B" = 61,
             "UW volunteer 5" = 37,
             "Family member C" = 69,
             "UW volunteer 6" = 60,
             "UW volunteer 7" = 76)

## by gene
chip_depths <- final_masked_depth %>% 
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  filter(tissue == "Urine cells") %>% 
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator_coding = sum(DP), .groups = "drop")

chip_depths_non_coding <- final_masked_depth %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator_coding = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)


## mutFreq by gene
mutFreq_counts <- 
  maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% c("Urine cells")) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, age) %>% 
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )



  
mutFreq_combined <- chip_depths %>%
  left_join(mutFreq_counts, 
            by = c("Samp" = "Tumor_Sample_Barcode",
                   "gene_name" = "Hugo_Symbol")) %>%
  mutate(
    n_muts   = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),
    mutFreq   = n_muts   / denominator_coding,
    mutBurden = mutReads / denominator_coding,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
    dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

## write mutation frequencies by gene/subject
write.csv(mutFreq_combined, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/urine_muts_by_gene.csv")

mutFreq_counts_non_coding <- 
  maf_masked_coding %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number) & Variant_Classification != "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% c("Urine cells")) %>%
  filter(!Variant_Classification == "Splice_Site") %>%
  mutate(age = age_map[Subject]) %>%
  #filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, age) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )

mutFreq_combined_non_coding <- chip_depths_non_coding %>%
  left_join(mutFreq_counts_non_coding, 
            by = c("Samp" = "Tumor_Sample_Barcode",
                   "gene_name" = "Hugo_Symbol")) %>%
  mutate(
    n_muts   = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),
    mutFreq   = n_muts   / denominator_coding,
    mutBurden = mutReads / denominator_coding,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

## write mutation frequencies by gene/subject
write.csv(mutFreq_combined_non_coding, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/urine_muts_by_gene_non_coding.csv")


###############################################################################
### Mutation frequency
###############################################################################


mutFreq_combined
## combine genes for plot
mutFreq_prep <- mutFreq_combined %>% 
  group_by(Subject, tissue, age, LFS, CTx) %>%
  summarise(denominator_coding = sum(denominator_coding), 
            n_muts = sum(n_muts),
            mutReads = sum(mutReads)) %>%
  mutate(mutFreq = n_muts / denominator_coding,
         mutBurden = mutReads / denominator_coding)
mutFreq_prep
## offsets are used for plotting manual jitter of individuals with same age
offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset)

mutFreq_prep <- offsets %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS" & CTx == "CTx" ~ "LFS/CTx",   # both LFS and CTx
      LFS == "LFS"             ~ "LFS/no-CTx",  # just LFS
      LFS == "non-LFS" & CTx == "CTx" ~ "non-LFS/CTx",
      TRUE                 ~ "non-LFS/no-CTx"      # if you want a default
    )
  )

patient_history_order = c("LFS/no-CTx", "LFS/CTx", "non-LFS/no-CTx", "non-LFS/CTx")
mutFreq_prep$shape_group = factor(mutFreq_prep$shape_group, levels = patient_history_order)

###############################################################################
### Mutation frequency plot
###############################################################################


mutFreq_prep <- mutFreq_prep %>% 
  mutate(Subject_abbr = if_else(Subject == "Family member A", "LFS02", "REL01"))

mutFreq_plot2 <- ggplot(mutFreq_prep, aes(x = Subject_abbr, y = mutFreq)) +
  geom_bar(stat = "identity") +
  labs(
    y = "TP53\nmutation frequency"
  ) +
  theme_minimal() + 
  theme(axis.title.x = element_blank())

show(mutFreq_plot2)

###############################################################################
### Mutation burden plot
###############################################################################


mutFreq_prep
mutBurden_plot2 <- ggplot(mutFreq_prep, aes(x = Subject_abbr, y = mutBurden)) +
  geom_bar(stat = "identity") +
  labs(
    y = "Overall\nmutation burden"
  ) +
  theme_minimal()  + 
  theme(axis.title.x = element_blank())


show(mutBurden_plot2)



################# non-Coding
mutFreq_counts_non_coding

################################################################################
########### Plot frequency and Burden combined
################################################################################

# mutBurden_plot2 <- mutBurden_plot2 + theme(legend.position = "none")
# mutFreq_plot2 <- mutFreq_plot2 + theme(legend.position = "none")
# 
# legend_shared <- get_legend(
#   mutFreq_plot2 + 
#     guides(
#       color = guide_legend(nrow = 2, title.position = "top"),  # force 1 row
#       shape = guide_legend(nrow = 2, title.position = "top")   # force 1 row
#     ) +
#     theme(legend.position = "top",
#           legend.title = element_blank(),
#           legend.text  = element_text(size = 8),
#           legend.margin = margin(r=1, l=1, b = -20, t = -25),
#           legend.key.size   = unit(0.5, "lines"),
#           legend.spacing.x  = unit(0.2, "cm"),
#           legend.spacing.y  = unit(0.2, "cm")
#     )
# )
# 
# mutBurden_plot2 <- mutBurden_plot2 + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   ylab(expression(atop(NA, atop(textstyle("Overall"), textstyle("mutation burden"))))) +
#   theme(legend.position = "none",
#         text = element_text(size = 8),
#         axis.title = element_text(size = 8),
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
#         axis.text  = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.text  = element_text(size = 8))
# 
# mutFreq_plot2 <- mutFreq_plot2 + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
#   ylab(expression(atop(NA, atop(textstyle("Overall"), textstyle("mutation frequency"))))) +
#   theme(legend.position  = "none",
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
#         axis.text  = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.text  = element_text(size = 8))
# 
# plots <- plot_grid(
#   plot_grid(mutFreq_plot2, mutBurden_plot2, ncol = 2, align = "v"), 
#   ncol = 1)
# 
# combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
#                             rel_heights = c(1.1,0.225))
# 
# show(combined_plots)
# ggsave("results/CHIP_full_mutation_frequency.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)





mutFreq_plot3 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq)) +
  geom_point(aes(color = LFS, fill = LFS, shape = plot_coding),
             size = GEOM_POINT_SIZE + 1, 
             stroke = 0.8) +
  geom_smooth(
    aes(x=age, y=mutFreq, color = LFS, linetype = plot_coding, group = interaction(LFS, plot_coding)),
    method = "lm", se = FALSE, linewidth = 0.7
  ) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(labels = fancy_scientific, limits = c(1e-7, 3.1e-6)) +
  scale_y_log10(limits = c(3e-7, 1e-5)) +
  scale_color_manual(values = LFS_colors, name = "Patient history") +
  scale_fill_manual(values = LFS_colors, name = "Patient history") +
  scale_shape_manual(values = coding_shapes, name = "Region") +
  labs(x = "Age", y = "Mutation frequency") +
  theme_minimal() + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text  = element_text(size = 8, margin=margin(r=3, l=2)),
        legend.margin = margin(r=1, l=1, b = 1, t = 1),
        legend.box.margin = margin(r=0, l=-30, b = 0, t = 0),
        legend.key.size   = unit(8, "pt"),
        legend.spacing.x  = unit(0.2, "cm"),
        legend.spacing.y  = unit(0.2, "cm"),
        #legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3), 
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8)
  ) + 
  guides(linetype = guide_legend(override.aes = list(color = 'black'), ncol=1),
         shape = guide_legend(ncol = 1),
         color = guide_legend(ncol =1))
mutFreq_plot3
