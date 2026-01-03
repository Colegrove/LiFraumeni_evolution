# ## Li-Fraumeni mutation frequencies all chip genes
#install.packages(c("ggplot2", "scales", "rlang"))
#### Run first

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
             "Patient" = 34,
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
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

final_masked_depth %>% print(width = Inf)

chip_depths_non_coding <- final_masked_depth %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

## mutFreq by gene
mutFreq_counts <- 
  maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
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
    mutFreq   = n_muts   / denominator,
    mutBurden = mutReads / denominator,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
    dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

## write mutation frequencies by gene/subject
write.csv(mutFreq_combined, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene.csv")
mutFreq_combined %>% filter(Hugo_Symbol == "TP53")
#### write non-coding
mutFreq_counts_non_coding <- 
  maf_masked_coding %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(!Variant_Classification == "Splice_Site") %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
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
    n_muts = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),
    mutFreq  = n_muts   / denominator,
    mutBurden = mutReads / denominator,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")

## write mutation frequencies by gene/subject
write.csv(mutFreq_combined_non_coding, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene_non_coding.csv")
mutFreq_combined_non_coding %>% filter(Hugo_Symbol == "TP53")
mutFreq_combined_coding <- mutFreq_combined %>% mutate(coding = "coding")
mutFreq_combined_non_coding <- mutFreq_combined_non_coding %>% mutate(coding = "non-coding-CHIP")
MF_CHIP_genic <- rbind(mutFreq_combined_coding, mutFreq_combined_non_coding)

write.csv(MF_CHIP_genic, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene_non_coding.csv")

########## Add mutagenesis region depths and mutation counts to non_coding genes
## first aggregate genes 
mutFreq_prep <- MF_CHIP_genic %>% 
  group_by(Subject, tissue, age, LFS, CTx, coding) %>%
  summarise(denominator = sum(denominator), 
            n_muts = sum(n_muts),
            mutReads = sum(mutReads)) %>%
  mutate(mutFreq = n_muts / denominator,
         mutBurden = mutReads / denominator)

## add mutagenesis regions depth
mutagenesis_depths <- final_masked_depth %>%
  filter(str_starts(Gene, "region")) %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(denominator = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples) %>% 
  mutate(coding = "non-coding-MUT")

## find mutagenesis mutation counts
mutFreq_counts_non_coding_mutagenesis <- 
  maf_masked_coding %>%
  filter(!is.na(MUT_region_StartPosition) | !is.na(MUT_region_EndPosition)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  group_by(Subject, Tissue) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )


MF_MUT <- mutFreq_counts_non_coding_mutagenesis %>%
  left_join(mutagenesis_depths, by = c("Subject" = "subject", "Tissue" = "tissue"))

apply_data <- mutFreq_prep %>% filter(coding == 'coding') %>% dplyr::select(Subject, tissue, age, LFS, CTx)
MF_MUT <- MF_MUT %>% left_join(apply_data, by = c("Tissue" = "tissue", 'Subject')) %>%
  mutate(mutFreq = n_muts/denominator,
         mutBurden = mutReads/denominator) %>%
  dplyr::select(-Samp) %>%
  dplyr::rename(tissue = Tissue)

mutFreq_prep <- rbind(mutFreq_prep, MF_MUT)
mutFreq_prep
## calculate total non-coding (CHIP + MUT)
total_non_coding <- mutFreq_prep %>% filter(coding != "coding") %>%
  group_by(Subject, tissue, age, LFS, CTx) %>%
  mutate(denominator_total = sum(denominator), 
         n_muts_total = sum(n_muts), 
         mutReads_total = sum(mutReads),
         mutFreq_total = n_muts_total/denominator_total,
         mutBurden_total = mutReads_total/denominator_total) %>%
  dplyr::select(-denominator, -n_muts, -mutReads, -mutFreq, -mutBurden) %>%
  dplyr::rename(denominator = denominator_total, n_muts = n_muts_total, 
                mutReads = mutReads_total, mutFreq = mutFreq_total, 
                mutBurden = mutBurden_total) %>%
  mutate(coding = "non-coding-total") %>%
  distinct()

total_non_coding
mutFreq_prep <- rbind(mutFreq_prep, total_non_coding)

mutFreq_prep 

###############################################################################
### Mutation frequency
###############################################################################




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

## function to turn "e-6" to "10^-6"
fancy_scientific <- function(l) { 
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text = l)
}

LFS_color <- "#882255"
nonLFS_color <- "#44aa99"

# Define colors by LFS status
shape_group_colors <- c(
  "non-LFS/no-CTx" = "#44AA99",  # non-LFS color
  "non-LFS/CTx"     = "#44AA99",  # same as non-LFS
  "LFS/no-CTx"     = "#882255",  # LFS color
  "LFS/CTx" = "#882255"   # same as LFS
)

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS/no-CTx" = 1,
  "non-LFS/CTx"     = 16,
  "LFS/no-CTx"     = 1,
  "LFS/CTx" = 16
)

## coding model

lm_model <- lm(mutBurden ~ age, data = mutFreq_prep %>% filter(coding == "coding"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_prep %>% filter(coding == "coding")
mutFreq_coding <- ggplot(mutFreq_prep %>% filter(coding == "coding"), aes(x = age_j, y = mutBurden)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'coding'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(1e-7, 1e-5)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Coding\nmutation burden"
  ) +
  # annotate(
  #   "text",
  #   x = 60, y = 1.8e-07,
  #   label = paste0("R^2 == ", round(r2, 3)),,
  #   parse = TRUE,
  #   size = 2.8
  # ) +
  annotate(
    "text",
    x = 60, y = 1.4e-07,
    label = paste0("italic(p) == ", signif(pval, 2)),,
    parse = TRUE,
    size = 2.8
  ) +
  theme_minimal() 

mutFreq_coding


################################################################################
########### 
################################################################################

mutFreq_coding <- mutFreq_coding + theme(legend.position = "none")
#mutFreq_non_coding <- mutFreq_non_coding + theme(legend.position = "none")

legend_shared <- get_legend(
  mutFreq_coding + 
    guides(
      color = guide_legend(nrow = 2, title.position = "top"),  # force 1 row
      shape = guide_legend(nrow = 2, title.position = "top")   # force 1 row
    ) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text  = element_text(size = 8),
          legend.margin = margin(r=1, l=1, b = -20, t = -25),
          legend.key.size   = unit(0.5, "lines"),
          legend.spacing.x  = unit(0.2, "cm"),
          legend.spacing.y  = unit(0.2, "cm")
    )
)

# mutFreq_non_coding <- mutFreq_non_coding + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   ylab(expression(atop(NA, atop(textstyle("Mutagenesis"), textstyle("mutation frequency"))))) +
#   theme(legend.position = "none",
#         text = element_text(size = 8),
#         axis.title = element_text(size = 8),
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
#         axis.text  = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.text  = element_text(size = 8))

mutFreq_coding <- mutFreq_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle("CHIP"), textstyle("mutation burden"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

# plots <- plot_grid(
#   plot_grid(mutFreq_coding, mutFreq_non_coding, ncol = 2, align = "v"), 
#   ncol = 1)

combined_plots <- plot_grid(mutFreq_coding, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

combined_plots
#ggsave(paste0("results/", file_out), combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
#ggsave(paste0("/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_2/", file_out), combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)







################################################################################
########### 
################################################################################


mutFreq_subject <- mutFreq_prep_CHIP_encode %>%
  mutate(age_decades = age / 10,
         depth_scaled = denominator / 10000000) %>%
  mutate(MF_scale = mutReads/depth_scaled)
mutFreq_prep_CHIP_encode

# With UW07 and no CTx
# model_freq_CHIP_all_samples_noCTx <- lm(
#   #n_muts ~ age_decades + depth_scaled + LFS_n,
#   #MF ~ age_decades + depth_scaled + LFS_n, 
#   MF_scale ~ age_decades + depth_scaled + LFS_n,
#   data = mutFreq_subject
# )
# summary(model_freq_CHIP_all_samples_noCTx)

#non_coding_set = 0 # non-coding-total
#non_coding_set = 1 # non-coding-MUT
non_coding_set = 2 # non-coding-CHIP

if(non_coding_set == 0){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-total"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-total")
  xlimits = c(-1.2, 1.5)
  xbreaks = seq(-1, 1.5, 0.5)
  #ylabel = "Non-coding\nCHIP + MUT MF"
  file_out = "MF_multi_coding_non_coding_total_ms.png"
}
if(non_coding_set == 1){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-MUT"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-MUT")
  xlimits = c(-1.6, 1.63)
  xbreaks = seq(-1.5, 1.5, 0.5)
  #ylabel = "MUT MF"
  file_out = "MF_multi_coding_non_coding_MUTonly_ms.png"
}
if(non_coding_set == 2){
  #lm_model <- lm(mutFreq ~ age, data = mutFreq_subject %>% filter(coding == "non-coding-CHIP"))
  mutFreq_subject_non_coding = mutFreq_subject %>% filter(coding == "non-coding-CHIP")
  xlimits = c(-0.5, 1.5)
  xbreaks = seq(-0.5, 1.5, 0.5)
  #ylabel = "non-coding\nCHIP MF"
  file_out = "MF_multi_coding_non_coding_CHIPonly_ms.png"
}

## coding
mutFreq_subject_coding <- mutFreq_subject %>%
  filter(coding == "coding")

mutFreq_subject_coding
model_freq_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_coding
)
summary(model_freq_coding)


model_freq_non_coding <- glm(
  #n_muts ~ age_decades + depth_scaled + LFS_n,
  #MF ~ age_decades + depth_scaled + LFS_n, 
  #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
  MF_scale ~ 1 + age_decades + LFS_n + CTx_n,
  data = mutFreq_subject_non_coding
)
summary(model_freq_non_coding)
mutFreq_subject_non_coding

coef_freq_CHIP_coding <- tidy(model_freq_coding, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(
    p.value < 0.0001 ~ "***",
    p.value < 0.01   ~ "**",
    p.value < 0.05   ~ "*",
    TRUE ~ ""
  ))

coef_freq_CHIP_coding
term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
CHIP_freq_plot_coding <- ggplot(coef_freq_CHIP_coding,
                                aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
  geom_point(color = "black", size = 1) +
  geom_errorbarh(height = 0, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-25, 45), breaks = seq(-25, 45, 5)) +
  geom_text(aes(label = sig, x = estimate),
            hjust = 0.5, vjust = 0.2, size = 5, color = "black") +
  labs(x = expression("Effect size (mutant bases / "  * 10^7 * " bases)"), y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.y = element_markdown(size = 8),
        axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))

CHIP_freq_plot_coding

#ggsave("results/MF_multi_coding_ms.png", CHIP_freq_plot_coding, width = 1.5, height = 1.5, units = "in", dpi = 300)

# ####### freq + burden plot
# 
# CHIP_freq_plot_coding <- CHIP_freq_plot_coding + theme(
#   plot.margin = margin(1,1,1,1)) + 
#   coord_cartesian(clip = "off")
# CHIP_freq_plot_non_coding <- CHIP_freq_plot_non_coding + theme(
#   plot.margin = margin(1,1,1,1)) + 
#   coord_cartesian(clip = "off")
# # CHIP_burden_plot <- CHIP_burden_plot + theme(
# #   plot.margin = margin(1,1,1,1)
# # )
# 
# 
# combined_plot <- CHIP_freq_plot_coding | CHIP_freq_plot_non_coding
# 
# combined_plot



combined_plot <- combined_plots / CHIP_freq_plot_coding

combined_plot <- combined_plot + theme(plot.margin = margin(0,0,10,0))
combined_plot
ggsave(paste0("results/", "MB_regression_supp.png"), combined_plot, width = 3.5, height = 3.5, units = "in", dpi = 300)
#ggsave(paste0("/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_2/", file_out), combined_plot, width = 3.5, height = 1, units = "in", dpi = 300)


# save supp table
# 
# model_to_row <- function(model, gene, coding) {
#   
#   coding_tidy <- broom::tidy(model)
#   
#   coding_wide <- coding_tidy %>%
#     dplyr::select(term, estimate, std.error, p.value) %>%
#     dplyr::mutate(
#       est_se = sprintf("%.4g (%.4g)", estimate, std.error),
#       p      = p.value
#     ) %>% 
#     dplyr::select(term, est_se, p) %>% 
#     dplyr::filter(term %in% c("(Intercept)", "age_decades", "LFS_n", "CTx_n")) %>% 
#     tidyr::pivot_wider(
#       names_from  = term,
#       values_from = c(est_se, p),
#       names_glue  = "{term}.{.value}"
#     ) %>% 
#     dplyr::mutate(
#       model  = "MF ~ 1 + decades + LFS + CTx",
#       gene   = gene,
#       coding = coding,
#       .before = 1
#     ) %>% 
#     dplyr::rename(
#       `(Intercept).Estimate` = `(Intercept).est_se`,
#       `decades.Estimate`     = `age_decades.est_se`,
#       `LFS.Estimate`         = `LFS_n.est_se`,
#       `CTx.Estimate`         = `CTx_n.est_se`,
#       `decades.p`            = `age_decades.p`,
#       `LFS.p`                = `LFS_n.p`,
#       `CTx.p`                = `CTx_n.p`
#     ) %>%
#     dplyr::select(
#       model, gene, coding,
#       `(Intercept).Estimate`, `(Intercept).p`,
#       `decades.Estimate`,     `decades.p`,
#       `LFS.Estimate`,         `LFS.p`,
#       `CTx.Estimate`,         `CTx.p`
#     )
#   
#   coding_wide
# }
# row_coding <- model_to_row(model_freq_coding,      gene = "CHIP", coding = "coding")
# row_noncoding <- model_to_row(model_freq_non_coding,   gene = "CHIP", coding = "non-coding")
# model_table <- bind_rows(row_coding, row_noncoding)
# 
# write.table(
#   model_table,
#   file = "results/MF_multivariate_CHIP_summary.csv",
#   sep = ",",
#   row.names = FALSE
# )







