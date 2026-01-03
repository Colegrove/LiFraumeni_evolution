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



## mutFreq by gene
mutFreq_counts <- 
  maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  mutate(mutType = paste0(Reference_Allele,">",Tumor_Seq_Allele2)) %>%
  group_by(Tumor_Sample_Barcode, Subject, Tissue, Hugo_Symbol, age, mutType) %>% 

  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  ) %>%
  mutate(
    n_muts   = replace_na(n_muts, 0),
    mutReads = replace_na(mutReads, 0),,
    LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  filter(str_detect(mutType, "^[ACGT]>[ACGT]$")) %>%
  group_by(LFS, mutType) %>%
  summarise(count = n()) %>% print(n=Inf)

mutType_bar <- ggplot(mutFreq_counts, 
                      aes(x = mutType, y = count, fill = LFS, color=LFS)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Mutation Type",
       y = "Count",
       fill = "Patient history",
       color = "Patient history") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    axis.text.y = element_text(size=8),
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = LFS_colors) +
  scale_color_manual(values = LFS_colors)
mutType_bar

df_prop <- mutFreq_counts %>%
  group_by(LFS) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

mutType_prop_plot <- ggplot(df_prop, 
                            aes(x = mutType, y = prop, fill = LFS, color = LFS)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits= c(0,0.2)) +
  labs(x = "Mutation Type",
       y = "Proportion",
       fill = "Patient history",
       color = "Patient history") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    axis.text.y = element_text(size=8),
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = LFS_colors) +
  scale_color_manual(values = LFS_colors)



    

mutFreq_counts_non_coding_mutagenesis <- 
  maf_masked_coding %>%
  filter(!is.na(MUT_region_StartPosition) | !is.na(MUT_region_EndPosition)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  mutate(mutType = paste0(Reference_Allele,">",Tumor_Seq_Allele2)) %>%
  print(width = Inf) %>%
  group_by(Tumor_Sample_Barcode, Subject, Tissue, Hugo_Symbol, age, mutType) %>% 
    
    summarise(
      n_muts   = n(),
      mutReads = sum(t_alt_count),
      .groups  = "drop"
    ) %>%
    mutate(
      n_muts   = replace_na(n_muts, 0),
      mutReads = replace_na(mutReads, 0),,
      LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS"),
      CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
    filter(str_detect(mutType, "^[ACGT]>[ACGT]$")) %>%
    group_by(LFS, mutType) %>%
    summarise(count = n()) %>% print(n=Inf)
      
  
  






mutType_bar <- ggplot(mutFreq_counts_non_coding_mutagenesis, 
                      aes(x = mutType, y = count, fill = LFS, color=LFS)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous() +
  labs(x = "Mutation Type",
       y = "Count",
       fill = "Patient history",
       color = "Patient history") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    axis.text.y = element_text(size=8),
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = LFS_colors) +
  scale_color_manual(values = LFS_colors)
mutType_bar

df_prop <- mutFreq_counts_non_coding_mutagenesis %>%
  group_by(LFS) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

mutType_prop_plot <- ggplot(df_prop, 
                            aes(x = mutType, y = prop, fill = LFS, color = LFS)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits= c(0,0.2)) +
  labs(x = "Mutation Type",
       y = "Proportion",
       fill = "Patient history",
       color = "Patient history") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=8),
    axis.text.y = element_text(size=8),
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = LFS_colors) +
  scale_color_manual(values = LFS_colors)











#coding = replace_na(coding, "coding"), 
    #age = age_map[subject]) %>%
    #dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")


## write mutation frequencies by gene/subject
write.csv(mutFreq_combined, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene.csv")

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
  mutate(coding = "non-coding-total")

mutFreq_prep <- rbind(mutFreq_prep, total_non_coding)



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

lm_model <- lm(mutFreq ~ age, data = mutFreq_prep %>% filter(coding == "coding"))
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_coding <- ggplot(mutFreq_prep %>% filter(coding == "coding"), aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'coding'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(1.3e-7, 4e-7)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Coding\nmutation frequency"
  ) +
  theme_minimal() 

mutFreq_coding
###############################################################################
### Non-coding plot
###############################################################################

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

mutFreq_non_coding <- ggplot(mutFreq_prep %>% filter(coding == "non-coding-total"), aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep %>% filter(coding == 'non-coding-total'),
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(limits = c(1.5e-6, 3e-5), labels=fancy_scientific) +
  scale_y_log10(limits = c(1.3e-7, 4.2e-7)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Non-coding\nmutation frequency"
  ) +
  theme_minimal()

mutFreq_non_coding
################################################################################
########### Plot frequency and Burden combined
################################################################################

mutFreq_coding <- mutFreq_coding + theme(legend.position = "none")
mutFreq_non_coding <- mutFreq_non_coding + theme(legend.position = "none")

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

mutFreq_non_coding <- mutFreq_non_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  ylab(expression(atop(NA, atop(textstyle("Non-coding"), textstyle("mutation frequency"))))) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_coding <- mutFreq_coding + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle("Coding"), textstyle("mutation frequency"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

plots <- plot_grid(
  plot_grid(mutFreq_coding, mutFreq_non_coding, ncol = 2, align = "v"), 
  ncol = 1)

combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

ggsave("results/MF_CHIP_ms2BC.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
ggsave("/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_2/MF_CHIP_ms2BC.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)

