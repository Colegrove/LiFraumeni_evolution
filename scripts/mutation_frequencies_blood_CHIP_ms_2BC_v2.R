# ## Li-Fraumeni mutation frequencies all chip genes

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
             "UW volunteer 2" = 25,
             "UW volunteer 3" = 25,
             "UW volunteer 4" = 25,
             "Patient" = 35,
             "Family member A" = 39,
             "Family member B" = 61,
             "UW volunteer 5" = 65,
             "Family member C" = 69,
             "UW volunteer 6" = 75,
             "UW volunteer 7" = 75)


## use newly annotated coding/non-coding
CHIP_depth_full_final <- final %>%
  filter(in_CDS) %>%
  group_by(Samp) %>% 
  summarise(denominator_coding = sum(DP))
CHIP_depth_full_final %>% print(n=Inf)

## gene based with new annotations
chip_depths <- final %>% 
  filter(in_CDS) %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue, gene_name) %>%
  summarise(denominator_coding = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)
chip_depths

## mutFreq by gene
mutFreq_counts <-
  filt_maf_CHIP %>%
  filter(coding_from_maf == "coding") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Tumor_Sample_Barcode, Hugo_Symbol, coding, age) %>%
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
    coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
    rename("subject" = 'Subject', "gene_name" = "Hugo_Symbol")
mutFreq_combined
## write mutation frequencies by gene/subject
write.csv(mutFreq_combined, "/Users/huntc10/Desktop/CHIP_muts_by_gene_transcript.csv")

######### burden by gene
mutBurden_prep <-
  filt_maf_CHIP %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(InBed == TRUE) ## only mutations within the bed targets to match depth denominator above

lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")
chip_depths <- chip_depths %>% dplyr::select(Samp, gene_name, denominator_coding)
### filter for plot
mutBurden_prep <- mutBurden_prep %>%
  left_join(chip_depths, by=c("Hugo_Symbol" = "gene_name", "Tumor_Sample_Barcode" = "Samp" )) %>%
  print(width = Inf) %>%
  group_by(coding, Subject, denominator_coding, age, Hugo_Symbol) %>%
  summarise(n_muts = n(),mutReads = sum(t_alt_count)) %>%
  mutate(mutFreq = n_muts/denominator_coding) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  print(n=Inf)

###############################################################################
### Mutation frequency
###############################################################################

#### mutation frequency # of mutations/denominator
family_patient_blood_samples <- c("PBMC", "Buffy coat")
mutFreq_prep <-
  filt_maf_CHIP %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(InBed == TRUE) ## only mutations within the bed targets to match depth denominator above

mutFreq_prep %>% filter(coding_from_maf == 'non-coding') %>% print(width = Inf)

mutFreq_prep
lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")


mutFreq_prep %>% print(width = Inf)
CHIP_depth_full

### filter for plot
mutFreq_prep <- mutFreq_prep %>%
  left_join(CHIP_depth_full) %>%
  left_join(CHIP_depth_full_final) %>% 
  print(width = Inf) %>%
  group_by(coding, Subject, denominator_coding, age) %>%
  summarise(n_muts = n()) %>%
  mutate(mutFreq = n_muts/denominator_coding) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  print()



## combine genes for plot
mutFreq_prep <- mutFreq_combined %>% 
  group_by(Subject, tissue, coding, age, LFS, CTx) %>%
  summarise(denominator_coding = sum(denominator_coding), 
            n_muts = sum(n_muts),
            mutReads = sum(mutReads)) %>%
  mutate(mutFreq = n_muts / denominator_coding,
         mutBurden = mutReads / denominator_coding)

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
mutFreq_prep
mutFreq_plot2 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep,
              se = FALSE, method = "lm", color = '#444444') +
  # geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
  #             se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_y_log10(limits = c(3.9e-7, 1.1e-6), labels=fancy_scientific) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Overall\nmutation frequency"
  ) +
  theme_minimal()

show(mutFreq_plot2)

###############################################################################
### Mutation burden
###############################################################################

###############################################################################
### Mutation burden plot
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

mutBurden_plot2 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutBurden)) +
  geom_smooth(data = mutFreq_prep,
              se = FALSE, method = "lm", color = '#444444') +
  # geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
  #             se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  scale_y_log10(limits = c(1.5e-6, 3e-5), labels=fancy_scientific) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Overall\nmutation burden"
  ) +
  theme_minimal()

show(mutBurden_plot2)


################################################################################
########### Plot frequency and Burden combined
################################################################################

mutBurden_plot2 <- mutBurden_plot2 + theme(legend.position = "none")
mutFreq_plot2 <- mutFreq_plot2 + theme(legend.position = "none")

legend_shared <- get_legend(
  mutFreq_plot2 + 
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

mutBurden_plot2 <- mutBurden_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  ylab(expression(atop(NA, atop(textstyle("Overall"), textstyle("mutation burden"))))) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_plot2 <- mutFreq_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle("Overall"), textstyle("mutation frequency"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

plots <- plot_grid(
  plot_grid(mutFreq_plot2, mutBurden_plot2, ncol = 2, align = "v"), 
  ncol = 1)

combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

show(combined_plots)
ggsave("results/CHIP_full_mutation_frequency.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
