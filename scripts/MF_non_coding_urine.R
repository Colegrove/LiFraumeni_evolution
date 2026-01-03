# ## Li-Fraumeni mutation frequencies TP53 coding and non-coding


### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
#family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
family_patient_blood_samples <- c("PBMC", "Buffy coat")


################################################################################
############################# Coding vs non-coding #############################
################################################################################

###################################################
########### mutation frequencies
###################################################


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

tp53_depth_coding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  group_by(Samp) %>% 
  summarise(denominator_coding = sum(DP))

tp53_depth_noncoding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(is.na(exon_number)) %>%
  group_by(Samp) %>% 
  summarise(denominator_noncoding = sum(DP))

tp53_depth_split <- tp53_depth_coding %>%
  left_join(tp53_depth_noncoding)
tp53_depth_split

#### mutation frequency # of mutations/denominator
family_patient_blood_samples <- c("PBMC", "Buffy coat")
mutFreq_prep <- 
  maf_masked_coding %>%
  #filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  #filter(Tissue %in% c("Urine cells")) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(age = age_map[Subject]) %>%
  mutate(plot_coding = if_else(!is.na(exon_number), "coding", "non-coding"))
mutFreq_prep %>% print(width = Inf)

mutFreq_prep
lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")
### filter for plot
mutFreq_prep <- mutFreq_prep %>%
  left_join(tp53_depth_split) %>%
  #group_by(plot_coding, Subject, Tissue, denominator_coding, denominator_noncoding, age, am_class) %>%
  group_by(plot_coding, Subject, Tissue, denominator_coding, denominator_noncoding, age) %>%
  summarise(n_muts = n()) %>%
  mutate(mutFreq = if_else(plot_coding == "coding", n_muts/denominator_coding, n_muts/denominator_noncoding)) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  print()


mutFreq_prep %>% filter(Tissue ==)
######################################
samplecodes <- ann_wide %>% dplyr::select(Subject, SampleCode)
mutFreq_prep
mutFreq_prep_bar <- mutFreq_prep %>%
  left_join(samplecodes) %>%
  filter(plot_coding == "coding") %>%
  # Create a combined label for clarity
  mutate(TissueLabel = paste(Tissue, SampleCode, sep = " – ")) %>%
  mutate(TissueLabel = fct_reorder(TissueLabel, mutFreq, .fun = identity, .desc = FALSE))

tissue_MFs <- ggplot(mutFreq_prep_bar, aes(x = reorder(TissueLabel, -mutFreq), y = mutFreq)) +
  geom_col(width = 0.7, color = "black") +
  labs(
    y = "MF"
  ) +
  #scale_y_log10() +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    legend.position = "none"
  )
ggsave("results/tp53_tissue_coding_MFs.png", tissue_MFs, width = 5, height = 3, units = "in", dpi = 300)

mutFreq_prep
######################################


mutFreq_prep %>% arrange(desc(mutFreq))
offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset) %>%
  dplyr::select(-n_muts, -mutFreq)

mutFreq_prep <- mutFreq_prep %>%
  left_join(offsets) %>%
  print() %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS"     & CTx == "CTx"    ~ "LFS/CTx",
      LFS == "LFS"     & CTx == "non-CTx" ~ "LFS/no-CTx",
      LFS == "non-LFS" & CTx == "CTx"    ~ "non-LFS/CTx",
      LFS == "non-LFS" & CTx == "non-CTx" ~ "non-LFS/no-CTx"
    )
  )
patient_history_order <- c("LFS/no-CTx", "LFS/CTx", "non-LFS/no-CTx", "non-LFS/CTx")
mutFreq_prep$shape_group <- factor(mutFreq_prep$shape_group, levels = patient_history_order)
###############################################################################
### Mutation frequency plot
###############################################################################

coding_shapes <- c("coding" = 17, "non-coding" = 1)  # triangle vs circle

# Colors for patient history
shape_group_colors <- c(
  "non-LFS/no-CTx" = "#44AA99",
  "non-LFS/CTx"    = "#44AA99",
  "LFS/no-CTx"     = "#882255",
  "LFS/CTx"        = "#882255"
)


# Shapes for region: circle vs triangle (both fillable)
coding_shapes <- c("coding" = 24, "non-coding" = 21)  # triangle vs circle

LFS_colors <- c("LFS" = "#882255", "non-LFS" = "#44AA99")




mutFreq_prep
mutFreq_prep <- mutFreq_prep %>% filter(am_class == "likely_pathogenic") %>% filter(plot_coding == "coding")


mutFreq_prep
## version 3 no vertical lines and coding/non-coding regression
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
  scale_y_log10(limits = c(2e-7, 1e-5)) +
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
ggsave("results/tp53_urine_non_coding.png", mutFreq_plot3, width = 3, height = 2.5, units = "in", dpi = 300)









