# ## Li-Fraumeni mutation frequencies and burden
## TP53 only


### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
family_patient_blood_samples <- c("PBMC", "Buffy coat")


################################################################################
############################ Coding LFS vs non-LFS #############################
################################################################################

###################################################
########### mutation frequencies
###################################################

## function to turn "e-6" to "10^-6"
fancy_scientific <- function(l) { 
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  parse(text = l)
}

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


tp53_depth_coding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(!is.na(exon_number)) %>%
  filter(!inRepeatMask | is.na(inRepeatMask)) %>%
  group_by(Samp) %>% 
  summarise(denominator_coding = sum(DP))

tp53_depth_noncoding <- final_masked_depth %>%
  filter(gene_name == "TP53") %>%
  filter(is.na(exon_number)) %>%
  filter(!inRepeatMask) %>%
  group_by(Samp) %>% 
  summarise(denominator_noncoding = sum(DP))



tp53_depth_split <- tp53_depth_coding %>%
  left_join(tp53_depth_noncoding)


#### mutation frequency # of mutations/denominator
family_patient_blood_samples <- c("PBMC", "Buffy coat")

mutFreq_counts <- 
  maf_masked_coding %>%
  filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53") %>%
  group_by(Tumor_Sample_Barcode, age, Subject) %>%
  summarise(
    n_muts   = n(),
    mutReads = sum(t_alt_count),
    .groups  = "drop"
  )
mutFreq_counts


lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")

mutFreq_counts
tp53_depth_split
### filter for plot
mutFreq_prep <- mutFreq_counts %>%
  left_join(tp53_depth_split, by= c("Tumor_Sample_Barcode"="Samp")) %>%
  mutate(mutFreq = n_muts/denominator_coding) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) 


offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset) %>%
  dplyr::select(-n_muts, -mutFreq)

mutFreq_prep_TP53 <- mutFreq_prep %>%
  left_join(offsets) %>%
  print() %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS" & CTx == "CTx" ~ "LFS/CTx",   # both LFS and CTx
      LFS == "LFS"             ~ "LFS/no-CTx",  # just LFS
      LFS == "non-LFS" & CTx == "CTx" ~ "non-LFS/CTx",
      TRUE                 ~ "non-LFS/no-CTx"      # if you want a default
    )
  )

patient_history_order = c("LFS/no-CTx", "LFS/CTx", "non-LFS/no-CTx", "non-LFS/CTx")
mutFreq_prep_TP53$shape_group = factor(mutFreq_prep_TP53$shape_group, levels = patient_history_order)

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

lm_model <- lm(mutFreq ~ age, data = mutFreq_prep)
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]

mutFreq_plot2 <- ggplot(mutFreq_prep_TP53, aes(x = age_j, y = mutFreq)) +
  geom_smooth(data = mutFreq_prep_TP53,
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  # geom_smooth(data = mutFreq_prep_TP53 %>% filter(LFS == "non-LFS" & coding == "coding"),
  #             se = FALSE, method = "lm", color = "#44AA99") +  # match non-LFS color
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(labels=fancy_scientific, limits = c(3e-7,3e-06)) +
  scale_y_log10() +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Mutation frequency"
  ) +
  theme_minimal()

show(mutFreq_plot2)


###################################################
########### mutation burden
###################################################
mutFreq_prep_TP53 <- mutFreq_prep_TP53 %>%
  mutate(mutBurden = mutReads/denominator_coding)

mutBurden_plot2 <- ggplot(mutFreq_prep_TP53, aes(x = age_j, y = mutBurden)) +
  geom_smooth(data = mutFreq_prep_TP53,
              se = FALSE, method = "lm", color = '#444444') +
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  # geom_smooth(data = mutBurden_prep_TP53 %>% filter(LFS == "non-LFS" & coding == "coding"),
  #             se = FALSE, method = "lm", color = "#44AA99") +  # match non-LFS color
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  #scale_y_log10(labels=fancy_scientific, limits = c(3e-7,1e-05)) +
  scale_y_log10() +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Mutation burden"
  ) +
  theme_minimal()

show(mutBurden_plot2)


### combine Frequency and Burden Plots

legend_shared <- get_legend(
  mutFreq_plot2 + 
    guides(
      color = guide_legend(nrow = 2, title.position = "top"),  # force 1 row
      shape = guide_legend(nrow = 2, title.position = "top")   # force 1 row
    ) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.text  = element_text(size = 8, margin=margin(r=3)),
          legend.margin = margin(r=1, l=1, b = -15, t = -25),
          legend.box.margin = margin(r=0, l=0, b = 0, t = 5),
          legend.key.size   = unit(8, "pt"),
          legend.spacing.x  = unit(0.2, "cm"),
          legend.spacing.y  = unit(0.2, "cm")
          #legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3)
    )
)

mutBurden_plot2 <- mutBurden_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))
mutFreq_plot2 <- mutFreq_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  theme(legend.position  = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8)
  )


mutBurden_plot2 <- mutBurden_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  ylab(expression(atop(NA, atop(textstyle(italic("TP53")), textstyle("mutation burden"))))) +
  theme(legend.position = "none",
        text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        #axis.title.y = element_text(size = 8, vjust = 0.2, hjust = 0.5),
        axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

mutFreq_plot2 <- mutFreq_plot2 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
  ylab(expression(atop(NA, atop(textstyle(italic("TP53")), textstyle("mutation frequency"))))) +
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
combined_plots
ggsave("results/MF_MB_TP53_supplement.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
