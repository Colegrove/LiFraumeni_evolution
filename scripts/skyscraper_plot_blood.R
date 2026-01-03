## Li-Fraumeni skyscraper plot


### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
#family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

skyscraper_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(coding == "coding") %>%
  #filter(Subject %in% family) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  
  mutate(
    color_group = case_when(
      am_class == "likely_benign" & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign" & t_alt_count > 1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count > 1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous" & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous" & t_alt_count > 1 ~ "ambiguous_LC"
    )
  ) 

skyscraper_prep2 <- skyscraper_prep %>%
  #group_by(Subject) %>%
  #mutate(n_tissues_for_subject = n_distinct(Tissue)) %>%
  #ungroup() %>%
  mutate(
    x_axis_var = if_else(Subject == "Patient", Tissue, Subject)
  )
skyscraper_prep2 <- skyscraper_prep2 %>%
  mutate(x_axis_var = factor(
    x_axis_var,
    #levels = c("Family member A", "Family member B", "Family member C", "Buffy coat", "Plasma", "Whole blood")
    #levels = c( "Buffy coat", "Family member A", "Family member B", "Family member C", "UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
    levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
    ))

skyscraper_prep2 %>% filter(Subject == "UW volunteer 7") %>% print(width = Inf)

###########################################
##### Alpha missense pathogenicity
###########################################
skyscraper <- skyscraper_prep2 %>%
  ggplot(
    aes(
      x = x_axis_var,
      y = SampCodingOrder,
      fill = color_group,
      label = t_alt_count
    )
  ) +
  geom_tile() +
  
  # uncomment
  #geom_text( size=font.subscript.size*25.4/72.27) +
  geom_text(size=16*25.4/72.27) +
  
  geom_tile(data = skyscraper_prep2,
            aes(
              x = x_axis_var,
              y = SampCodingOrder
            ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  scale_y_continuous(expand = c(0,0),
                     breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                     labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                     limits = c(0,25),
                     name = "Number of <i>TP53</i> mutations") +
  Col.amClass.fill +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8),
        panel.spacing = unit(0,"pt"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.25,0.7),
  )

show(skyscraper)

skyscraper_slide <- skyscraper +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5,size = 20),
    axis.text.y = element_text(size = 20),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
    legend.text = (element_text(size = 20)),
    legend.title = (element_text(size=20))
  ) + geom_text(size=16*25.4/72.27)

show(skyscraper_slide)
ggsave("results/skyscraper_plot_buffycoat_controls.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)

###########################################
##### DMS rfs scores
###########################################
skyscraper_prep2 %>% filter(Subject == "Family member C") %>% print(width = Inf)
DMS_effect <- DMS %>% filter(!is.na(effect))
skyscraper_prep2_merged <- skyscraper_prep2 %>%
  left_join(DMS_effect %>% dplyr::select(effect, rfs_median), by = c("HGVSp_Short" = 'effect'))

x_levels <- unique(skyscraper_prep2$x_axis_var)
skyscraper_prep2_data <- skyscraper_prep2_merged %>%
  filter(!is.na(rfs_median)) %>% 
  group_by(Subject) %>%
  arrange(desc(rfs_median)) %>%
  mutate(SampCodingOrder_rfs = row_number()) %>%
  mutate(x_axis_var = factor(x_axis_var, levels = x_levels))

skyscraper <- skyscraper_prep2_data %>%
  ggplot(
    aes(
      x = x_axis_var,
      y = SampCodingOrder_rfs,
      fill = rfs_median,
      label = t_alt_count
    )
  ) +
  geom_tile() +
  
  # uncomment
  #geom_text( size=font.subscript.size*25.4/72.27) +
  geom_text(size=16*25.4/72.27) +
  
  geom_tile(data = skyscraper_prep2_merged,
            aes(
              x = x_axis_var,
              y = SampCodingOrder
            ), inherit.aes = F, color = "black", fill = NA, linewidth = .1) + 
  scale_y_continuous(expand = c(0,0),
                     breaks=c(0.5,5.5,10.5,15.5,20.5,25.5,30.5,35.5,40.5,45.5,50.5,55.5,60.5),
                     labels=c(0,5,10,15, 20, 25, 30,35,40,45,50,55,60),
                     limits = c(0,25),
                     name = "Number of <i>TP53</i> DBD mutations") +
  scale_fill_gradientn(
    colours = c("blue", "white", "red"),
    limits = c(-1, 1),
    oob = squish,
    name = "RFS Median"
  ) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5,size = 8),
        panel.spacing = unit(0,"pt"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.25,0.7),
  )

show(skyscraper)

skyscraper_slide <- skyscraper +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5,size = 20),
    axis.text.y = element_text(size = 20),
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7),
    legend.text = (element_text(size = 20)),
    legend.title = (element_text(size=20))
  ) + geom_text(size=16*25.4/72.27)

show(skyscraper_slide)
ggsave("results/skyscraper_plot_buffycoat_controls.png", skyscraper_slide, width = 17, height = 13, units = "in", dpi = 300)

#################################
##### Add age and LFS information
#################################

subjects <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
#ages <- c("20-30", "20-30", "20-30", "20-30", "30s", "39", "61", "60s", "69", "70s", "70s")
ages <- c("25", "25", "25", "25", "35", "39", "61", "65", "69", "75", "75")
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
  

annotation_plot <- ggplot(age_lf_long, aes(x = x_axis_var, y = y)) +
  geom_text(aes(label = value), size = 4, vjust = 0) +
  scale_y_continuous(
    breaks = c(1.1, 1.6, 2.1),
    labels = c("CTx", "LFS", "age"),
    #expand = c(0, 0),
    limits = c(1, 2.3)
  ) +
  scale_x_discrete(expand = c(0.07, 0)) +
  theme_minimal(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
show(annotation_plot)


combined_plot <- cowplot::plot_grid(
  skyscraper + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = c(0.2,0.7),),
  annotation_plot,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot)

ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 7, height = 7, units = "in", dpi = 300)

########################### combined for slide

annotation_plot_slide <- ggplot(age_lf_long, aes(x = x_axis_var, y = y)) +
  geom_text(aes(label = value), size = 8, vjust = 0) +
  scale_y_continuous(
    breaks = c(1.1, 1.6, 2.1),
    labels = c("CTx", "LFS", "age"),
    #expand = c(0, 0),
    limits = c(1, 2.3)
  ) +
  scale_x_discrete(expand = c(0.07, 0)) +
  theme_minimal(base_size = 8) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 30),
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0, "pt")
  )
show(annotation_plot_slide)

combined_plot_slide <- cowplot::plot_grid(
  skyscraper_slide + theme(plot.margin = margin(0, 0, 2, 0, "pt"), legend.position = c(0.2,0.7),),
  annotation_plot_slide,
  ncol = 1,
  rel_heights = c(1, 0.15),
  greedy = TRUE
)
show(combined_plot_slide)

ggsave("results/skyscraper_plot_annotated.png", combined_plot_slide, width = 22, height = 10, units = "in", dpi = 300)

#ggsave("results/skyscraper_plot_annotated.png", combined_plot, width = 11, height = 6, units = "in", dpi = 300)


################################################################################
############################ Grouped pathogenicity #############################
################################################################################

########## proportion of pathogenic/benign

pathogenicity <- skyscraper_prep2 %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  filter(!(Subject == "Patient" | Subject == "UW volunteer 7")) %>%
  mutate(pathogenicity_am = case_when(
    color_group %in% c("likely_pathogenic", "likely_pathogenic_LC") ~ "pathogenic",
    color_group %in% c("likely_benign", "likely_benign_LC") ~ "benign",
    TRUE ~ color_group)) %>%
  group_by(LFS, pathogenicity_am) %>%
  summarise(count = n()) %>%
  spread(pathogenicity_am, count, fill = 0)

chi_test <- chisq.test(pathogenicity[,-1]) # exclude LFS column
chi_test

skyscraper_prep2_data %>% print(n=Inf, width = Inf)

mutation_proportions <- pathogenicity %>%
  mutate(total = pathogenic + benign,  # Calculate total mutations for each LFS group
         prop_pathogenic = pathogenic / total,
         prop_benign = benign / total)

p_value = chi_test$p.value

# Plot the proportions of pathogenic mutations for each group
pathogenicity_bar <- ggplot(mutation_proportions, aes(x = LFS, y = prop_pathogenic)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(y = "Proportion of TP53 mutations") +
  theme_minimal() +
  #scale_fill_manual(values = c("LFS" = "red", "non-LFS" = "blue")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.title.x = element_blank()) + 
  annotate("text", x = 1.5, y = 0.85, label = paste("chi-squared p-value = ", round(p_value, 3)), size = 5, color = "black")

mutation_proportions_long <- mutation_proportions %>%
  pivot_longer(cols = c("prop_pathogenic", "prop_benign"), 
               names_to = "mutation_type", 
               values_to = "proportion") %>%
  mutate(mutation_type = recode(mutation_type, 
                                "prop_pathogenic" = "Pathogenic", 
                                "prop_benign" = "Benign"))

# Plot the stacked bar plot
pathogenicity_bar <- ggplot(mutation_proportions_long, aes(x = LFS, y = proportion, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  # Stacked bar plot
  labs(y = "Proportion of TP53 mutations", fill = "Mutation Type") +  # Label the axes and fill legend
  theme_minimal() +
  scale_fill_manual(values = c("Pathogenic" = "red", "Benign" = "blue")) +  # Color for each mutation type
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.title.x = element_blank()) + 
  annotate("text", x = 1.5, y = 1.05, label = paste("p-value = ", round(p_value, 3)), size = 5, color = "black")  # Add p-value text

pathogenicity_bar
ggsave("results/pathogenicity_bar_chi.png", pathogenicity_bar, width = 4, height = 3, units = "in", dpi = 300)


########## pathogenic vs benign large clones

pathogenicity_size <- skyscraper_prep2 %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  #filter(!(Subject == "Patient" | Subject == "UW volunteer 7")) %>%
  mutate(pathogenicity_am = case_when(
    color_group %in% c("likely_pathogenic", "likely_pathogenic_LC") ~ "pathogenic",
    color_group %in% c("likely_benign", "likely_benign_LC") ~ "benign",
    TRUE ~ color_group)) %>%
  filter(pathogenicity_am != "ambiguous") 
  
plot_vaf <- ggplot(pathogenicity_size, aes(x = pathogenicity_am, y = VAF, fill = pathogenicity_am)) +
  geom_jitter(aes(color = pathogenicity_am), width = 0.2, height = 0, size = 2, alpha = 0.7) +  # Set color for dots
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  labs(x = "Mutation Type", y = "Variant Allele Frequency (VAF)") +
  scale_color_manual(values = c("benign" = "blue", "pathogenic" = "red")) +  # Color for dots
  scale_fill_manual(values = c("benign" = "blue", "pathogenic" = "red")) +  # Color for fill in boxplots
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",  # Remove legend from this plot
        axis.title.x = element_blank(),
        axis.title.x = element_blank())  # Optional: remove x-axis label if desired

# Right plot: t_alt_count by pathogenicity_am, with color distinction
plot_t_alt_count <- ggplot(pathogenicity_size, aes(x = pathogenicity_am, y = t_alt_count, fill = pathogenicity_am)) +
  geom_jitter(aes(color = pathogenicity_am), width = 0.2, height = 0, size = 2, alpha = 0.7) +  # Set color for dots
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  labs(x = "Mutation Type", y = "# of duplex reads") +
  scale_color_manual(values = c("benign" = "blue", "pathogenic" = "red")) +  # Color for dots
  scale_fill_manual(values = c("benign" = "blue", "pathogenic" = "red")) +  # Color for fill in boxplots
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        axis.title.x = element_blank())  # Place legend at the top

# Combine the two plots side by side, sharing a single legend
pathogenicity_sizes <- grid.arrange(plot_vaf, plot_t_alt_count, ncol = 2)  # Place legend at the top
ggsave("results/pathogenicity_sizes.png", pathogenicity_sizes, width = 4, height = 3, units = "in", dpi = 300)


###################################################
########### mutation frequencies
###################################################

# age_map <- c("UW volunteer 1" = "25",
#              "UW volunteer 2" = "25",
#              "UW volunteer 3" = "25",
#              "UW volunteer 4" = "25",
#              "Patient" = "35",
#              "Family member A" = "39",
#              "Family member B" = "61",
#              "UW volunteer 5" = "65",
#              "Family member C" = "69",
#              "UW volunteer 6" = "75",
#              "UW volunteer 7" = "75")
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

### depths at each position 
tp53_depth_full <- DP_filt %>%
  filter(Chr == "chr17") %>%
  arrange(Pos) %>%
  filter(Pos < 70672726) %>% # filter out mutagenesis panel on chr17 to only keep TP53
  group_by(Samp) %>%
  summarise(denominator = sum(DP)) %>%
  mutate(denominator_coding = denominator) %>%
  mutate(denominator_noncoding = denominator)

tp53_depth_coding <- DP_filt %>%
  filter(Chr == "chr17") %>%
  arrange(Pos) %>%
  filter(Pos < 70672726) %>% # filter out mutagenesis panel on chr17 to only keep TP53
  filter(coding == "coding") %>%
  print(width = Inf) %>%
  group_by(Samp) %>%
  summarise(denominator_coding = sum(DP))

tp53_depth_noncoding <- DP_filt %>%
  filter(Chr == "chr17") %>%
  arrange(Pos) %>%
  filter(Pos < 70672726) %>% # filter out mutagenesis panel on chr17 to only keep TP53
  filter(coding == "non-coding") %>%
  print(width = Inf) %>%
  group_by(Samp) %>%
  summarise(denominator_noncoding = sum(DP))

tp53_depth_split <- tp53_depth_coding %>%
  left_join(tp53_depth_noncoding)



#### mutation frequency # of mutations/denominator
family_patient_blood_samples <- c("PBMC", "Buffy coat")
mutFreq_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53")


lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")


### filter for plot
mutFreq_prep <- mutFreq_prep %>%
  #left_join(tp53_depth_full) %>%
  left_join(tp53_depth_split) %>%
  print(width = Inf)%>%
  group_by(coding, Subject, denominator_coding, denominator_noncoding, age) %>%
  summarise(n_muts = n()) %>%
  #mutate(mutFreq = n_muts/denominator) %>%
  mutate(mutFreq = if_else(coding == "coding", n_muts/denominator_coding, n_muts/denominator_noncoding)) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  print()

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
      LFS == "LFS" & CTx == "CTx" ~ "LFS+CTx",   # both LFS and CTx
      LFS == "LFS"             ~ "LFS",  # just LFS
      LFS == "non-LFS" & CTx == "CTx" ~ "CTx",
      TRUE                 ~ "non-LFS"      # if you want a default
    )
  )
patient_history_order = c("LFS", "LFS+CTx", "non-LFS", "CTx")
mutFreq_prep$shape_group = factor(mutFreq_prep$shape_group, levels = patient_history_order)

### Mutation frequency plot
mutFreq_plot <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq, color = coding)) +
  geom_line(aes(group = Subject), color = "grey") +
  geom_point(aes(shape=factor(shape_group)),size = 1, alpha = 1, stroke = 1.2) +
  #geom_point(alpha = 1) +
  geom_smooth(aes(group=coding), se = FALSE, method = "lm") +
  scale_y_log10(limits = c(1e-7, 1e-05)) +
  scale_x_continuous(limits = c(20,80), breaks= (seq(20,80,5))) +
  scale_shape_manual(values = c("non-LFS" = 1, "LFS" = 2, "LFS+CTx" = 17, "CTx" = 16)) +
  #scale_y_log10(limits = c(4e-8, 4e-06)) +
  labs(x = "Age", y = "Mutation frequency (# mutations/sequencing depth)", color = "Coding region", shape = "Patient history") +
  theme_minimal()
show(mutFreq_plot)



#### mutation burden # of mutant reads/denominator
mutBurden_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53")

#### filter for plot
mutBurden_prep <- mutBurden_prep %>%
  #left_join(tp53_depth_full) %>%
  left_join(tp53_depth_split) %>%
  group_by(coding, Subject, denominator_coding, denominator_noncoding, age) %>%
  summarise(mutReads = sum(t_alt_count)) %>%
  mutate(mutBurden = if_else(coding == "coding", mutReads/denominator_coding, mutReads/denominator_noncoding)) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx"))
  #mutate(mutBurden = mutReads/denominator)

mutBurden_prep <- mutBurden_prep %>%
  left_join(offsets) %>%
  mutate(
    shape_group = case_when(
      LFS == "LFS" & CTx == "CTx" ~ "LFS+CTx",   # both LFS and CTx
      LFS == "LFS"             ~ "LFS",  # just LFS
      LFS == "non-LFS" & CTx == "CTx" ~ "CTx",
      TRUE                 ~ "non-LFS"      # if you want a default
    )
  )

patient_history_order = c("LFS", "LFS+CTx", "non-LFS", "CTx")
mutBurden_prep$shape_group = factor(mutBurden_prep$shape_group, levels = patient_history_order)


#### Mutation burden Plot
mutBurden_plot <- ggplot(mutBurden_prep, aes(x = age_j, y = mutBurden, color = coding)) +
  geom_line(aes(group = Subject), color = "grey") +
  geom_point(aes(shape=factor(shape_group)),size = 1, alpha = 1, stroke = 1.2) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_y_log10(limits = c(1e-7, 1e-05)) +
  scale_x_continuous(limits = c(20,80), breaks= (seq(20,80,5))) +
  scale_shape_manual(values = c("non-LFS" = 1, "LFS" = 2, "LFS+CTx" = 17, "CTx" = 16)) +
  #scale_y_log10(limits = c(4e-8, 4e-06)) +
  labs(x = "Age", y = "Mutation Burden (# mutant reads/sequencing depth)", color = "Coding Region", shape = "Patient history") +
  theme_minimal()
show(mutBurden_plot)


### combine Frequency and Burden Plots

mutBurden_plot <- mutBurden_plot + theme(legend.position = "none")
mutFreq_plot <- mutFreq_plot + 
    guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
    theme(legend.position  = c(.5,.8), legend.direction = "vertical", legend.box = "horizontal")
combined_plots <- plot_grid(
  #legend, 
  plot_grid(mutFreq_plot, mutBurden_plot, ncol = 2, align = "v"), 
  ncol = 1)
  #rel_heights = c(0.4, 1))
show(combined_plots)  
ggsave("results/tp53_coding_noncoding_mutation_frequency.png", combined_plots, width = 7, height = 4, units = "in", dpi = 300)
#ggsave("results/tp53_full_mutation_frequency.png", combined_plots, width = 5, height = 4, units = "in", dpi = 300)


###################################################
########### protein-affecting
###################################################

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



#### mutation frequency # of mutations/denominator
family_patient_blood_samples <- c("PBMC", "Buffy coat")
mutFreq_protein_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53") %>%
  print(n=Inf)

#mutFreq_protein_prep %>% filter(Subject == "Patient") %>% print(width = Inf)
mutFreq_protein_prep <- mutFreq_protein_prep %>%
  #left_join(tp53_depth_full) %>%
  left_join(tp53_depth_split) %>%
  group_by(Mutation_type, Subject, denominator_coding, denominator_noncoding, age) %>%
  mutate(protein_affecting = coding) %>%
  mutate(protein_affecting = if_else(Mutation_type == 'Silent', "non-coding", protein_affecting)) %>%
  mutate(protein_affecting = if_else(protein_affecting == "coding", "affecting", "non-affecting")) %>%
  group_by(Subject, protein_affecting, denominator_coding, denominator_noncoding, age) %>%
  summarise(n_muts = n()) %>%
  #mutate(mutFreq = if_else(protein_affecting == "affecting", n_muts/denominator_coding, n_muts/denominator_noncoding)) %>%
  mutate(mutFreq = if_else(protein_affecting == "affecting", n_muts/(denominator_coding*2/3), n_muts/(denominator_noncoding + denominator_coding*1/3))) %>%
  print()

mutFreq_protein_plot <- ggplot(mutFreq_protein_prep, aes(x = age, y = mutFreq, color = protein_affecting, group=protein_affecting)) +
  geom_point(alpha = 1) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_y_log10(limits = c(6e-8, 2e-05)) +
  #scale_y_log10(limits = c(4e-8, 4e-06)) +
  labs(x = "Age", y = "Mutation Frequency (# mutations/sequencing depth)", color = "Protein affecting") +
  theme_minimal()
show(mutFreq_protein_plot)

#### mutation burden # of mutant reads/denominator
mutBurden_protein_prep <-
  filt_maf %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53")

mutBurden_protein_prep <- mutBurden_protein_prep %>%
  #left_join(tp53_depth_full) %>%
  left_join(tp53_depth_split) %>%
  group_by(Mutation_type, Subject, denominator_coding, denominator_noncoding, age) %>%
  
  mutate(protein_affecting = coding) %>%
  mutate(protein_affecting = if_else(Mutation_type == 'Silent', "non-coding", protein_affecting)) %>%
  mutate(protein_affecting = if_else(protein_affecting == "coding", "affecting", "non-affecting")) %>%
  group_by(Subject, protein_affecting, denominator_coding, denominator_noncoding, age) %>%
  summarise(mutReads = sum(t_alt_count)) %>%
  #mutate(mutBurden = if_else(protein_affecting == "affecting", mutReads/denominator_coding, mutReads/denominator_noncoding)) %>%
  mutate(mutBurden = if_else(protein_affecting == "affecting", mutReads/(denominator_coding*2/3), mutReads/(denominator_noncoding + denominator_coding*1/3))) %>%
  print()

mutBurden_protein_prep %>% filter(Subject == "Patient")

mutBurden_protein_plot <- ggplot(mutBurden_protein_prep, aes(x = age, y = mutBurden, color = protein_affecting, group=protein_affecting)) +
  geom_point(alpha = 1) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_y_log10(limits = c(6e-8, 2e-05)) +
  #scale_y_log10(limits = c(4e-8, 4e-06)) +
  labs(x = "Age", y = "Mutation Burden (# mutant reads/denominator)", color = "Protein affecting") +
  theme_minimal()
show(mutBurden_protein_plot)

mutFreq_protein_plot <- mutFreq_protein_plot + theme(legend.position = "none")
mutBurden_protein_plot <- mutBurden_protein_plot + theme(legend.position = "none")
legend <- get_legend(mutBurden_protein_plot + theme(legend.position = "top"))
combined_protein_plots <- plot_grid(
  legend, 
  plot_grid(mutFreq_protein_plot, mutBurden_protein_plot, ncol = 2, align = "v"), 
  ncol = 1, 
  rel_heights = c(0.1, 1))
show(combined_protein_plots)  
ggsave("results/tp53_split_denominator_protein_affecting.png", combined_protein_plots, width = 5, height = 4, units = "in", dpi = 300)
#ggsave("results/tp53_same_denominator_protein_affecting.png", combined_protein_plots, width = 5, height = 4, units = "in", dpi = 300)






#####################################
###### non-coding mutation regions (do mutations cluster, especially at indels?)
#####################################

filt_maf %>% filter(is.na(Subject))

filt_maf %>% filter(Hugo_Symbol == "TP53") %>% filter(Subject == "Family member B") %>% filter(coding == "non-coding")

non_coding_TP53 <- filt_maf %>%
  filter(Hugo_Symbol == "TP53") %>%
  #filter(Start_Position >7673000 & Start_Position <7676000) %>%
  arrange(Start_Position) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(coding == "non-coding" ) %>%
  group_by(Subject, coding, Tissue, Start_Position, Variant_Type, t_alt_count) %>%
  summarise(Unique_Mutations = n(), .groups = "drop") %>%
  print(width = Inf, n=Inf)

tp53_bed_targets <- read_delim("./inputs/BEDs/pan00367.targets.hs38DH.bed", col_names = c("Chr", "Start_Position", "End_Position", "Exon"))
tp53_bed_targets

# non_coding_locations <- ggplot(non_coding_TP53, aes(x = Start_Position, y = Subject)) +
#   geom_jitter(height = 0.2, size = 2, alpha = 0.7) +
#   theme_minimal() +
#   xlim(7669000, 7678000) +
#   labs(x = "Position", y = "Subject", title = "Non-Coding Mutation positions") +
#   theme(axis.text.y = element_text(size = 8))
# 
# show(non_coding_locations)
# ggsave("results/tp53_non_coding_locations.png", non_coding_locations, width = 10, height = 4, units = "in", dpi = 300)


tp53_bed_targets <- read_delim("./inputs/BEDs/pan00367.targets.hs38DH.bed", col_names = c("Chr", "x_min", "x_max", "Exon"))
tp53_bed_targets



non_coding_plot <- ggplot() +
  # Exon rectangles in background
  geom_rect(
    data = tp53_bed_targets,
    aes(xmin = x_min, xmax = x_max, ymin = -Inf, ymax = Inf),
    fill = "#66ccee", alpha = .9, 
    inherit.aes = FALSE
  ) +
  # Mutation points per subject
  geom_jitter(
    data = non_coding_TP53,
    aes(x = Start_Position, y = Subject, shape = Variant_Type, color = Variant_Type),
    width = 0, height = 0.2, size = 2, alpha = 0.8
  ) +
  scale_shape_manual(
    values = c(
      SNP = 16,   # filled circle
      DEL = 15,   # filled square
      INS = 17    # filled square
    )
  ) +
  scale_color_manual(
    values = c(
      SNP = "black",
      DEL = "#44aa99",
      INS = "#aa4499"
    )
  ) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), 
        legend.position = "top") +
  labs(x = "Genomic Position") + 
  coord_cartesian(xlim = c(7669000, 7677000))
  #coord_cartesian(xlim = c(7675800, 7676000)) ## family member b
  #coord_cartesian(xlim = c(7673300, 7673600)) ## family member c
  #coord_cartesian(xlim = c(7674000, 7676000)) ## patient
show(non_coding_plot)
ggsave("results/TP53_non_coding_gene_locations.png", non_coding_plot, width = 7, height = 4, units = "in", dpi = 300)




filt_maf %>%
  filter(Hugo_Symbol == "TP53") %>%
  #filter(Subject == "UW volunteer 6") %>%
  filter(Subject == "Patient") %>%
  #filter(Start_Position >7673000 & Start_Position <7676000) %>%
  arrange(Start_Position) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(coding == "non-coding" ) %>%
  #group_by(Subject, coding, Tissue, Start_Position) %>%
  #summarise(Unique_Mutations = n(), .groups = "drop") %>%
  print(width = Inf, n=Inf)

non_coding_TP53 %>% arrange(Start_Position) %>% print(n=Inf)









################################################################################
#################### new annotations skyscraper
################################################################################

library(patchwork)


ann_wide <- tibble::tibble(
  Subject    = c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
                 "UW volunteer 5","UW volunteer 6","UW volunteer 7",
                 "Patient","Family member A","Family member C","Family member B"),
  SampleCode = c("UW01","UW02","UW03","UW04","UW05","UW06","UW07",
                 "LFS01","LFS02","LFS03","REL01"),
  Age        = c(25,25,25,25,65,75,75,35,39,69,61),
  LFS        = c("No","No","No","No","No","No","No","Yes","Yes","Yes","No"),
  CTx        = c("No","No","No","No","No","No","Yes","Yes","No","No","No"),
  Depth      = c("High","High","High","High","High","High","High","High","High","High","High")
) %>%
  mutate(
    LFS   = factor(LFS, levels = c("No","Yes")),
    CTx   = factor(CTx, levels = c("No","Yes")),
    Depth = factor(Depth, levels = c("Low","High"))
  )

# --- your skyscraper data (unchanged prep) ---
family <- c("Family member A", "Family member B", "Family member C")
mstp   <- c("UW volunteer 1","UW volunteer 2","UW volunteer 3","UW volunteer 4",
            "UW volunteer 5","UW volunteer 6","UW volunteer 7")
family_patient_blood         <- c("Family member A","Family member B","Family member C","Patient")
family_patient_blood_samples <- c("PBMC","Buffy coat")

skyscraper_prep <- filt_maf %>%
  filter(Tissue != "Urine cells", coding == "coding") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(
    color_group = case_when(
      am_class == "likely_benign"     & t_alt_count == 1 ~ "likely_benign",
      am_class == "likely_benign"     & t_alt_count >  1 ~ "likely_benign_LC",
      am_class == "likely_pathogenic" & t_alt_count == 1 ~ "likely_pathogenic",
      am_class == "likely_pathogenic" & t_alt_count >  1 ~ "likely_pathogenic_LC",
      am_class == "ambiguous"         & t_alt_count == 1 ~ "ambiguous",
      am_class == "ambiguous"         & t_alt_count >  1 ~ "ambiguous_LC"
    )
  )

skyscraper_prep2 <- skyscraper_prep %>%
  mutate(x_axis_var = if_else(Subject == "Patient", Tissue, Subject))

# --- single shared x-order by Age (Patient tissues inherit Patient's age) ---
axis_table <- skyscraper_prep2 %>%
  distinct(x_axis_var, Subject) %>%
  mutate(subject_key = if_else(Subject == "Patient", "Patient", x_axis_var)) %>%
  left_join(ann_wide, by = c("subject_key" = "Subject")) %>%
  arrange(Age, subject_key)

x_levels <- axis_table$x_axis_var
skyscraper_prep2 <- skyscraper_prep2 %>%
  mutate(x_axis_var = factor(x_axis_var, levels = x_levels))

# --- skyscraper (top) with x labels/ticks removed ---
skyscraper <- skyscraper_prep2 %>%
  ggplot(aes(x = x_axis_var, y = SampCodingOrder, fill = color_group, label = t_alt_count)) +
  geom_tile() +
  geom_text(size = 10*25.4/72.27) +
  geom_tile(aes(x = x_axis_var, y = SampCodingOrder), inherit.aes = FALSE,
            data = skyscraper_prep2, color = "black", fill = NA, linewidth = .1) +
  scale_y_continuous(
    expand = c(0,0),
    breaks = c(0.5,5.5,10.5,15.5,20.5),
    labels = c(0,5,10,15,20),
    limits = c(0,25),
    name = "Number of <i>TP53</i> mutations"
  ) +
  Col.amClass.fill +
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),      # <-- remove labels
    axis.ticks.x = element_blank(),      # <-- remove ticks
    panel.spacing = unit(0,"pt"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.x.top = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=8, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    legend.position = c(0.25,0.7)
  )

skyscraper_slide <- skyscraper +    # keep your slide sizing, but no x labels
  theme(
    axis.text.y = element_text(size = 20),
    strip.text.x.top = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    axis.title.y = element_markdown(family = "sans", size=20, face = "bold", color = "black"),
    legend.text  = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) + geom_text(size=16*25.4/72.27)

# --- annotation rows (bottom x labels = SampleCode; Patient PBMC/Buffy coat both "LFS01") ---
axis_table <- axis_table %>% mutate(x_axis_var = factor(x_axis_var, levels = x_levels))

# label map using SampleCode; duplicate "LFS01" for both patient columns
label_map <- axis_table %>%
  mutate(lbl = if_else(subject_key == "Patient", "LFS01", SampleCode)) %>%
  select(x_axis_var, lbl)

base_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(0, 5.5, 0, 5.5))

p_age <- ggplot(axis_table, aes(x_axis_var, " ", fill = Age)) +
  geom_tile(height = 1) +
  scale_fill_gradient(low = "grey90", high = "grey20") +
  labs(y = "Age", fill = "Age") +
  base_theme + theme(axis.text.x = element_blank(),
                     axis.text.y = element_text(color = "grey20"))

p_lfs <- ggplot(axis_table, aes(x_axis_var, " ", fill = LFS)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = c("Yes" = "grey20", "No" = "grey90")) +
  labs(y = "LFS", fill = "LFS") +
  base_theme + theme(axis.text.x = element_blank(),
                     axis.text.y = element_text(color = "grey20"))

p_ctx <- ggplot(axis_table, aes(x_axis_var, " ", fill = CTx)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = c("Yes" = "grey20", "No" = "grey90")) +
  labs(y = "CTx", fill = "CTx") +
  base_theme + theme(axis.text.x = element_blank(),
                     axis.text.y = element_text(color = "grey20"))

p_depth <- ggplot(axis_table, aes(x_axis_var, " ", fill = Depth)) +
  geom_tile(height = 1) +
  scale_fill_manual(values = c("High" = "black", "Low" = "grey90")) +
  labs(y = "Depth", fill = "Depth") +
  base_theme +
  scale_x_discrete(labels = setNames(label_map$lbl, label_map$x_axis_var)) +  # <-- UW01, LFS02, etc.
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 8),
        axis.text.y = element_text(color = "grey20"))

# --- stack: skyscraper on top, then annotation rows ---
final_plot <- (skyscraper / p_age / p_lfs / p_ctx / p_depth) +
  plot_layout(heights = c(6, 0.9, 0.6, 0.6, 0.8), guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.spacing.y = grid::unit(0.15, "lines")
  )

final_plot





