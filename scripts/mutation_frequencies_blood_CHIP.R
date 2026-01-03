# ## Li-Fraumeni mutation frequencies


### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")

family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
#family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

CHIP_genes <- c("NRAS", "BRINP3", "DNMT3A", "IDH1", "GATA2", "KIT", "TET2", "NPM1", 
                "EZH2", "RAD21", "HNRNPK", "PTEN", "SMC3", "WT1", "KMT2A", "CBL", "KRAS", 
                "PTPN11", "FLT3", "IDH2", "MYH11", "CEBPA", "ASXL1", "RUNX1", "U2AF1", "SMC1A", 
                "STAG2", "PHF6", "TP53")


###################################################
########### Pull sequencing depths in coding regions
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


### depths at each position 
CHIP_depth_full <- DP_filt_CHIP %>%
  filter(coding == "coding") %>% ## includes only .bed file targets (including targeted introns)
  arrange(Pos) %>%
  group_by(Samp) %>% 
  summarise(denominator = sum(DP)) %>%
  mutate(denominator_coding = denominator)
CHIP_depth_full %>% print(n=Inf)

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
  print(width = Inf) %>%
  group_by(coding, Subject, denominator_coding, age) %>%
  summarise(n_muts = n()) %>%
  mutate(mutFreq = n_muts/denominator_coding) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx")) %>%
  print()

## offsets are used for plotting manual jitter of individuals with same age
offsets <- mutFreq_prep %>%
  ungroup() %>%
  group_by(Subject) %>%
  mutate(age_offset = runif(1,-2.1,2.1), 
         age_j = age + age_offset) %>%
  dplyr::select(-n_muts, -mutFreq)
offsets

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

## color by LFS status, regression is generated using only non-LFS 
mutFreq_plot <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq, color = LFS)) +
  #geom_line(aes(group = Subject), color = "grey") +
  geom_point(aes(shape=factor(shape_group)),size = 3, alpha = 1, stroke = 1.2) +
  geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"), aes(group=LFS), se = FALSE, method = "lm") +
  scale_x_continuous(limits = c(20,80), breaks= (seq(20,80,5))) +
  scale_y_log10(labels=fancy_scientific) +
  #scale_y_log10(limits = c(3e-7, 3e-05)) +
  scale_shape_manual(values = c("non-LFS" = 1, "LFS" = 2, "LFS+CTx" = 17, "CTx" = 16)) +
  labs(x = "Age", y = "Mutation frequency (# mutations/sequencing depth)", color = "LFS status") +
  theme_minimal()
show(mutFreq_plot)

LFS_color <- "#882255"
nonLFS_color <- "#44aa99"

# Define colors by LFS status
shape_group_colors <- c(
  "non-LFS" = nonLFS_color,  # non-LFS color
  "CTx"     = nonLFS_color,  # same as non-LFS
  "LFS"     = LFS_color,  # LFS color
  "LFS+CTx" = LFS_color   # same as LFS
)

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS" = 1,
  "CTx"     = 16,
  #"LFS"     = 2,
  #"LFS+CTx" = 17
  "LFS"     = 1,
  "LFS+CTx" = 16
)

mutFreq_plot2 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq)) +
  geom_point(aes(color = shape_group, shape = shape_group), size = 3, alpha = 1, stroke = 1.2) +
  geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
              se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  scale_y_log10(labels=fancy_scientific) +
  #scale_y_log10(limits = c(3e-7, 3e-5)) +
  #scale_y_log10(limits = c(4e-7, 2e-6)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Mutation frequency (# mutations / sequencing depth)"
  ) +
  theme_minimal()

show(mutFreq_plot2)

###############################################################################
### Mutation burden
###############################################################################

#### mutation burden # of mutant reads/denominator
mutBurden_prep <-
  filt_maf_CHIP %>%
  filter(Tissue != "Urine cells") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  filter(InBed == TRUE) %>% ## only mutations within the bed targets to match depth denominator above
  print(width = Inf)

mutFreq_prep
mutBurden_prep

#### filter for plot
mutBurden_prep <- mutBurden_prep %>%
  left_join(CHIP_depth_full) %>%
  group_by(coding, Subject, denominator_coding, age) %>%
  summarise(mutReads = sum(t_alt_count)) %>%
  mutate(mutBurden = mutReads/denominator_coding) %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "non-CTx"))

mutBurden_prep
mutFreq_prep

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


###############################################################################
### Mutation burden plot
###############################################################################


mutBurden_plot <- ggplot(mutBurden_prep, aes(x = age_j, y = mutBurden, color = LFS)) +
  #geom_line(aes(group = Subject), color = "grey") +
  geom_point(aes(shape=factor(shape_group)),size = 3, alpha = 1, stroke = 1.2) +
  geom_smooth(data = mutBurden_prep %>% filter(LFS == "non-LFS"),se = FALSE, method = "lm") +
  scale_y_log10(labels=fancy_scientific) +
  #scale_y_log10(limits = c(1e-7, 1e-05)) +
  #scale_y_log10(limits = c(3e-7, 3e-05)) +
  scale_x_continuous(limits = c(20,80), breaks= (seq(20,80,5))) +
  scale_shape_manual(values = c("non-LFS" = 1, "LFS" = 2, "LFS+CTx" = 17, "CTx" = 16)) +
  labs(x = "Age", y = "Mutation Burden (# mutant reads/sequencing depth)", color = "LFS status") +
  theme_minimal()
show(mutBurden_plot)

# Define colors by LFS status
shape_group_colors <- c(
  "non-LFS" = nonLFS_color,  # non-LFS color
  "CTx"     = nonLFS_color,  # same as non-LFS
  "LFS"     = LFS_color,  # LFS color
  "LFS+CTx" = LFS_color   # same as LFS
)

# Define shapes by patient history
shape_group_shapes <- c(
  "non-LFS" = 1,
  "CTx"     = 16,
  #"LFS"     = 2,
  #"LFS+CTx" = 17
  "LFS"     = 1,
  "LFS+CTx" = 16
)

mutBurden_plot2 <- ggplot(mutBurden_prep, aes(x = age_j, y = mutBurden)) +
  geom_point(aes(color = shape_group, shape = shape_group), size = 3, alpha = 1, stroke = 1.2) +
  geom_smooth(data = mutBurden_prep %>% filter(LFS == "non-LFS"),
              se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
  scale_y_log10(labels=fancy_scientific) +
  #scale_y_log10(limits = c(3e-7, 3e-5)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Mutation Burden (# mutant reads / sequencing depth)"
  ) +
  theme_minimal()

show(mutBurden_plot2)



################################################################################
########### Plot frequency and Burden combined
################################################################################

mutBurden_plot2 <- mutBurden_plot2 + theme(legend.position = "none")
mutFreq_plot2 <- mutFreq_plot2 + 
    guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
    theme(legend.position  = c(.3,.85), legend.direction = "vertical", legend.box = "horizontal", legend.title = element_blank(),
          legend.background = element_rect(color = "black"))
combined_plots <- plot_grid(
  #legend, 
  plot_grid(mutFreq_plot2, mutBurden_plot2, ncol = 2, align = "v"), 
  ncol = 1)
  #rel_heights = c(0.4, 1))
show(combined_plots)  
ggsave("results/CHIP_full_mutation_frequency.png", combined_plots, width = 6, height = 5, units = "in", dpi = 300)
ggsave("results/CHIP_full_mutation_frequency.png", combined_plots, width = 8, height = 4, units = "in", dpi = 300)


################################################################################
########### protein-affecting
################################################################################

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
  filt_maf_CHIP %>%
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
  filt_maf_CHIP %>%
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

non_coding_CHIP <- filt_maf_CHIP %>%
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




CHIP_genes <- filt_maf_CHIP %>%
  #filter(Hugo_Symbol == "TP53") %>%
  #filter(Subject == "UW volunteer 6") %>%
  #filter(Subject == "Patient") %>%
  #filter(Start_Position >7673000 & Start_Position <7676000) %>%
  arrange(Start_Position) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(coding == "coding" ) %>%
  #group_by(Subject, coding, Tissue, Start_Position) %>%
  #summarise(Unique_Mutations = n(), .groups = "drop") %>%
  pull(Hugo_Symbol)
unique(CHIP_genes)




non_coding_TP53 %>% arrange(Start_Position) %>% print(n=Inf)
