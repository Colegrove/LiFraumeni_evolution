# ## Li-Fraumeni mutation frequencies all mutagenesis regions


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
########### Pull sequencing depths in mutagenesis regions
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
final_masked_depth %>% print(width = Inf)

mutagenesis_depths <- final_masked_depth %>%
  filter(str_starts(Gene, "region")) %>% 
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(denominator_mutagenesis = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

chip_depths_non_coding <- final_masked_depth %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number)) %>%
  left_join(sample_map %>% dplyr::select(sample, subject, tissue), by = c("Samp" = "sample")) %>%
  group_by(Samp, subject, tissue) %>%
  summarise(denominator_genic = sum(DP), .groups = "drop") %>% 
  filter(tissue %in% family_patient_blood_samples)

non_coding_depth <- mutagenesis_depths %>%
  left_join(chip_depths_non_coding) %>%
  mutate(total_non_coding_denominator = denominator_mutagenesis + denominator_genic)
non_coding_depth


# ## mutFreq by gene
# mutFreq_counts <- 
#   maf_masked_coding %>%
#   filter(!is.na(exon_number) | Variant_Classification == "Splice_Site") %>%
#   filter(!inRepeatMask | Variant_Classification == "Splice_Site") %>%
#   filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
#   filter(Tissue %in% family_patient_blood_samples) %>%
#   mutate(age = age_map[Subject]) %>%
#   filter(Hugo_Symbol %in% CHIP_genes) %>%
#   group_by(Tumor_Sample_Barcode, Hugo_Symbol, age) %>% 
#   summarise(
#     n_muts   = n(),
#     mutReads = sum(t_alt_count),
#     .groups  = "drop"
#   )
# 
# 
# mutFreq_combined <- chip_depths %>%
#   left_join(mutFreq_counts, 
#             by = c("Samp" = "Tumor_Sample_Barcode",
#                    "gene_name" = "Hugo_Symbol")) %>%
#   mutate(
#     n_muts   = replace_na(n_muts, 0),
#     mutReads = replace_na(mutReads, 0),
#     mutFreq   = n_muts   / denominator_coding,
#     mutBurden = mutReads / denominator_coding,
#     LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
#     CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
#     #coding = replace_na(coding, "coding"), 
#     age = age_map[subject]) %>%
#     dplyr::rename("Subject" = 'subject', "Hugo_Symbol" = "gene_name")
# 
# ## write mutation frequencies by gene/subject
# write.csv(mutFreq_combined, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene.csv")

## mutations frequency non_coding genic
mutFreq_counts_non_coding_genic <- 
  maf_masked_coding %>%
  filter(!is.na(gene_name)) %>%
  filter(is.na(exon_number) & Variant_Classification != "Splice_Site") %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  filter(!Variant_Classification == "Splice_Site") %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol %in% CHIP_genes) %>%
  group_by(Subject, Tissue) %>%
  summarise(
    n_muts_genic   = n(),
    mutReads_genic = sum(t_alt_count),
    .groups  = "drop"
  )
testBed_MUT_add <- testBed_MUT %>%
  dplyr::rename(MUT_region = Gene)

maf_masked_coding_mut <- maf_masked_coding %>%
  left_join(testBed_MUT_add %>% dplyr::select(-InBed, -inMask), by = c("Chromosome" = "Chr", "Start_Position"= "Pos")) %>%
  left_join(testBed_MUT_add%>% dplyr::select(-InBed, -inMask), by = c("Chromosome" = "Chr", "End_Position"= "Pos"), 
            suffix=c("_StartPosition","_EndPosition")) %>% print(width = Inf)

## mutation frequency mutagenesis
mutFreq_counts_non_coding_mutagenesis <- 
  maf_masked_coding_mut %>%
  filter(!is.na(MUT_region_StartPosition) | !is.na(MUT_region_EndPosition)) %>%
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  group_by(Subject, Tissue) %>%
  summarise(
    n_muts_MUT   = n(),
    mutReads_MUT = sum(t_alt_count),
    .groups  = "drop"
  )
mutFreq_counts_non_coding_mutagenesis

non_coding_depth
mutFreq_counts_non_coding_mutagenesis
mutFreq_counts_non_coding_genic
mutFreq_combined_non_coding <- non_coding_depth %>%
  left_join(mutFreq_counts_non_coding_mutagenesis, by = c("subject" = "Subject", "tissue" = "Tissue")) %>%
  left_join(mutFreq_counts_non_coding_genic, by = c("subject" = "Subject", "tissue" = "Tissue")) %>% 
  mutate(
    n_muts_total = n_muts_MUT + n_muts_genic,
    mutFreq_MUT   = n_muts_MUT   / denominator_mutagenesis,
    mutFreq_Total   = n_muts_total   / total_non_coding_denominator,
    #mutBurden = mutReads / denominator_coding,
    LFS = if_else(subject %in% lfs_subjects, "LFS", "non-LFS"),
    CTx = if_else(subject %in% ctx_subjects, "CTx", "non-CTx"), 
    #coding = replace_na(coding, "coding"), 
    age = age_map[subject]) %>%
  dplyr::rename("Subject" = 'subject')
mutFreq_combined_non_coding
## write mutation frequencies by gene/subject
write.csv(mutFreq_combined_non_coding, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/mutagenesis_muts_by_gene_non_coding.csv")


###############################################################################
### Mutation frequency
###############################################################################


mutFreq_combined_non_coding
mutFreq_prep <- mutFreq_combined_non_coding 
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
is.function(fancy_scientific)
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
lm_model <- lm(mutFreq_MUT ~ age, data = mutFreq_prep)
model_summary <- summary(lm_model)
r2 <- model_summary$r.squared
pval <- coef(model_summary)[2, 4]
r2
pval

mutFreq_plot2 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq_MUT)) +
  geom_smooth(data = mutFreq_prep,
              se = FALSE, method = "lm", color = '#444444') +
  # geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
  #             se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  #scale_y_log10(limits = c(1e-9, 4e-5), labels=as_labeller(function(x) fancy_scientific(x))) +
  scale_y_log10() +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "Mutagenesis\nmutation frequency"
  ) +
  theme_minimal() 

show(mutFreq_plot2)

###############################################################################
### Mutation burden plot
###############################################################################
mutFreq_prep %>% print(width = Inf)
mutFreq_plot3 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq_Total)) +
  geom_smooth(data = mutFreq_prep,
              se = FALSE, method = "lm", color = '#444444') +
  # geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
  #             se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
  geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
  #scale_y_log10(limits = c(1e-9, 4e-5), labels=as_labeller(function(x) fancy_scientific(x))) +
  scale_y_log10() +
  scale_color_manual(values = shape_group_colors, name = "Patient history") +
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  labs(
    x = "Age",
    y = "All non-coding\nmutation frequency"
  ) +
  theme_minimal() 

show(mutFreq_plot3)


################################################################################
########### Plot frequency and frequency
################################################################################

mutFreq_plot2 <- mutFreq_plot2 + theme(legend.position = "none")
mutFreq_plot3 <- mutFreq_plot3 + theme(legend.position = "none")

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

mutFreq_plot3 <- mutFreq_plot3 + 
  scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
  ylab(expression(atop(NA, atop(textstyle("All non-coding"), textstyle("mutation frequency"))))) +
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
  ylab(expression(atop(NA, atop(textstyle("Mutagenesis"), textstyle("mutation frequency"))))) +
  theme(legend.position  = "none",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
        axis.text  = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8))

plots <- plot_grid(
  plot_grid(mutFreq_plot2, mutFreq_plot3, ncol = 2, align = "v"), 
  ncol = 1)

combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
                            rel_heights = c(1.1,0.225))

show(combined_plots)
ggsave("results/non_coding_full_mutation_frequency.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)




non_coding_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene_non_coding.csv"
path = non_coding_path

all_counts <- tibble(read.csv(path)) %>% print()

  
mf_ratio_prep <- all_counts %>%
dplyr::select(Subject, LFS, Hugo_Symbol, coding, mutFreq, age) %>%
# wide: one row per Subject/gene, columns for coding & non-coding-CHIP
pivot_wider(
  names_from  = coding,
  values_from = mutFreq
) %>%
# compute ratio: coding / non-coding-CHIP
mutate(
  MF_ratio = coding / `non-coding-CHIP`
) %>% filter(Hugo_Symbol == "DNMT3A")

mf_ratio_prep
mf_ratio_plot <- ggplot(mf_ratio_prep,
                        aes(x = age, y = MF_ratio, color = LFS)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.9) +
  #facet_wrap(~ Hugo_Symbol) +
  #scale_y_log10() +  # optional but usually helpful for ratios
  scale_color_manual(values = c("non-LFS" = "#44AA99", "LFS" = "#882255")) +
  labs(
    x = NULL,
    y = "coding MF / non-coding-CHIP MF",
    color = "LFS status"
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.position = "top",
    legend.title    = element_blank()
  )

mf_ratio_plot



                                ###########################################################################################
                        ###############################################################################
                                  ########################################################################################################
                        ###########################################################################################
                        #####################################################################################################################
                              ##############################################################################
                        ###########################################################################################
                                      ##############################################################################
                              ###########################################################################################



coding = FALSE # true for coding false for non-coding

non_coding_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/mutagenesis_muts_by_gene_non_coding.csv"

path = non_coding_path


mut_counts <- tibble(read.csv(path)) %>% dplyr::select(-X) %>%
  mutate(CTx = ifelse(is.element(Subject, c("Patient", "UW volunteer 7")), 1, 0), 
         LFS = ifelse(is.element(Subject, c("Patient", "Family member A",  "Family member C")), 1, 0))

mut_counts %>% print(width = Inf)

depth_adjust = 10000000
mut_dat <- mut_counts %>% 
  mutate(decades = age / 10, 
         depth_MUT = denominator_mutagenesis/depth_adjust,
         depth_total = total_non_coding_denominator/depth_adjust) 

mut_dat <- mut_dat %>%
  mutate(MF_MUT = n_muts_MUT/depth_MUT,
         MF_total = n_muts_total/depth_total)
mut_dat %>% print(width = Inf)

#with and without CTx status
models <- c(#"MF_MUT ~ 1 + decades + depth + LFS",
            #"MF_MUT ~ 1 + decades + depth + LFS + CTx",
            "MF_MUT ~ 1 + decades + LFS + CTx",
            #"MF_total ~ 1 + decades + depth + LFS",
            #"MF_total ~ 1 + decades + depth + LFS + CTx",
            "MF_total ~ 1 + decades + LFS + CTx")

mut_dat %>% print(width = Inf)
allRegressions <- 
  foreach(familyToCheck = c("gaussian", "poisson"), .combine = "rbind")%do%{
    foreach(model = models, .combine = "rbind")%do%{
      foreach(excludeCON07 = c("TRUE", "FALSE"), .combine = "rbind")%do%{
          mod <- summary(glm(model, data = dat, family = familyToCheck))
          
          call <- paste(mod$call)
          coef <- mod$coefficients
          #This is just so we can standardize across Gaus and Pois links                   
          colnames(coef) <- c("Estimate", "Std. Error", "Stat", "p")
          regressionResults <- bind_cols(model_type = call[1], call = model, 
                                         model = call[3], 
                                         family = familyToCheck,
                                         noCON07 = excludeCON07,
                                         #gene = geneToCheck,
                                         variables = rownames(coef), 
                                         coef) 
        
      }
    }
  }

cols <- c( "#389688",  "#731433" )
names(cols) <- c( 0, 1)
allRegressions

# # chip_dat_exclude <- chip_dat %>% left_join(ave_counts) %>%
# #   filter(ave_muts_per_sample >= 2)
# 
# freq_depth <- mut_dat %>% 
#   ggplot(aes(x = depth, y = MF)) +
#   #ggplot(aes(x = depth, y = n_muts)) + 
#   geom_point(aes(color = paste(LFS))) +
#   facet_wrap(~factor(gene, levels = geneList), scales = "free") +
#   geom_smooth(method='lm') +
#   #labs(x = "Reads/base (in 10000000s)", y = "Number of unique mutations", col = "LFS") +
#   labs(x = "Reads/base (in 10000000s)", y = "MF", col = "LFS") +
#   theme_classic() + 
#   scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# freq_age <- mut_dat %>% 
#   ggplot(aes(x = decades, y = MF)) +
#   #ggplot(aes(x = decades, y = n_muts)) + 
#   geom_point(aes(color = paste(LFS))) +
#   facet_wrap(~factor(gene, levels = geneList), scales = "free") +
#   geom_smooth(method='lm') +
#   #labs(x = "Age (in decades)", y = "Number of unique mutations", col = "LFS") +
#   labs(x = "Age (in decades)", y = "MF", col = "LFS") +
#   theme_classic() + 
#   scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# freq_age
# height <- 10
# width <- 12
# 
# if(coding){
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_depth.pdf', freq_depth, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age.pdf', freq_age, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_depth.pdf', burden_depth, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_age.pdf', burden_age, height = height, width = width)
# } else {
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_depth_non_coding.pdf', freq_depth, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age_non_coding.pdf', freq_age, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_depth_non_coding.pdf', burden_depth, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_age_non_coding.pdf', burden_age, height = height, width = width)
# }

### plot for supplement
# supp_theme <- function() {
#   theme(
#     axis.title.x = element_text(size = 8),
#     axis.text.x  = element_text(size = 8, angle = 0, vjust = 0, hjust = 0.5),
#     axis.title.y = element_text(size = 8),
#     axis.text.y  = element_text(size = 8),
#     strip.text   = element_text(size = 8),
#     legend.text  = element_text(size = 8),
#     legend.title = element_text(size = 8),
#     legend.position = "top",
#     legend.box.margin = margin(-8, 0, -8, 0),
#     legend.margin = margin(0, 0, 0, 0)
#   )
# }
# 
# freq_age_supp <- freq_age + supp_theme()
# freq_depth_supp <- freq_depth + supp_theme()
# burden_depth_supp <- burden_depth + supp_theme()
# burden_age_supp <- burden_age + supp_theme()
# 
# height <- 6
# width <- 6
# if(coding){
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_depth.png', freq_depth_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age.png', freq_age_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_depth.png', burden_depth_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_age.png', burden_age_supp, height = height, width = width)
# } else {
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_depth_non_coding.png', freq_depth_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age_non_coding.png', freq_age_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_depth_non_coding.png', burden_depth_supp, height = height, width = width)
#   ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_age_non_coding.png', burden_age_supp, height = height, width = width)
# }


## Ok, let's start by just looking at Gaussian, all individuals included no CTx
base_frequency_MUT <- allRegressions %>% filter(family == "gaussian", 
                                            #call == "n_muts ~ 1 + decades + depth + LFS", 
                                            call == "MF_MUT ~ 1 + decades + LFS + CTx",
                                            noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest') 



base_frequency_total <- allRegressions %>% filter(family == "gaussian", 
                                         #call == "n_mutReads ~ 1 + decades + depth + LFS", 
                                         call == "MF_total ~ 1 + decades + LFS + CTx",
                                         noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest')

base_frequency_total


write.csv( base_frequency, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_non_coding.csv", row.names = FALSE )
write.csv( base_burden, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_non_coding.csv", row.names = FALSE )
write.csv( base_frequency_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_non_coding_CTx.csv", row.names = FALSE )
write.csv( base_burden_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_non_coding_CTx.csv", row.names = FALSE )

