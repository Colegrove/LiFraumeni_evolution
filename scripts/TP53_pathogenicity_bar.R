## Li-Fraumeni tp53 pathogenic comparison vs benign and LFS status VAF  


### family member and control samples
family <- c("Family member A", "Family member B", "Family member C")
mstp <- c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "UW volunteer 5", "UW volunteer 6", "UW volunteer 7")
lfs_subjects <- c("Patient", "Family member A", "Family member C")
ctx_subjects <- c("Patient", 'UW volunteer 7')
family_patient_blood <- c("Family member A", "Family member B", "Family member C", "Patient")
#family_patient_blood_samples <- c("PBMC", "Buffy coat", "Plasma", "Whole blood")
family_patient_blood_samples <- c("PBMC", "Buffy coat")

maf_masked_coding %>% print(width = Inf)
skyscraper_prep <-
  maf_masked_coding %>%
  filter(Hugo_Symbol == "TP53") %>%
  filter(!is.na(am_pathogenicity)) %>%
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
skyscraper_prep$am_pathogenicity
skyscraper_prep2 <- skyscraper_prep %>%
  mutate(
    x_axis_var = if_else(Subject == "Patient", Tissue, Subject),
    x_axis_var = factor(
      x_axis_var,
      levels = c("UW volunteer 1", "UW volunteer 2", "UW volunteer 3", "UW volunteer 4", "Buffy coat", "Family member A", "Family member B", "UW volunteer 5", "Family member C", "UW volunteer 6", "UW volunteer 7")
    )
  )


################################################################################
############################ Grouped pathogenicity #############################
################################################################################

########## proportion of pathogenic/benign

pathogenicity <- skyscraper_prep2 %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  mutate(CTx = if_else(Subject %in% ctx_subjects, "CTx", "no-CTx")) %>%
  #filter(!(Subject == "Patient" | Subject == "UW volunteer 7")) %>%
  mutate(pathogenicity_am = case_when(
    color_group %in% c("likely_pathogenic", "likely_pathogenic_LC") ~ "pathogenic",
    color_group %in% c("likely_benign", "likely_benign_LC", "ambiguous") ~ "benign",
    TRUE ~ color_group)) %>%
  group_by(LFS, CTx, pathogenicity_am) %>%
  summarise(count = n()) %>%
  spread(pathogenicity_am, count, fill = 0)

pathogenicity

## chi-squared
# chi_test <- chisq.test(pathogenicity[,c("benign", "pathogenic")])
# chi_test
# p_value = chi_test$p.value

## fishers
fisher_test <- fisher.test(pathogenicity[, c("benign", "pathogenic")])
fisher_test
p_value <- fisher_test$p.value

mutation_proportions <- pathogenicity %>%
  mutate(total = pathogenic + benign,
         prop_pathogenic = pathogenic / total,
         prop_benign = benign / total)

mutation_proportions

mutation_proportions_long <- mutation_proportions %>%
  mutate(label = paste(LFS,CTx, sep="\n")) %>%
  pivot_longer(cols = c("prop_pathogenic", "prop_benign"), 
               names_to = "mutation_type", 
               values_to = "proportion") %>%
  mutate(mutation_type = recode(mutation_type, 
                                "prop_pathogenic" = "pathogenic", 
                                "prop_benign" = "benign/ambiguous"))
  
order_levels <- mutation_proportions_long %>%
  filter(mutation_type == "pathogenic") %>%
  arrange(proportion) %>% 
  pull(label)
mutation_proportions_long <- mutation_proportions_long %>%
  mutate(label = factor(label, levels = order_levels))

# Plot the stacked bar plot

pathogenicity_bar <- ggplot(mutation_proportions_long, aes(x = label, y = proportion, fill = mutation_type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +  
  labs(y = "Proportion of TP53 mutations", fill = "Mutation Type") +  
  theme_minimal() +
  scale_fill_manual(values = c("pathogenic" = "#E08143", "benign/ambiguous" = "#92A8D6")) +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=8), 
        axis.text.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size =8),
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(8, "pt"),
        legend.margin = margin(b=-8, t=0),
        plot.margin = margin(l=2, t=2, r=2, b=2),
        legend.key.spacing.y = unit(1, "pt"))  +
  labs(y = expression(italic("TP53")~"mutation proportion")) +
  guides(fill = guide_legend(nrow=2)) +
  annotate("text",
    x = 2.5, y = 1.1,
    label = deparse(bquote(italic(p) == .(round(p_value, 3)))),
    parse = TRUE,
    size = 8*25.4/72.27,
    color = "black"
  )
pathogenicity_bar
ggsave("results/pathogenicity_bar_fisher_ms3H.png", pathogenicity_bar, width = 2.5, height = 2, units = "in", dpi = 300)
#ggsave("/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_3/tp53_pathogenicity_ms3H.png", combined_plot_legend_crop, width = 2.5, height = 2, units = "in", dpi = 300)





########## pathogenic vs benign large clones

fancy_scientific <- function(l) { 
  cat("input breaks")
  print(l)
  #l <- l[!is.na(l)]
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  out <- parse(text = l)
  print(out)
  return(out)
}

pathogenicity_size <- skyscraper_prep2 %>%
  mutate(LFS = if_else(Subject %in% lfs_subjects, "LFS", "non-LFS")) %>%
  #filter(!(Subject == "Patient" | Subject == "UW volunteer 7")) %>%
  mutate(pathogenicity_am = case_when(
    color_group %in% c("likely_pathogenic", "likely_pathogenic_LC") ~ "pathogenic",
    color_group %in% c("likely_benign", "likely_benign_LC") ~ "benign",
    TRUE ~ color_group)) %>%
  filter(pathogenicity_am != "ambiguous") 

## mann-whitney
mann_whit <- wilcox.test(VAF ~ pathogenicity_am, data = pathogenicity_size)
pval <- signif(mann_whit$p.value, 3)

plot_vaf <- ggplot(pathogenicity_size, aes(x = pathogenicity_am, y = VAF, fill = pathogenicity_am)) +
  geom_jitter(aes(color = pathogenicity_am), width = 0.2, height = 0, size = 2, alpha = 0.95) +  # Set color for dots
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  labs(x = "Mutation Type", y = "VAF") +
  scale_y_log10(limits = c(4e-5, 8e-4), breaks = c(4e-5, 2e-4, 8e-4)) +
  #scale_y_log10(labels=fancy_scientific, breaks = c(1e-5, 3e-5, 1e-4, 3e-4, 1e-3)) +
  scale_color_manual(values = c("benign" = "#92A8D6", "pathogenic" = "#E08143")) +  # Color for dots
  scale_fill_manual(values = c("benign" = "#92A8D6", "pathogenic" = "#E08143")) +  # Color for fill in boxplots
  theme_minimal() +
  theme(#axis.text.x = element_text(size = 8, angle=45),
        axis.text.x = element_text(size = 8, angle=0),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=8),
        legend.position = "none",
        axis.title.x = element_blank(),
        #plot.margin = margin(b=-13, t=6,r=2,l=2),
        plot.margin = margin(10,1,1,1)) +
  annotate("text",
           x = 1,
           y = 7e-4,
           label = deparse(bquote(italic(p) == .(signif(pval, 3)))),
           parse= TRUE,
           size = 8*25.4/72.27,
           hjust = 0.5) 



ggsave("results/pathogenicity_sizes_ms.png", plot_vaf, width = 2, height = 1.5, units = "in", dpi = 300)
#ggsave("results/pathogenicity_sizes_ms.png", plot_vaf, width = 1.5, height = 1.5, units = "in", dpi = 300)
