require(foreach)
require(tidyverse)
require(MASS)



### script to calculate and plot single regressions for MF and MB by gene
## First generate CHIP_muts_by_gene_non_coding.csv from xxxx

# mutation_frequencies_blood_CHIP_ms_2BC_v3.R


#coding = FALSE # true for coding false for non-coding

#coding_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene.csv"
non_coding_path <- "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/CHIP_muts_by_gene_non_coding.csv"


path = non_coding_path


# call in main 
#coding = FALSE

if(coding){
  coding_filt = "coding"
}else{
  coding_filt = "non-coding-CHIP"
}


chip_counts <- tibble(read.csv(path)) %>%
  filter(coding == coding_filt) %>%
  dplyr::select(Subject, denominator, age, Hugo_Symbol, 
                n_muts, mutReads, LFS, CTx) %>% 
  dplyr::rename(depth = denominator, 
         n_mutReads = mutReads, 
         gene = Hugo_Symbol) %>% 
  mutate(CTx = ifelse(is.element(Subject, c("Patient", "UW volunteer 7")), 1, 0), 
         LFS = ifelse(is.element(Subject, c("Patient", "Family member A",  "Family member C")), 1, 0))

## append the average number of mutations per sample for exclusion
ave_counts <- chip_counts %>%
  group_by(gene) %>%
  summarise(total_muts = sum(n_muts, na.rm = TRUE), .groups = "drop") %>%
  mutate(ave_muts_per_sample = total_muts / 11)
total_row <- tibble(
  gene = "Total",
  total_muts = sum(ave_counts$total_muts, na.rm = TRUE),
  ave_muts_per_sample = sum(ave_counts$total_muts, na.rm = TRUE) / 11
)
ave_counts <- bind_rows(ave_counts, total_row) %>% dplyr::select(-total_muts)
ave_counts %>% print(n = Inf)
geneList <- c(unique(chip_counts$gene), "Total")
totals <- chip_counts %>% group_by(Subject, age, LFS, CTx) %>% 
  summarize(depth = sum(depth), n_muts = sum(n_muts), 
            n_mutReads = sum(n_mutReads)) %>%
  mutate(gene = "Total") %>%
  ungroup()

depth_adjust = 10000000
chip_dat <- bind_rows(chip_counts, totals) %>% arrange(Subject, gene) %>% 
  mutate(decades = age / 10, depth = depth/depth_adjust) 

chip_dat <- chip_dat %>%
  mutate(MF = n_muts/depth,
         MB = n_mutReads/depth)

#with and without CTx status
models <- c(#"n_muts ~ 1 + decades + depth + LFS", 
            #"n_mutReads ~ 1 + decades + depth + LFS",
            #"n_muts ~ 1 + decades + depth + LFS + CTx",
            #"n_mutReads ~ 1 + decades + depth + LFS + CTx",
            #"MF ~ 1 + decades + depth + LFS",
            #"MF ~ 1 + decades + depth + LFS + CTx",
            "MF ~ 1 + decades + LFS + CTx",
            #"MB ~ 1 + decades + depth + LFS",
            #"MB ~ 1 + decades + depth + LFS + CTx",
            "MB ~ 1 + decades + LFS + CTx")


allRegressions <- 
  foreach(familyToCheck = c("gaussian", "poisson"), .combine = "rbind")%do%{
    foreach(model = models, .combine = "rbind")%do%{
      foreach(excludeCON07 = c("TRUE", "FALSE"), .combine = "rbind")%do%{
        foreach(geneToCheck = geneList, .combine = "rbind")%do%{
          
          if(excludeCON07){ 
            dat = chip_dat %>% filter(Subject != "UW volunteer 7")
          }else{
            dat = chip_dat
          }
          dat <- dat %>% filter(gene == geneToCheck)
          
          mod <- summary(glm(model, data = dat, family = familyToCheck))
          
          call <- paste(mod$call)
          coef <- mod$coefficients
          #This is just so we can standardize across Gaus and Pois links                   
          colnames(coef) <- c("Estimate", "Std. Error", "Stat", "p")
          regressionResults <- bind_cols(model_type = call[1], call = model, 
                                         model = call[3], 
                                         family = familyToCheck,
                                         noCON07 = excludeCON07,
                                         gene = geneToCheck,
                                         variables = rownames(coef), 
                                         coef) 
        }
      }
    }
  }

cols <- c( "#389688",  "#731433" )
names(cols) <- c( 0, 1)


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
chip_dat <- chip_dat %>%
  mutate(shape_group = case_when(LFS == 1 & CTx == 0 ~ "LFS/no-CTx",
                                 LFS == 1 & CTx == 1 ~ "LFS/CTx",
                                 LFS == 0 & CTx == 0 ~ "non-LFS/no-CTx",
                                 LFS == 0 & CTx == 1 ~ "non-LFS/CTx"))

chip_dat_exclude <- chip_dat %>% left_join(ave_counts) %>%
  filter(ave_muts_per_sample >= 2)


pvals_age <- chip_dat_exclude %>%
  group_by(gene) %>%
  do({
    m <- lm(MF ~ decades, data = .)
    tidy(m)
  }) %>%
  ungroup() %>%
  filter(term == "decades") %>%
  mutate(
    gene  = factor(gene, levels = geneList),
    label = paste0("italic(p) == ", signif(p.value, 2))
  )


freq_depth <- chip_dat_exclude %>% 
  ggplot(aes(x = depth, y = MF)) +
  #ggplot(aes(x = depth, y = n_muts)) + 
  geom_point(aes(color = shape_group, shape=shape_group)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") + 
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  facet_wrap(~factor(gene, levels = geneList), scales = "free") +
  geom_smooth(method='lm') +
  #labs(x = "Reads/base (in 10000000s)", y = "Number of unique mutations", col = "LFS") +
  labs(x = "Reads/base (in 10000000s)", y = "MF", col = "LFS") +
  theme_classic() + 
  #scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
freq_depth
freq_age <- chip_dat_exclude %>% 
  ggplot(aes(x = decades, y = MF)) +
  #ggplot(aes(x = decades, y = n_muts)) + 
  geom_point(aes(color = shape_group, shape=shape_group)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") + 
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  facet_wrap(~factor(gene, levels = geneList), scales = "free") +
  geom_smooth(method='lm') +
  #labs(x = "Age (in decades)", y = "Number of unique mutations", col = "LFS") +
  labs(x = "Age (in decades)", y = expression("MF (mutations / " * 10^7 * " bases)"), col = "LFS") +
  theme_classic() + 
  geom_text(
    data = pvals_age,
    aes(x = Inf, y = Inf, label = label),
    parse = TRUE,
    hjust = 1.1, vjust = 1.5,
    size  = 2.8,          # ~8 pt
    inherit.aes = FALSE
  ) +
  #scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

burden_depth <- chip_dat_exclude %>% 
  ggplot(aes(x = depth, y = MB)) +
  #ggplot(aes(x = depth, y = n_mutReads)) + 
  geom_point(aes(color = shape_group, shape=shape_group)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") + 
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  facet_wrap(~factor(gene, levels = geneList), scales = "free") +
  geom_smooth(method='lm') +
  #labs(x = "Reads/base (in 10000000s)", y = "Number of mutant reads", col = "LFS") +
  labs(x = "Reads/base (in 10000000s)", y = "MB", col = "LFS") +
  theme_classic() + 
  #scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

burden_age <- chip_dat_exclude %>% 
  ggplot(aes(x = decades, y = MB)) +
  #ggplot(aes(x = decades, y = n_mutReads)) + geom_point() +
  geom_point(aes(color = shape_group, shape=shape_group)) +
  scale_color_manual(values = shape_group_colors, name = "Patient history") + 
  scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
  facet_wrap(~factor(gene, levels = geneList), scales = "free") +
  geom_smooth(method='lm') +
  #labs(x = "Age (in decades)", y = "Number of mutant reads", col = "LFS") +
  labs(x = "Age (in decades)", y = "MB", col = "LFS") +
  theme_classic() + 
  geom_text(
    data = pvals_age,
    aes(x = Inf, y = Inf, label = label),
    parse = TRUE,
    hjust = 1.1, vjust = 1.5,
    size  = 2.8,
    inherit.aes = FALSE
  ) +
  #scale_color_manual(values = cols, labels = c("0" = "non-LFS", "1" = "LFS")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


### plot for supplement
supp_theme <- function() {
  theme(
    axis.title.x = element_text(size = 8),
    axis.text.x  = element_text(size = 8, angle = 0, vjust = 0, hjust = 0.5),
    axis.title.y = element_text(size = 8),
    axis.text.y  = element_text(size = 8),
    strip.text   = element_text(size = 8),
    legend.text  = element_text(size = 8),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box.margin = margin(-8, 0, -8, 0),
    legend.margin = margin(0, 0, 0, 0)
  )
}

freq_age_supp <- freq_age + supp_theme() +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(size = 12),
    plot.tag.position = c(0, 1)
  )
freq_depth_supp <- freq_depth + supp_theme()
burden_depth_supp <- burden_depth + supp_theme()
burden_age_supp <- burden_age + supp_theme()


height <- 6
width <- 6
if(coding){
  #ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_depth.png', freq_depth_supp, height = height, width = width)
  ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age.png', freq_age_supp, height = height, width = width)
  ggsave(filename = "/Users/huntc10/Library/CloudStorage/OneDrive-UW/Li-Fraumeni/Manuscript_figures/Fig_S2/MF_by_gene_supp_S2.png", freq_age_supp, height = height, width = width)
  
  #ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_depth.png', burden_depth_supp, height = height, width = width)
  ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutBurden_by_age.png', burden_age_supp, height = height, width = width)
} else {
  ggsave(filename = '/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/gene_mutFreq_by_age_non_coding.png', freq_age_supp, height = height, width = width)
}


## Ok, let's start by just looking at Gaussian, all individuals included no CTx
base_frequency <- allRegressions %>% filter(family == "gaussian", 
                                            #call == "n_muts ~ 1 + decades + depth + LFS", 
                                            call == "MF ~ 1 + decades + depth + LFS",
                                            noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest') %>%
  left_join(ave_counts)

base_burden <- allRegressions %>% filter(family == "gaussian", 
                                         #call == "n_mutReads ~ 1 + decades + depth + LFS", 
                                         call == "MB ~ 1 + decades + depth + LFS",
                                         noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest') %>%
  left_join(ave_counts)

## Include CTx as a covariate
base_frequency_CTx <- allRegressions %>% filter(family == "gaussian", 
                                            #call == "n_muts ~ 1 + decades + depth + LFS + CTx", 
                                            #call == "MF ~ 1 + decades + depth + LFS + CTx",
                                            call == "MF ~ 1 + decades + LFS + CTx",
                                            noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest') %>%
  left_join(ave_counts)

base_burden_CTx <- allRegressions %>% filter(family == "gaussian", 
                                         #call == "n_mutReads ~ 1 + decades + depth + LFS + CTx", 
                                         #call == "MB ~ 1 + decades + depth + LFS + CTx",
                                         call == "MB ~ 1 + decades + LFS + CTx",
                                         noCON07 == FALSE) %>% 
  mutate( Estimate = paste0(signif(Estimate, 4), 
                            " (", signif(`Std. Error`, 4) ,")"), 
          p = signif(p, 3)) %>%
  dplyr::select(-model_type, -model, -family, -noCON07, -`Std. Error`, -Stat)  %>% 
  pivot_wider(names_from = variables,
              values_from = c(Estimate, p),
              names_glue = "{variables}.{.value}", 
              names_vary = 'slowest') %>%
  left_join(ave_counts)


if(coding){
  #write.csv( base_frequency, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency.csv", row.names = FALSE )
  #write.csv( base_burden, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden.csv", row.names = FALSE )
  write.csv( base_frequency_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_CTx.csv", row.names = FALSE )
  write.csv( base_burden_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_CTx.csv", row.names = FALSE )
  #write.csv( base_frequency_BH, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_BH.csv", row.names = FALSE )
  #write.csv( base_burden_BH, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_BH.csv", row.names = FALSE )
  #write.csv( base_frequency_BH_filt, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_BH_filt.csv", row.names = FALSE )
  #write.csv( base_burden_BH_filt, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_BH_filt.csv", row.names = FALSE )
}else{
  #write.csv( base_frequency, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_non_coding.csv", row.names = FALSE )
  #write.csv( base_burden, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_non_coding.csv", row.names = FALSE )
  write.csv( base_frequency_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_frequency_non_coding_CTx.csv", row.names = FALSE )
  #write.csv( base_burden_CTx, "/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-10-20-multiple_regression/base_burden_non_coding_CTx.csv", row.names = FALSE )
}
