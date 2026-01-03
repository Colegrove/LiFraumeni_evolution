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
  filter(Subject %in% family_patient_blood | Subject %in% mstp) %>%
  filter(Tissue %in% family_patient_blood_samples) %>%
  mutate(age = age_map[Subject]) %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(age = age_map[Subject]) %>%
  mutate(plot_coding = if_else(!is.na(exon_number), "coding", "non-coding"))
mutFreq_prep %>% print(width = Inf)


lfs_subjects = c("Patient", "Family member A", "Family member C")
ctx_subjects = c("Patient", "UW volunteer 7")
### filter for plot
mutFreq_prep <- mutFreq_prep %>%
  left_join(tp53_depth_split) %>%
  group_by(plot_coding, Subject, denominator_coding, denominator_noncoding, age) %>%
  summarise(n_muts = n()) %>%
  mutate(mutFreq = if_else(plot_coding == "coding", n_muts/denominator_coding, n_muts/denominator_noncoding)) %>%
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

# coding_shapes <- c("coding" = 17, "non-coding" = 1)  # triangle vs circle
# 
# # Colors for patient history
# shape_group_colors <- c(
#   "non-LFS/no-CTx" = "#44AA99",
#   "non-LFS/CTx"    = "#44AA99",
#   "LFS/no-CTx"     = "#882255",
#   "LFS/CTx"        = "#882255"
# )
# 
# 
# # Shapes for region: circle vs triangle (both fillable)
# coding_shapes <- c("coding" = 24, "non-coding" = 21)  # triangle vs circle
# 
# LFS_colors <- c("LFS" = "#882255", "non-LFS" = "#44AA99")
# 
# ## version 3 no vertical lines and coding/non-coding regression
# mutFreq_plot3 <- ggplot(mutFreq_prep, aes(x = age_j, y = mutFreq)) +
#   geom_point(aes(color = LFS, fill = LFS, shape = plot_coding),
#              size = GEOM_POINT_SIZE + 1, 
#              stroke = 0.8) +
#   geom_smooth(
#     aes(x=age, y=mutFreq, color = LFS, linetype = plot_coding, group = interaction(LFS, plot_coding)),
#     method = "lm", se = FALSE, linewidth = 0.7
#   ) +
#   scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
#   #scale_y_log10(labels = fancy_scientific, limits = c(1e-7, 3.1e-6)) +
#   scale_y_log10() +
#   scale_color_manual(values = LFS_colors, name = "Patient history") +
#   scale_fill_manual(values = LFS_colors, name = "Patient history") +
#   scale_shape_manual(values = coding_shapes, name = "Region") +
#   labs(x = "Age", y = "Mutation frequency") +
#   theme_minimal() + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   theme(legend.position = "top",
#         legend.title = element_blank(),
#         legend.text  = element_text(size = 8, margin=margin(r=3, l=2)),
#         legend.margin = margin(r=1, l=1, b = 1, t = 1),
#         legend.box.margin = margin(r=0, l=-30, b = 0, t = 0),
#         legend.key.size   = unit(8, "pt"),
#         legend.spacing.x  = unit(0.2, "cm"),
#         legend.spacing.y  = unit(0.2, "cm"),
#         #legend.background = element_rect(color = "black", fill = "white", linewidth = 0.3), 
#         axis.text.x = element_text(size=8),
#         axis.text.y = element_text(size=8),
#         axis.title.x = element_text(size=8),
#         axis.title.y = element_text(size=8)
#   ) + 
#   guides(linetype = guide_legend(override.aes = list(color = 'black'), ncol=1),
#          shape = guide_legend(ncol = 1),
#          color = guide_legend(ncol =1))
# 
# 
# mutFreq_plot3
# ggsave("results/tp53_LFS_status_mutation_frequency_non_coding_ms.png", mutFreq_plot3, width = 3, height = 2.5, units = "in", dpi = 300)

################################################################################
############################## coding -> non-coding lines
################################################################################
#MF_comparison_prep <- mutFreq_prep %>% dplyr::select(-age_offset, -age_j)
#write.csv(MF_comparison_prep, "/Volumes/feder-vol1/project/li_fraumeni/dat/TP53_MF_comparison.csv")


# MF_coding_comparison <- ggplot(MF_comparison_prep, 
#                               aes(x = factor(plot_coding, levels = c("non-coding", "coding")),
#                                   y = mutFreq)) +
#   geom_line(aes(group = Subject, color = LFS),
#             linewidth = 0.7, alpha = 0.7) +
#   geom_point(aes(color = LFS, fill = LFS, shape = plot_coding),
#              size = GEOM_POINT_SIZE + 1, stroke = 0.8) +
#   scale_y_log10(limits = c(3e-8, 1.1e-6)) +
#   scale_color_manual(values = LFS_colors, name = "Patient history") +
#   scale_fill_manual(values = LFS_colors, name = "Patient history") +
#   scale_shape_manual(values = coding_shapes) +
#   labs(y = "TP53\nmutation frequency") +
#   theme_minimal() +
#   theme(
#     legend.position = "top",
#     legend.title = element_blank(),
#     legend.text  = element_text(size = 8, margin = margin(r=3, l=2)),
#     legend.margin = margin(r=1, l=1, b=1, t=1),
#     legend.box.margin = margin(r=0, l=-30, b=0, t=0),
#     legend.key.size   = unit(8, "pt"),
#     legend.spacing.x  = unit(0.2, "cm"),
#     legend.spacing.y  = unit(0.2, "cm"),
#     axis.text.x = element_text(size=8),
#     axis.text.y = element_text(size=8),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size=8)
#   ) +
#   guides(
#     shape = guide_legend(ncol = 1),
#     color = guide_legend(ncol = 1)
#   )
# 
# MF_coding_comparison
# ggsave("results/MF_coding_comparison_tp53_ms.png", MF_coding_comparison, width = 2.5, height = 2, units = "in", dpi = 300)

################################################################################
############################## 
################################################################################

# shape_group_colors <- c(
#   "non-LFS/no-CTx" = "#44AA99",  # non-LFS color
#   "non-LFS/CTx"     = "#44AA99",  # same as non-LFS
#   "LFS/no-CTx"     = "#882255",  # LFS color
#   "LFS/CTx" = "#882255"   # same as LFS
# )
# 
# # Define shapes by patient history
# shape_group_shapes <- c(
#   "non-LFS/no-CTx" = 1,
#   "non-LFS/CTx"     = 16,
#   "LFS/no-CTx"     = 1,
#   "LFS/CTx" = 16
# )

ctx_shapes <- c(
  "non-CTx" = 1,   # open circle
  "CTx"    = 16   # filled circle
)
LFS_colors <- c(
  "non-LFS" = "#44AA99",  # non-LFS color
  "LFS"     = "#882255"  # LFS color
)

MF_comparison_prep <- mutFreq_prep %>% dplyr::select(-age_offset, -age_j)
MF_ratio_prep <- MF_comparison_prep %>%
  dplyr::select(Subject, age, LFS, CTx, shape_group, mutFreq, plot_coding) %>%
  pivot_wider(
    names_from  = plot_coding,
    values_from = mutFreq
  ) %>%
  mutate(
    MF_ratio = coding / `non-coding`
  ) 
MF_ratio_prep
MF_coding_ratio <- ggplot(MF_ratio_prep,
                          aes(x = LFS, y = MF_ratio, color = LFS, shape = CTx)) +
  geom_jitter(width = 0.15, height = 0, size = GEOM_POINT_SIZE + 0.5, stroke = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_y_log10(limits = c(3e-1, 30)) +
  scale_color_manual(values = LFS_colors, name = "Patient history") +
  scale_shape_manual(values = ctx_shapes) +
  labs(
    x = NULL,
    y = "TP53 coding MF / non-coding MF"
  ) +
  theme_minimal() +
  theme(
    legend.position   = "top",
    legend.title      = element_blank(),
    legend.text       = element_text(size = 8, margin = margin(r=3, l=2)),
    legend.margin     = margin(r=1, l=1, b=1, t=1),
    legend.box.margin = margin(r=0, l=-30, b=0, t=0),
    legend.key.size   = unit(8, "pt"),
    legend.spacing.x  = unit(0.2, "cm"),
    legend.spacing.y  = unit(0.2, "cm"),
    axis.text.x       = element_text(size = 8),
    axis.text.y       = element_text(size = 8),
    axis.title.y      = element_text(size = 8)
  )

MF_coding_ratio
ggsave("results/MF_coding_ratio_tp53_ms.png", MF_coding_ratio, width = 2.5, height = 2, units = "in", dpi = 300)

################################################################################
############################## 
################################################################################

# 
# 
# mutFreq_prep %>% filter(plot_coding == "non-coding") %>% print(n = Inf)
# mutBurden_plot2 <- ggplot(mutFreq_prep %>% filter(plot_coding == "non-coding"), aes(x = age_j, y = mutFreq)) +
#   geom_smooth(data = mutFreq_prep %>% filter(plot_coding == "non-coding"),
#               se = FALSE, method = "lm", color = '#444444') +
#   # geom_smooth(data = mutFreq_prep %>% filter(LFS == "non-LFS"),
#   #             se = FALSE, method = "lm", color = nonLFS_color) +  # match non-LFS color
#   geom_point(aes(color = shape_group, shape = shape_group), size = GEOM_POINT_SIZE, alpha = 1, stroke = 1.2) +
#   scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, 5)) +
#   #scale_y_log10(limits = c(1.5e-6, 3e-5), labels=fancy_scientific) +
#   scale_y_log10() +
#   scale_color_manual(values = shape_group_colors, name = "Patient history") +
#   scale_shape_manual(values = shape_group_shapes, name = "Patient history") +
#   labs(
#     x = "Age",
#     y = "TP53 non-coding\nmutation frequency"
#   ) +
#   theme_minimal()
# 
# show(mutBurden_plot2)
# mutBurden_plot2 <- mutBurden_plot2 + theme(legend.position = "none")
# ggsave("results/tp53_non_coding_regression_ms.png", mutBurden_plot2, width = 1.5, height = 1.5, units = "in", dpi = 300)
# 
# 
# 
# 
# mutBurden_plot2 <- mutBurden_plot2 + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   ylab(expression(atop(NA, atop(textstyle("TP53 non-coding"), textstyle("mutation frequency"))))) +
#   theme(legend.position = "none",
#         text = element_text(size = 8),
#         axis.title = element_text(size = 8),
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8, margin = margin(l=-10), hjust = 0.65),
#         axis.text  = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.text  = element_text(size = 8))
# 
# mutFreq_plot2 <- mutFreq_plot2 + 
#   scale_x_continuous(breaks = c(25, 35, 45, 55, 65, 75)) +
#   guides(color = guide_legend(ncol = 1, title.position = "top"), shape = guide_legend(ncol = 1, title.position = "top")) +
#   ylab(expression(atop(NA, atop(textstyle("overall"), textstyle("mutation frequency"))))) +
#   theme(legend.position  = "none",
#         axis.title.x = element_text(size = 8),
#         axis.title.y = element_text(size = 8, margin= margin(l=-10), hjust = 0.65),
#         axis.text  = element_text(size = 8),
#         legend.title = element_text(size = 8),
#         legend.text  = element_text(size = 8))
# 
# plots <- plot_grid(
#   plot_grid(mutFreq_plot2, mutBurden_plot2, ncol = 2, align = "v"), 
#   ncol = 1)
# 
# combined_plots <- plot_grid(plots, legend_shared, ncol = 1, 
#                             rel_heights = c(1.1,0.225))
# 
# show(combined_plots)
# ggsave("results/tp53_non_coding_mutation_frequency.png", combined_plots, width = 3.5, height = 1.5, units = "in", dpi = 300)
# 
# 
# 
# 
# 
# 
# 
# 
# ############## multiple regression
# mutFreq_prep
# mutFreq_prep_CHIP_encode <- mutFreq_prep %>%
#   mutate(
#     LFS_n = if_else(LFS == "LFS", 1, 0),
#     CTx_n = if_else(CTx == "CTx", 1, 0)
#   )
# 
# mutFreq_prep_CHIP_encode
# 
# mutFreq_subject <- mutFreq_prep_CHIP_encode %>% 
#   filter(plot_coding == "non-coding") %>%
#   mutate(age_decades = age / 10,
#          depth_scaled = denominator_noncoding / 10000000) %>%
#   mutate(MF = n_muts/denominator_noncoding,
#          MF_scale = n_muts/depth_scaled)
#          #MB_scale = mutReads/depth_scaled)
# 
# mutFreq_subject
# # With UW07 and no CTx
# # model_freq_CHIP_all_samples_noCTx <- lm(
# #   #n_muts ~ age_decades + depth_scaled + LFS_n,
# #   #MF ~ age_decades + depth_scaled + LFS_n, 
# #   MF_scale ~ age_decades + depth_scaled + LFS_n,
# #   data = mutFreq_subject
# # )
# # summary(model_freq_CHIP_all_samples_noCTx)
# 
# ## with ctx
# model_freq_CHIP_all_samples_CTx <- lm(
#   #n_muts ~ age_decades + depth_scaled + LFS_n,
#   #MF ~ age_decades + depth_scaled + LFS_n, 
#   #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
#   MF_scale ~ age_decades + LFS_n + CTx_n,
#   data = mutFreq_subject
# )
# summary(model_freq_CHIP_all_samples_CTx)
# 
# 
# coef_freq_CHIP <- tidy(model_freq_CHIP_all_samples_CTx, conf.int = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(sig = case_when(
#     p.value < 0.0001 ~ "***",
#     p.value < 0.01   ~ "**",
#     p.value < 0.05   ~ "*",
#     TRUE ~ ""
#   ))
# 
# term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
# coef_freq_CHIP
# 
# coef_freq_CHIP <- tibble::tribble(
#   ~term,       ~estimate, ~std.error, ~p.value,
#   "age_decades",   0.3025,   0.06588,   0.00251,
#   "LFS_n",         0.8457,   0.2613,    0.0143,
#   "CTx_n",        -0.5406,   0.3139,    0.129
# ) %>%
#   mutate(
#     conf.low = estimate - 1.96 * std.error,
#     conf.high = estimate + 1.96 * std.error,
#     sig = case_when(
#       p.value < 0.001 ~ "***",
#       p.value < 0.01  ~ "**",
#       p.value < 0.05  ~ "*",
#       TRUE ~ ""
#     )
#   )
# 
# coef_freq_CHIP
# 
# CHIP_freq_plot <- ggplot(coef_freq_CHIP,
#                          aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
#   geom_point(color = "black", size = 1) +
#   geom_errorbarh(height = 0, color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_y_discrete(labels = term_labels) +
#   scale_x_continuous(limits = c(-1.2, 1.5), breaks = seq(-1, 1.5, 0.5)) +
#   geom_text(aes(label = sig, x = estimate),
#             hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
#   labs(x = "Effect size", y = NULL) +
#   theme_classic(base_size = 8) +
#   theme(axis.text.y = element_markdown(size = 8),
#         axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))
# 
# CHIP_freq_plot
# 
# 
# ####### freq + burden plot
# 
# CHIP_freq_plot <- CHIP_freq_plot + theme(
#   plot.margin = margin(1,1,1,1)
# )
# 
# CHIP_burden_plot <- CHIP_burden_plot + theme(
#   plot.margin = margin(1,1,1,1)
# )
# 
# 
# combined_plot <- CHIP_freq_plot | CHIP_freq_plot
# 
# 
# 
# ggsave("results/covariates_combined_CHIP_ms.png", combined_plot, width = 3.5, height = 1.5, units = "in", dpi = 300)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### tp53 only
# mutFreq_prep
# mutFreq_prep_CHIP_encode <- mutFreq_prep %>%
#   mutate(
#     LFS_n = if_else(LFS == "LFS", 1, 0),
#     CTx_n = if_else(CTx == "CTx", 1, 0)
#   )
# 
# mutFreq_prep_CHIP_encode
# 
# mutFreq_subject <- mutFreq_prep_CHIP_encode %>% 
#   filter(plot_coding == "non-coding") %>%
#   mutate(age_decades = age / 10,
#          depth_scaled = denominator_noncoding / 10000000) %>%
#   mutate(MF = n_muts/denominator_noncoding,
#          MF_scale = n_muts/depth_scaled)
# #MB_scale = mutReads/depth_scaled)
# 
# mutFreq_subject
# # With UW07 and no CTx
# # model_freq_CHIP_all_samples_noCTx <- lm(
# #   #n_muts ~ age_decades + depth_scaled + LFS_n,
# #   #MF ~ age_decades + depth_scaled + LFS_n, 
# #   MF_scale ~ age_decades + depth_scaled + LFS_n,
# #   data = mutFreq_subject
# # )
# # summary(model_freq_CHIP_all_samples_noCTx)
# 
# ## with ctx
# model_freq_CHIP_all_samples_CTx <- lm(
#   #n_muts ~ age_decades + depth_scaled + LFS_n,
#   #MF ~ age_decades + depth_scaled + LFS_n, 
#   #MF_scale ~ age_decades + depth_scaled + LFS_n + CTx_n,
#   MF_scale ~ age_decades + LFS_n + CTx_n,
#   data = mutFreq_subject
# )
# summary(model_freq_CHIP_all_samples_CTx)
# 
# 
# coef_freq_CHIP <- tidy(model_freq_CHIP_all_samples_CTx, conf.int = TRUE) %>%
#   filter(term != "(Intercept)") %>%
#   mutate(sig = case_when(
#     p.value < 0.0001 ~ "***",
#     p.value < 0.01   ~ "**",
#     p.value < 0.05   ~ "*",
#     TRUE ~ ""
#   ))
# 
# term_labels <- c("age_decades" = "Age", "depth_scaled" = "Depth", "LFS_n" = "LFS", "CTx_n" = "CTx")
# coef_freq_CHIP
# 
# coef_freq_CHIP <- tibble::tribble(
#   ~term,       ~estimate, ~std.error, ~p.value,
#   "age_decades",   0.2339,    0.1191,     0.0903,
#   "LFS_n",         1.495,     0.4724,     0.0158,
#   "CTx_n",        -1.052,     0.5676,     0.106
# ) %>%
#   mutate(
#     conf.low = estimate - 1.96 * std.error,
#     conf.high = estimate + 1.96 * std.error,
#     sig = case_when(
#       p.value < 0.001 ~ "***",
#       p.value < 0.01  ~ "**",
#       p.value < 0.05  ~ "*",
#       TRUE ~ ""
#     )
#   )
# 
# coef_freq_CHIP
# 
# 
# CHIP_freq_plot <- ggplot(coef_freq_CHIP,
#                          aes(y = term, x = estimate, xmin = conf.low, xmax = conf.high, color = sig)) +
#   geom_point(color = "black", size = 1) +
#   geom_errorbarh(height = 0, color = "black") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_y_discrete(labels = term_labels) +
#   scale_x_continuous(limits = c(-2.2, 2.5), breaks = seq(-3, 3, 1)) +
#   geom_text(aes(label = sig, x = estimate),
#             hjust = 0.5, vjust = 0.1, size = 5, color = "black") +
#   labs(x = "Effect size", y = NULL) +
#   theme_classic(base_size = 8) +
#   theme(axis.text.y = element_markdown(size = 8),
#         axis.text.x = element_text(size = 8, angle=45, vjust = 0.5))
# 
# CHIP_freq_plot
# 
# 
# ####### freq + burden plot
# 
# CHIP_freq_plot <- CHIP_freq_plot + theme(
#   plot.margin = margin(1,1,1,1)
# )
# 
# 
# 
# combined_plot <- CHIP_freq_plot | CHIP_freq_plot
# 
# 
# 
# ggsave("results/covariates_combined_CHIP_ms_tp53.png", combined_plot, width = 3.5, height = 1.5, units = "in", dpi = 300)
# 
# 
# 
# 
# 
# 
# 
