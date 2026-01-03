#####
## Omega dn/ds plots
#####



#### bayes
data_bayes <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_bayes.tsv") %>% 
  filter(impact == "nonsynonymous")

bayes_dnds <- ggplot(data_bayes, aes(x = gene, y = mean_dnds)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = perc_25_dnds, ymax = perc_75_dnds), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  #scale_y_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(0, 3)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
bayes_dnds
ggsave("results/dnds_bayes.png", bayes_dnds, width = 5, height = 4, units = "in", dpi = 300)

data_counts_bayes <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_bayes.tsv")


### mle
data_mle <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_mle.tsv") %>%
  filter(impact == "nonsynonymous")

mle_dnds <- ggplot(data_mle, aes(x = gene, y = dnds)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene", y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
mle_dnds
ggsave("results/dnds_mle.png", mle_dnds, width = 5, height = 4, units = "in", dpi = 300)


################################################################################
########## mutation counts (same for bayes or mle)
################################################################################

bayes_counts_prop <- ggplot(data_counts_bayes, aes(x = gene, y = mutations, fill = impact)) +
  geom_col(position = "fill", color = "black") +
  scale_fill_manual(values = c("nonsynonymous" = "firebrick", "synonymous" = "steelblue")) +
  labs(x = "Gene", y = "Mutation proportion", fill = "Impact") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
bayes_counts_prop
ggsave("results/dnds_bayes_counts_prop.png", bayes_counts_prop, width = 5, height = 4, units = "in", dpi = 300)

bayes_counts <- ggplot(data_counts_bayes, aes(x = gene, y = mutations, fill = impact)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("nonsynonymous" = "firebrick", "synonymous" = "steelblue")) +
  labs(x = "Gene", y = "Number of mutations", fill = "Impact") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
bayes_counts
ggsave("results/dnds_bayes_counts.png", bayes_counts, width = 5, height = 4, units = "in", dpi = 300)


# data_counts <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/output_mut_table.tsv")
# df_summary <- data_counts %>%
#   mutate(total_mutations = rowSums(dplyr::select(., -GENE, -IMPACT, -CONTEXT_MUT))) %>%
#   group_by(GENE, IMPACT, CONTEXT_MUT) %>%
#   summarise(total_mutations = sum(total_mutations), .groups = "drop") 
# 
# df_summary <- df_summary %>% filter(GENE == "CBL")
# 
# data_counts <- ggplot(df_summary, aes(x = CONTEXT_MUT, y = total_mutations)) +
#   geom_col() +
#   facet_wrap(~ GENE, scales = "free_x") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(x = "Context",
#        y = "Total mutations")
# data_counts
# ggsave("results/context_counts.png", data_counts, width = 8, height = 4, units = "in", dpi = 300)

################################################################################
####### compare normal data to downsampled results
################################################################################

data_bayes_downsample_10000 <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_bayes_constant.tsv") %>% 
  filter(impact == "nonsynonymous")
data_bayes_downsample_1000 <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_bayes_constant.tsv") %>% 
  filter(impact == "nonsynonymous")

data_mle_downsample_10000 <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_mle_constant.tsv") %>%
  filter(impact == "nonsynonymous")
data_mle_downsample_1000 <- read_delim("/Volumes/feder-vol1/project/li_fraumeni/scripts/2025-07-18-omega_dnds/dn_ds_output_mle_constant.tsv") %>%
  filter(impact == "nonsynonymous")


data_bayes_downsample_10000 <- data_bayes_downsample_10000 %>% mutate(depth = "10000")
data_bayes_downsample_1000 <- data_bayes_downsample_1000 %>% mutate(depth = "1000")
data_bayes <- data_bayes %>% mutate(depth = "actual")
data_mle_downsample_10000 <- data_mle_downsample_10000 %>% mutate(depth = "10000")
data_mle_downsample_1000 <- data_mle_downsample_1000 %>% mutate(depth = "1000")
data_mle <- data_mle %>% mutate(depth = "actual")

downsample_bayes <- rbind(data_bayes_downsample_10000, data_bayes_downsample_1000, data_bayes)
downsample_mle <- rbind(data_mle_downsample_10000, data_mle_downsample_1000, data_mle)

downsample_mle
mle_dnds <- ggplot(downsample_mle, aes(x = depth, y = dnds, fill = depth)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 6)) +
  facet_wrap(~ gene) +
  labs(y = "dN/dS", fill = "Depth") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold")
  )
mle_dnds

bayes_dnds <- ggplot(downsample_bayes, aes(x = depth, y = dnds, fill = depth)) +
  geom_col() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 6)) +
  facet_wrap(~ gene) +
  labs(y = "dN/dS", fill = "Depth") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(face = "bold")
  )
bayes_dnds

bayes_dnds <- ggplot(downsample_bayes, aes(x = depth, y = mean_dnds, fill = depth)) +
  geom_col() +
  geom_errorbar(aes(ymin = perc_25_dnds, ymax = perc_75_dnds), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  scale_y_continuous(limits = c(0, 2)) +
  facet_wrap(~ gene) +
  labs(y = "dN/dS") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
bayes_dnds
