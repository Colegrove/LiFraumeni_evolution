library(ggplot2)

set.seed(2)

n <- 20
vaf_low  <- runif(round(n * 0.7),  0, 0.001)
vaf_high <- runif(round(n * 0.3),  0.001, 0.03)
vaf <- c(vaf_low, vaf_high)[1:n]

df <- data.frame(
  Mutation = paste0("M", seq_len(n)),
  VAF      = vaf
)

# User‐adjustable detection threshold (in % VAF)
threshold <- 0.00001 

# Create label strings
und_label <- paste0("Undetectable (<", threshold, ")")
det_label <- paste0("Detectable (≥", threshold, ")")


# 2) Classify by threshold
df$status <- ifelse(df$VAF >= threshold, det_label, und_label)
df$status <- factor(df$status, levels = c(und_label, det_label))

# 3) Define named colors
cols <- c("grey70", "#E69F00")
names(cols) <- c(und_label, det_label)

# 4) Plot
ggplot(df, aes(x = Mutation, y = VAF, color = status, size = VAF)) +
  # threshold line
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
  
  # points
  geom_point() +
  
  # apply named colors
  scale_color_manual(values = cols) +
  
  # size ~ VAF
  scale_size_continuous(range = c(2, 8), guide = FALSE) +
  
  # y-axis up to 3%
  scale_y_continuous(
    breaks = seq(0, 0.03, 0.01),
    labels = paste0(seq(0, 0.03, 0.01)),
    limits = c(0, 0.03)
  ) +
  
  labs(
    #x     = "Mutation",
    y     = "Variant Allele Frequency (VAF)",
    color = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5),
    #axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank()
  )





drop_trailing_zeros <- function(x) {
  s <- format(round(x, 4), nsmall = 4, scientific = FALSE, trim = TRUE)
  sub("\\.?0+$", "", s)
}

ggplot(df, aes(x = Mutation, y = VAF, color = status)) +
  geom_hline(yintercept = threshold, linetype = "dashed") +
  geom_point(aes(size = VAF)) +
  scale_color_manual(values = cols) +
  scale_size_continuous(range = c(2, 8), guide = FALSE) +
scale_y_log10(
  breaks       = log_breaks(n=4),
  #minor_breaks = NULL,
  #labels = drop_trailing_zeros
) +
  
  labs(
    y     = "Variant Allele Frequency (VAF)",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x        = element_blank(),
    panel.grid.major.x = element_blank()
  )





filt_maf %>% filter(Hugo_Symbol == "TP53") %>% filter(Subject == "Patient") %>% filter(coding == "coding") %>% arrange(VAF) %>% print(width = Inf)

