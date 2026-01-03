library(tidyverse)
library(sigfit)
library(Cairo)

#install.packages("pkgbuild")
#devtools::install_github("kgori/sigfit", build_opts = c("--no-resave-data", "--no-manual"))
data("cosmic_signatures_v3.2")

#######
### Patient mutagenesis
#######

## import matrix from sigProfiler matrix generator
sigPath <- "/Volumes/feder-vol1/project/li_fraumeni/dat/mutationalSignatures/output/SBS/li_fraumeni.SBS96.region"
df <- read_delim(sigPath, delim = "\t")
SBS_96_counts_all <- df %>%
  pivot_longer(cols = -MutationType, names_to = "sample", values_to = "Count")

SBS_96_counts_all %>% filter(sample == "Skin") %>%print(n=Inf)

## convert matrix for use in sigFit
mut_mat_wide <- pivot_wider(SBS_96_counts_all, names_from = "MutationType", values_from = "Count", values_fill = list(Count = 0))
mutation_catalog <- as.matrix(mut_mat_wide[,-1])  # Exclude sample column
rownames(mutation_catalog) <- mut_mat_wide$sample

## convert to other format
# Define function to convert mutation types
convert_trinucleotide <- function(trinucleotide) {
  # Extract the base before, the mutation, and the base after
  match <- regmatches(trinucleotide, regexec("([ACGT])\\[([ACGT])>([ACGT])\\]([ACGT])", trinucleotide))
  print(match)

  ref_first <- match[[1]][2]  # First base (e.g., A)
  ref_mid <- match[[1]][3]    # Middle base (before mutation) (e.g., C)
  alt_mid <- match[[1]][4]    # Middle base (after mutation) (e.g., A)
  ref_last <- match[[1]][5]   # Last base (e.g., A)
  
  cat("ref_first")
  print(ref_first)
  cat("ref_mid")
  print(ref_mid)
  # Construct the before and after mutation sequences
  before_mut <- paste0(ref_first, ref_mid, ref_last)
  after_mut <- paste0(ref_first, alt_mid, ref_last)
    
  out <- paste0(before_mut, ">", after_mut)
  cat("out")
  print(out)
    
}

#colnames(mut_mat_wide) <- sapply(colnames(mut_mat_wide), convert_trinucleotide)

colnames(mutation_catalog) <- sapply(colnames(mutation_catalog), convert_trinucleotide)
sigfit_column_order <- colnames(cosmic_signatures_v3.2)
current_colnames <- colnames(mutation_catalog)
order_index <- match(sigfit_column_order, current_colnames)

mutation_catalog <- mutation_catalog[, order_index]

## fit to COSMIC signatures
#fit <- sigfit::fit_signatures(counts = mut_mat_wide, signatures = cosmic_signatures_v3.2)
fit <- sigfit::fit_signatures(counts = mutation_catalog, signatures = cosmic_signatures_v3.2)
plot_spectrum(mutation_catalog, "./results/spectrum.pdf")
plot_exposures(fit, pdf_path = "./results/exposures.pdf")


#plot_spectrum(mut_mat_wide, "./results/spectrum.pdf")
#plot_exposures(fit, pdf_path = "./results/exposures.pdf")




