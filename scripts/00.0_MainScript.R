## Li-Fraumeni evolution main script

## Load Modules ####
library("tidyverse")
library("cowplot")
library("ggtext")
library("tidyr")
library("rtracklayer")
library("dplyr")
library("ggplot2")
library("gridtext")
library("tibble")
library("foreach")
library("MASS")
library("broom")
library("patchwork")
library("stringr")
library("vroom")
library("Rsamtools")
library("Biostrings")



# Load inputs list
inputs <- read_delim("processing_config.txt", delim="=", col_names = FALSE, comment = "#", skip_empty_rows = TRUE) %>%
  pivot_wider(names_from = X1, values_from = X2) %>% 
  type_convert()

## sample information
sample_data <- read_csv(inputs$extra_data_file)

## Load Functions####
source("scripts/01.1.04_loadMaf.R")
source("scripts/01.1.07_loadBed.R")
source("scripts/skyscraper_color.R")

# variant classification table renames subvariants to "intron", "exon", etc.
variant_clasification_table = read_delim(inputs$varClassTrans) 

## Load AlphaMissense file ####
alphamissense <- read_delim(inputs$AlphaMissense, delim = "\t", skip = 3) %>% 
  dplyr::rename('CHROM' = "#CHROM")

## Load IntOGen file
intOGen <- read_delim(inputs$intOGen, delim = "\t") %>% 
  dplyr::select(chr,pos,alt, boostDM_score, boostDM_class) %>% 
  mutate(chr = paste0("chr",chr))

# Read in BED for MUT
targ_bed12_MUT <- loadBed(inputs$bed12_MUT)
targ_bed6_MUT <- loadBed(inputs$bed6_MUT)

# Convert genome.mut to depth file
source("scripts/02.1.0_Mut2Depth.R")
# Read in Depth
source("scripts/02.1.1_LoadDepth.R")
## filter out samples that do not contain depth requirements 
source("scripts/02.1.2_FilterDepth.R")

# Read in MAF
source("scripts/02.4.1_LoadMaf.R")
source("scripts/02.4.2_FiltMaf.R")

# assign coding regions and repeat masking filter
source("scripts/apply_coding_depth.R")

# filter variants that existed at high VAF in other samples
source("scripts/snp_cross_contamination.R")

######## Depth analysis Fig 1

## S1 and Table S4
source("scripts/depth_by_panel_S1.R")
# Fig 1
source("scripts/depth_visualization_ms_1D.R")

######## Blood analysis
## MF CHIP regression: 2BD
source("scripts/MF_blood_coding_noncoding_ms_2BC_v4.R")
## multiple regression: 2CE
source("scripts/ms_2DE_coefficient.R")
# Fig S3
source("scripts/MB_blood_coding_supp_v4.R")
## TP53 regression: 3EF
source("scripts/MF_tp53_coding_noncoding_ms_3BC_v4.R")
## TP53 multiple regression: 3GH
source("scripts/ms_3DE_coefficient_v2.R")

## plot single regression genes coding and non-coding
## Figure S2, Table SX
coding = FALSE
source("scripts/regression_genes_plot_supp_v2.R")
coding = TRUE
source("scripts/regression_genes_plot_supp_v2.R")

## Table S3
source("scripts/supp_table_sequencing_overlap.R")

# Figure 3A
source("scripts/skyscraper_blood_ms_3A.R")
# Figure 2A
source("scripts/mutation_count_CHIP_ms_2A.R")

## TP53 MF coding/non-coding

######## Tissue analysis
### 4A and 4H
source("scripts/skyscraper_tissues_ms_4A.R")
source("scripts/variant_types_tissues_ms_4A.R")
source("scripts/mutation_tissue_overlap_ms_4A.R")
### S6
source("scripts/blood_contamination_tissues_supp.R")

## Figure 4C and Table S10
source("scripts/coding_non_coding_tissues_MF_ms_4D.R")

## Figure S7
source("scripts/mutation_tissue_overlap_all_supp.R")

## 3C/D & 4F/G
source("scripts/TP53_binding_domain_ms.R")

## 3K & S4
source("scripts/apply_all_alphamissense.R")
source("scripts/Alphamissense_grouped_ms_3_CHIP_supp.R")

## Fig 4D
source("scripts/mutsigs_tissues_ms.R")
## Fig 2F
source("scripts/mutsigs_blood_ms.R")
## Fig 2G
source("scripts/SBSG_spectra_plot.R")
## Fig 2H
source("scripts/all_possible_target_mutations.R")
source("scripts/all_possible_muts_annotations.R")
source("scripts/dnds_blood_ms_2E.R")
## Fig 3B
source("scripts/lollipop_blood_ms_3_v2.R")
## Fig 4E
source("scripts/lollipop_tissue_ms_4_v2.R")
## Fig 3I
source("scripts/MF_ratio_tp53_ms_3F.R")
## Fig 3J
source("scripts/dnds_blood_TP53_grouped_ms_3.R")
## Fig 4B
source("scripts/DNV_tissues_ms_4B.R")
## Fig 4I
source("scripts/dnds_classic_tissues_grouped_ms.R")

## Fig 4K
source("scripts/close_muts_181.R")
source("scripts/phasing_tissue_ms_4G.R")

## Fig S5
source("scripts/181_LFS_frequency_supp.R")



