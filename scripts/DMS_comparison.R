## Hunter Colegrove
## 10 Jul 2025
## Li-Fraumeni

## Compare DMS data from TP53 DBD (https://www.nature.com/articles/s41588-024-02039-4#Sec27)
## to alphamissense pathogenicity


DMS_df <- read_delim("inputs/DMS_TP53/DMS_TP53_DBD.txt")

columns_to_filter <- c("hg38_genomic", "hg38_protein", "library_id", "location", "type_g", "effect", 
                       "codon", "codon_ref", "codon_alt", "aa_ref", "aa_alt", "rfs_median")
#columns_to_merge <- c("", )
DMS <- DMS_df %>% dplyr::select(all_of(columns_to_filter)) %>% print(width = Inf, n=100)

DMS
alphamissense
AM <- alphamissense %>% mutate(effect = paste0("p.", protein_variant))

DMS_AM <- DMS %>% left_join(AM, by="effect")
DMS_AM <- DMS_AM %>% filter(!(is.na(am_pathogenicity))) %>% print(width = Inf)


background_df_sections <- tibble(
  xmin = c(-Inf, -Inf, 0.33, 0.33, 0.564, 0.564),
  xmax = c(0.33, 0.33, 0.564, 0.564, Inf, Inf),
  ymin = c(-Inf, 0, -Inf, 0, -Inf, 0),
  ymax = c(0, Inf, 0, Inf, 0, Inf),
  category = factor(c(
    "AM Benign & DMS Low RFS", "AM Benign & DMS High RFS",
    "AM Ambiguous & DMS Low RFS", "AM Ambiguous & DMS High RFS",
    "AM Pathogenic & DMS Low RFS", "AM Pathogenic & DMS High RFS"
  ), levels = c(
    "AM Benign & DMS Low RFS", "AM Benign & DMS High RFS",
    "AM Ambiguous & DMS Low RFS", "AM Ambiguous & DMS High RFS",
    "AM Pathogenic & DMS Low RFS", "AM Pathogenic & DMS High RFS"
  ))
)

background_colors <- c(
  "AM Benign & DMS Low RFS"     = "#a6cee3",
  "AM Benign & DMS High RFS"    = "#1f78b4", 
  "AM Ambiguous & DMS Low RFS"  = "#fdbf6f", 
  "AM Ambiguous & DMS High RFS" = "#ff7f00",  
  "AM Pathogenic & DMS Low RFS" = "#fdbbd2",  
  "AM Pathogenic & DMS High RFS"= "#e31a1c"  
)

ggplot(DMS_AM, aes(x = am_pathogenicity, y = rfs_median)) + 
  geom_rect(data = background_df_sections,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category),
            alpha = 1, color = NA, inherit.aes = FALSE) +
  
  geom_point() + 
  scale_fill_manual(values = background_colors) 
