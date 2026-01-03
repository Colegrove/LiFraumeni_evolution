################################################################################
#### create analysis .bed (overlapping probe and cds regions)
################################################################################

## First run bed_file_intersect_coverage.sh to generate cds_baits_intersect.bed
## this finds the actual position ranges overlapping between cds and baits
## generates coding_targets.bed for use as analysis targets (coding/non-coding)

intersect_path <- "/Volumes/feder-vol1/project/li_fraumeni/
scripts/2025-10-02-bed_file_prep/cds_baits_intersect.bed"

columns <- c("cds_chr", "cds_start", "cds_end", "transcript_ID", "id", "strand",
             "probe_chr", "probe_start", "probe_end", "bp_overlap")
overlap_df <- read_delim(intersect_path, col_names = columns) %>%
  mutate(overlap_start = pmax(cds_start, probe_start),
         overlap_end = pmin(cds_end, probe_end))

out_file <- file.path("/Volumes/feder-vol1/project/li_fraumeni/scripts/",
            "2025-10-02-bed_file_prep/",
            "coding_targets.bed")

write_tsv(
  overlap_df %>% select(cds_chr, overlap_start, overlap_end),
  out_file,
  col_names = FALSE
)
