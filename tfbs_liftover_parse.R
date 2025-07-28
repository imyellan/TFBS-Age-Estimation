library(dplyr)
library(readr)
library(purrr)
library(dtplyr)
library(furrr)

# import arguments
psl_dir <- commandArgs(trailingOnly = TRUE)[1]
TF <- commandArgs(trailingOnly = TRUE)[2]
ancestral <- commandArgs(trailingOnly = TRUE)[3]

if (ancestral == "ancestral") {
  ancestral <- "_ancestral"
} else {
  ancestral <- NULL
}

# import other data
# print(paste0("Rscript /home/iyellan/scripts/tfbs_liftover_parse.R /home/iyellan/scratch/tfbs_ages/liftover_results ", TF))
psl_files <- list.files(psl_dir, full.names = T, pattern = paste0(TF, "\\..*.txt"))

if (is.null(ancestral)) {
  ages_df <- read_csv("~/scratch/tfbs_ages/human_cactus_genome_ages_timetree.csv") %>%
    transmute(
      species = gsub(" ", "_", scientific_name_b),
      species_age = precomputed_age
    )
} else if (ancestral == "_ancestral") {
  ages_df <- read_tsv("~/scratch/tfbs_ages/hum_ancestral_genomes.txt",
    col_names = c("species", "clade", "precomputed_age")
  ) %>%
    select(-clade)
}

TF_ages_df <- read_csv("~/scratch/tfbs_ages/all_codebook_TFs_ages_DBDs_darks.csv") %>%
  transmute(TF = `Gene name`, TF_age = max_age)

psl_colnames <- c(
  "tfbs_nm", "matches", "misMatches", "repMatches", "nCount", "qNumInsert",
  "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName",
  "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd",
  "blockCount", "blockSizes", "qStarts", "tStarts"
)

col_type_spec <- "ciiiiiiiicciiiciiiiccc"

import_psl <- function(psl_fil) {
  print(psl_fil)
  # if psl_fil has > 0 rows, import
  if (file.info(psl_fil)$size >= 1) {
    psl <- read_tsv(psl_fil,
      col_names = psl_colnames,
      col_types = col_type_spec,
      col_select = c(
        tfbs_nm, matches, repMatches, tBaseInsert,
        tName, tStart, tEnd, strand, blockCount
      )
    ) %>%
      mutate(
        TF = TF,
        species = gsub(".*\\.(.*)\\.txt", "\\1", basename(psl_fil)) # ,
        # qStarts = as.character(qStarts), tStarts = as.character(tStarts)
      )
    psl
  }
}

# if (!file.exists(file.path(psl_dir, paste0(TF, "_tfbs_liftover_ages.csv.gz")))) {
plan(multicore, workers = (parallel::detectCores() - 1))
psl_df <- psl_files %>% future_map_dfr(import_psl)
plan(sequential)
write_csv(psl_df, file.path(psl_dir, paste0(TF, ancestral, "_tfbs_liftover.csv.gz")))

# write targets to bed file for later
lft_bed_dir <- paste0("~/scratch/tfbs_ages/liftover_beds", ancestral)
dir.create(lft_bed_dir, showWarnings = F)
psl_df %>%
  transmute(species, tName, tStart, tEnd, tfbs_nm,
    score = 0,
    strand = sapply(strand, function(x) strsplit(x, "")[[1]][2])
  ) %>%
  group_by(species) %>%
  group_walk(~ write_tsv(
    x = .x %>% unique(),
    file = file.path(
      lft_bed_dir,
      paste0(
        TF, "_", .y, ancestral,
        ".bed"
      )
    ),
    col_names = F
  ))

psl_df <- lazy_dt(psl_df) %>%
  mutate(
    total_all_matches = matches + repMatches,
    tfbs_length = as.numeric(gsub(".*_", "", tfbs_nm)) -
      as.numeric(gsub(".*_([0-9]+)_.*", "\\1", tfbs_nm)),
    target_length = tEnd - tStart,
    frac_matches = total_all_matches / tfbs_length,
    target_length_near_tfbs =
    # two conditions: 1 block the length of the tfbs, or
    # 2 adjacent blocks (thanks to the liftover algorithm) which can be
    # detected by the blockCount and tBaseInsert
      (abs(tfbs_length - target_length) == 0 & blockCount == 1) |
        (abs(tfbs_length - target_length) == 1 & blockCount == 2 &
          tBaseInsert == 1)
  ) %>%
  select(
    species, TF, tfbs_nm, tfbs_length, frac_matches,
    target_length_near_tfbs
  ) %>%
  unique()

# psl_df_summ <-
#   psl_df %>% mutate(aln_length = tStart - tEnd) %>%
#   # group_by(TF, species, tfbs_nm, tName, tStart) %>%
#   transmute(total_matches = sum(matches),
#             total_mismatches = sum(misMatches),
#             total_repMatches = sum(repMatches),
#             total_nCount = sum(nCount),
#             total_all_matches = total_matches + total_repMatches,
#             tfbs_length = as.numeric(gsub(".*_", "", tfbs_nm)) -
#               as.numeric(gsub(".*_([0-9]+)_.*", "\\1", tfbs_nm)),
#             frac_matches = total_all_matches/tfbs_length) %>% unique()

psl_df_ages <- left_join(psl_df, ages_df) %>%
  relocate(TF, species, .before = tfbs_nm) # %>%
# as_tibble()
write_csv(
  psl_df_ages %>% as_tibble(),
  file.path(psl_dir, paste0(TF, ancestral, "_tfbs_liftover_ages.csv.gz"))
)
# } else {
# psl_df_ages <- read_csv(file.path(psl_dir, paste0(TF, "_tfbs_liftover_ages.csv.gz"))) #%>%
# if(!"species_age" %in% colnames(psl_df_ages)){
#   psl_df_ages <- psl_df_ages %>% mutate(species_age = age) %>%
#     select(-age)
# }
# }

# if (!file.exists(file.path(psl_dir, paste0(TF, "_tfbs_liftover_ages_tmp.csv.gz")))) {
psl_df_max_ages_tmp <-
  psl_df_ages %>%
  # lazy_dt() %>%
  # filter(tfbs_nm == "chr16_30957983_30957995") %>%
  left_join(
    lazy_dt(TF_ages_df)
  ) %>%
  # slice_sample(n = 1000000) %>%
  # filter(frac_matches >= 0.75) %>%
  mutate(
    frac_0 = target_length_near_tfbs,
    frac_075 = .data$frac_matches >= 0.75 & target_length_near_tfbs,
    frac_090 = .data$frac_matches >= 0.90 & target_length_near_tfbs,
    frac_095 = .data$frac_matches >= 0.95 & target_length_near_tfbs,
    frac_100 = .data$frac_matches >= 1 & target_length_near_tfbs,
    target_length_near_tfbs_num = as.numeric(target_length_near_tfbs)
  ) %>%
  # slice to wittle down multi-mapping hits to the single "best" hit per species
  group_by(tfbs_nm, species) %>%
  # partition(cluster) %>%
  slice_max(target_length_near_tfbs_num, with_ties = T) %>% # first prioritize nearness to TFBS length
  slice_max(frac_matches, n = 1, with_ties = F) %>% # then prioritize match %
  as_tibble()
write_csv(
  psl_df_max_ages_tmp,
  file.path(psl_dir, paste0(
    TF, ancestral,
    "_tfbs_liftover_ages_tmp.csv.gz"
  ))
)
# } else{
# psl_df_max_ages_tmp <- read_csv(file.path(psl_dir, paste0(TF, "_tfbs_liftover_ages_tmp.csv.gz")))
# }
