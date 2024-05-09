library(dplyr)
library(readr)
library(purrr)
library(ape)
# library(multidplyr)
library(dtplyr)
library(stringr)
library(tidyr)
library(furrr)
# library(Biostrings)
# library(universalmotif)
source("~/tfbs_age_estimation/summarize_tfbs_ages_fns.R")

# import arguments
psl_dir <- commandArgs(trailingOnly = TRUE)[1]
TF <- commandArgs(trailingOnly = TRUE)[2]

mammal_dates <- c(
  Eutheria = 99.1887, Boreoeutheria = 94, Euarchontoglires = 87.2,
  Euarchonta = 84.73, Primatomorpha = 79.36470, Primates = 73.75492,
  Simiiformes = 42.9, Catarrhini = 28.82, Hominoidea = 19.45714,
  Hominidae = 15.2, Homininae = 8.6, Hominini = 6.4, Homo = 0
)

# import other data
psl_df_max_ages_tmp <-
  read_csv(file.path(psl_dir, 
                     paste0(TF, "_ancestral_tfbs_liftover_ages_tmp.csv.gz"))) %>% 
  lazy_dt()

ages_df <- read_tsv("~/scratch/tfbs_ages/hum_ancestral_genomes.txt", 
                    col_names = c("ancestral_species", "clade", "precomputed_age"))
tree <- read.tree("~/scratch/tfbs_ages/241-mammalian-2020v2.nwk")
pwm_trip_threshes <- read_tsv("~/scratch/tfbs_ages/PWM_Triple_Thresholds.txt",
                              col_names = c("TF", "PWM_ID", "threshold")
)
trip_opt_thresh <- pwm_trip_threshes %>%
  filter(TF == !!TF) %>%
  pull(threshold)
PWM_scores <- read_tsv(
  list.files("~/scratch/tfbs_ages/tmp_beds",
             pattern = paste0(TF, "_renamed.bed"),
             full.names = T
  ),
  col_names = c("chr", "start", "end", "tfbs_nm", "PWM_score", "strand")
)
PWM <- universalmotif::read_matrix(
  list.files("~/scratch/tfbs_ages/motifs/pfms_ali", 
             pattern = paste0(TF, "_"), full.names = T), 
  type = "PCM", sep = "\t", headers = F)
pwm_mat <- PWM@motif
in_mat_colsums <- colSums(pwm_mat + 0.01*0.25)
in_mat_freq <- (pwm_mat + 0.01*0.25)/in_mat_colsums
moods_pwm <- log(in_mat_freq/0.25)

human_motif_hits <-
  Biostrings::readDNAStringSet(list.files("~/scratch/tfbs_ages/tmp_fas",
                                          pattern = paste0(TF, "_"), full.names = T
  ))
human_motif_hits <- data.frame(
  tfbs_nm = names(human_motif_hits),
  hum_seq = as.data.frame(human_motif_hits)$x
) %>%
  unique()
animal_motif_hits <-
  list.files("~/scratch/tfbs_ages/liftover_beds_ancestral/",
             pattern = paste0(TF, "_.*fa$"), full.names = T) %>%
  map(import_animal_fas) %>% list_rbind()

## deal with issue of matches being undercounted due to softmasking
### extract sequences from genome
motif_hit_df <- full_join(human_motif_hits, animal_motif_hits) %>%
  mutate(
    tfbs_nm.species =
      paste0(tfbs_nm, ".", species, ".", lftover_seq)
  )

hum_length <- nchar(head(motif_hit_df$hum_seq, n = 1))
# start.time <- Sys.time()

## scan with motif
motif_hit_df <- motif_hit_df %>% 
  rowwise() %>%
  mutate(lftover_pwm_score = 
           moods_replica(moods_pwm, lftover_seq, hum_length)) %>%
  ungroup()

## then simple identity calculation
motif_hit_df <- motif_hit_df %>%
  mutate(
    seq_compare = mapply(compare_strings, hum_seq, lftover_seq),
    frac_matches_new = 1 - (str_count(seq_compare, "\\?") / nchar(seq_compare))
  ) %>%
  mutate(lftover_score_over_thresh = lftover_pwm_score >= trip_opt_thresh)

psl_df_max_ages_tmp <- psl_df_max_ages_tmp  %>%
  left_join(lazy_dt(motif_hit_df) %>%
              select(species, tfbs_nm, frac_matches_new, lftover_score_over_thresh) %>%
              group_by(species, tfbs_nm) %>%
              slice_max(unname(frac_matches_new), with_ties = F)) %>%
  mutate(
    frac_0 = ifelse(!is.na(frac_matches_new), 
                    frac_matches_new >= 0 & target_length_near_tfbs, frac_0),
    frac_075 = ifelse(!is.na(frac_matches_new), 
                      frac_matches_new >= 0.75 & target_length_near_tfbs, frac_075),
    frac_090 = ifelse(!is.na(frac_matches_new), 
                      frac_matches_new >= 0.90 & target_length_near_tfbs, frac_090),
    frac_095 = ifelse(!is.na(frac_matches_new), 
                      frac_matches_new >= 0.95 & target_length_near_tfbs, F),
    frac_100 = ifelse(!is.na(frac_matches_new), 
                      frac_matches_new == 1 & target_length_near_tfbs, frac_100)
  ) %>%
  as_tibble()


## Summarize ages
# cluster <- new_cluster(parallel::detectCores() - 1)
# cluster_copy(cluster, c("get_adjusted_age", "type_species", "mammal_dates", "tree", "ages_df"))
# cluster_library(cluster, c("dplyr", "ape"))

## 
plan(multicore)
ancestral_genome_ages <- psl_df_max_ages_tmp %>%
  as_tibble() %>% 
  group_by(tfbs_nm) %>%
  group_split() %>%
  future_map(~traverse_ancestral(
    .x, unique(.x$tfbs_nm),
    grep("frac_[10]|lftover_score_over_thresh",
         colnames(psl_df_max_ages_tmp), value = T)
    )) %>%
  list_rbind() %>%
  pivot_wider(names_from = "frac_thresh", 
              values_from = c(ancestral_age, ancestral_clade))
plan(sequential)
psl_df_max_ages <- psl_df_max_ages_tmp %>% left_join(ancestral_genome_ages)

# psl_df_max_ages <-
#   psl_df_max_ages_tmp %>%
#   group_by(tfbs_nm) %>%
#   partition(cluster) %>%
#   summarize(
#     n_species_cons_0 = dplyr::n_distinct(.data$species[frac_0]),
#     n_species_cons_075 = dplyr::n_distinct(.data$species[frac_075]),
#     n_species_cons_090 = dplyr::n_distinct(.data$species[frac_090]),
#     n_species_cons_100 = dplyr::n_distinct(.data$species[frac_100]),
#     n_species_cons_pwm = dplyr::n_distinct(.data$species[lftover_score_over_thresh]),
#    
#      max_age_0 = max(.data$species_age[frac_0], na.rm = T),
#     max_age_075 = max(.data$species_age[frac_075], na.rm = T),
#     max_age_090 = max(.data$species_age[frac_090], na.rm = T),
#     max_age_100 = max(.data$species_age[frac_100], na.rm = T),
#     max_age_pwm = max(.data$species_age[lftover_score_over_thresh], na.rm = T),
#     max_age_species_0 =
#       paste(unique(.data$species[which(.data$species_age == max_age_0 & frac_0)]),
#             collapse = ";"
#       ),
#     max_age_species_075 =
#       paste(unique(.data$species[which(.data$species_age == max_age_075 & frac_075)]),
#             collapse = ";"
#       ),
#     max_age_species_090 =
#       paste(unique(.data$species[which(.data$species_age == max_age_090 & frac_090)]),
#             collapse = ";"
#       ),
#     max_age_species_100 =
#       paste(unique(.data$species[which(.data$species_age == max_age_100 & frac_100)]),
#             collapse = ";"
#       ),
#     max_age_species_pwm =
#       paste(
#         unique(.data$species[which(.data$species_age == max_age_pwm &
#                                      lftover_score_over_thresh)]),
#         collapse = ";"
#       ),
#     
#     species_0_all =
#       paste(unique(.data$species[which(.data$species_age <= max_age_0 & frac_0)]),
#             collapse = ";"
#       ),
#     species_075_all =
#       paste(unique(.data$species[which(.data$species_age <= max_age_075 & frac_075)]),
#             collapse = ";"
#       ),
#     species_090_all =
#       paste(unique(.data$species[which(.data$species_age <= max_age_090 & frac_090)]),
#             collapse = ";"
#       ),
#     species_100_all =
#       paste(unique(.data$species[which(.data$species_age <= max_age_100 & frac_100)]),
#             collapse = ";"
#       ),
#     species_pwm_all =
#       paste(
#         unique(.data$species[which(.data$species_age <= max_age_pwm &
#                                      lftover_score_over_thresh)]),
#         collapse = ";"
#       ),
#     
#     # not using TF ages anymore
#     adjusted_age_0 = unlist(mapply(get_adjusted_age, species_0_all, max_age_0)),
#     adjusted_age_075 = unlist(mapply(get_adjusted_age, species_075_all, max_age_075)),
#     adjusted_age_090 = unlist(mapply(get_adjusted_age, species_090_all, max_age_090)),
#     adjusted_age_100 = unlist(mapply(get_adjusted_age, species_100_all, max_age_100)),
#     adjusted_age_pwm = unlist(mapply(get_adjusted_age, species_pwm_all, max_age_pwm)),
#     species_0_all_adjusted_age =
#       paste(
#         unique(.data$species[which(.data$species_age <= adjusted_age_0 & frac_0 &
#                                      (TF_age >= species_age))]),
#         collapse = ";"
#       ),
#     species_075_all_adjusted_age =
#       paste(
#         unique(.data$species[which(.data$species_age <= adjusted_age_075 & frac_075 &
#                                      (TF_age >= species_age))]),
#         collapse = ";"
#       ),
#     species_090_all_adjusted_age =
#       paste(
#         unique(.data$species[which(.data$species_age <= adjusted_age_090 & frac_090 &
#                                      (TF_age >= species_age))]),
#         collapse = ";"
#       ),
#     species_100_all_adjusted_age =
#       paste(
#         unique(.data$species[which(.data$species_age <= adjusted_age_100 & frac_100 &
#                                      (TF_age >= species_age))]),
#         collapse = ";"
#       ),
#     species_pwm_all_adjusted_age =
#       paste(
#         unique(.data$species[which(.data$species_age <= adjusted_age_pwm &
#                                      lftover_score_over_thresh)]),
#         collapse = ";"
#       )
#   ) %>%
#   collect() %>%
#   unique() %>%
#   # fix infinite values
#   mutate(across(everything(), ~ ifelse(is.infinite(.), NA, .)))

## write files
write_csv(psl_df_max_ages, 
          file.path(psl_dir, paste0(TF, "_ancestral_tfbs_liftover_max_ages.csv.gz")))

