library(dplyr)
library(readr)
library(purrr)
library(ape)
library(multidplyr)
library(dtplyr)
library(stringr)
# library(Biostrings)
# library(universalmotif)

# import arguments
psl_dir <- commandArgs(trailingOnly = TRUE)[1]
TF <- commandArgs(trailingOnly = TRUE)[2]

# for adjusted date determination
mammal_dates <- c(
  Eutheria = 99.1887, Boreoeutheria = 94, Euarchontoglires = 87.2,
  Euarchonta = 84.73, Primatomorpha = 79.36470, Primates = 73.75492,
  Simiiformes = 42.9, Catarrhini = 28.82, Hominoidea = 19.45714,
  Hominidae = 15.2, Homininae = 8.6, Hominini = 6.4, Homo = 0
)
type_species <- c(
  Eutheria = "Loxodonta_africana", Boreoeutheria = "Sorex_araneus",
  Euarchontoglires = "Mus_musculus", Euarchonta = "Tupaia_tana",
  Primatomorpha = "Galeopterus_variegatus", Primates = "Propithecus_coquereli",
  Simiiformes = "Callithrix_jacchus", Catarrhini = "Macaca_mulatta",
  Hominoidea = "Nomascus_leucogenys", Hominidae = "Pongo_abelii",
  Homininae = "Gorilla_gorilla", Hominini = "Pan_troglodytes",
  Homo = "Homo_sapiens"
)

get_adjusted_age <- function(in_string, max_age) {
  # browser()
  if (is.infinite(max_age)) {
    return(-Inf)
  }
  for (taxa in names(type_species)[mammal_dates <= max_age]) {
    small_clade_flag <- T
    tax_spec <- extract.clade(
      tree, getMRCA(tree, c("Homo_sapiens", type_species[taxa]))
    )$tip.label
    if (taxa != "Homo") {
      tax_spec <- tax_spec[tax_spec != "Homo_sapiens"]
      num_spec <- length(tax_spec)
      matching_spec <- names(unlist(sapply(tax_spec, function(x) grep(x, in_string))))
      n_matching_spec <- length(matching_spec)
      # to account for clades with a small number of species sharing the MRCA with humans
      clade_spec <- ages_df %>%
        filter(species_age == mammal_dates[taxa]) %>%
        pull(species) %>%
        gsub(" ", "_", .)
      small_clade_flag <- any(matching_spec %in% clade_spec)
      # small_clade_flag <- case_when(
      #   taxa == "Eutheria" ~any(matching_spec %in% euth),
      #   taxa == "Euarchonta" ~ any(grepl("Tupaia", matching_spec)),
      #   taxa == "Primatomorpha" ~ any(grepl("Galeopterus", matching_spec)),
      #   taxa == "Hominoidea" ~ any(grepl("Nomascus", matching_spec)),
      #   taxa == "Hominidae" ~ any(grepl("Pongo", matching_spec)),
      #   taxa == "Homininae" ~ any(grepl("Gorilla", matching_spec)),
      #   taxa == "Hominini" ~ any(grepl("Pan", matching_spec)),
      #   TRUE ~ TRUE)
      # since there are only two chimp species, treating presence in one as sufficient
      if (taxa == "Hominini") {
        n_matching_spec <- 2
      }
    } else {
      tax_spec <- "Homo_sapiens"
      n_matching_spec <- 1
      num_spec <- 1
    }
    # browser()
    if (n_matching_spec / num_spec >= 0.6 & small_clade_flag) {
      return(mammal_dates[taxa])
      break
    }
  }
}

import_animal_fas <- function(fa_path) {
  if (file.info(fa_path)$size > 0) {
    species <- gsub(".*_([A-Za-z]+_[A-Za-z]+).fa", "\\1", basename(fa_path))
    dna_str <- Biostrings::readDNAStringSet(fa_path)
    data.frame(
      species = species,
      tfbs_nm = names(dna_str),
      lftover_seq = as.data.frame(dna_str)$x
    ) %>% unique()
  }
}

compare_strings <- function(str1, str2) {
  if (!is.na(str1) & !is.na(str2) & nchar(str1) == nchar(str2)) {
    Biostrings::compareStrings(str1, str2)
  } else {
    return("")
  }
}

# import other data
psl_df_max_ages_tmp <- read_csv(file.path(psl_dir, paste0(TF, "_tfbs_liftover_ages_tmp.csv.gz"))) %>% lazy_dt()
ages_df <- read_csv("~/scratch/tfbs_ages/human_cactus_genome_ages_timetree.csv")
ages_df <- ages_df %>%
  transmute(
    species = gsub(" ", "_", scientific_name_b),
    species_age = precomputed_age
  )
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
PWM <- universalmotif::read_cisbp(list.files("~/scratch/tfbs_ages/motifs",
  pattern = paste0(TF, "_"), full.names = T
))
human_motif_hits <-
  Biostrings::readDNAStringSet(list.files("~/scratch/tfbs_ages/tmp_fas",
    pattern = paste0(TF, "_"), full.names = T
  ))
human_motif_hits <- data.frame(
  tfbs_nm = names(human_motif_hits),
  hum_seq = as.data.frame(human_motif_hits)$x
) %>%
  unique()
animal_motif_hits <- list.files("~/scratch/tfbs_ages/liftover_beds/",
  pattern = paste0(TF, "_.*fa$"), full.names = T
) %>%
  map(import_animal_fas) %>%
  list_rbind()

## deal with issue of matches being undercounted due to softmasking
### extract sequences from genome
motif_hit_df <- full_join(human_motif_hits, animal_motif_hits) %>%
  mutate(
    tfbs_nm.species =
      paste0(tfbs_nm, ".", species, ".", lftover_seq)
  )
liftover_pwm_scores <-
  as.data.frame(universalmotif::scan_sequences(
    universalmotif::convert_type(PWM, "PWM", pseudocount = 1e-3),
    Biostrings::DNAStringSet(
      motif_hit_df %>% filter(nchar(lftover_seq) == nchar(hum_seq)) %>%
        select(tfbs_nm.species, lftover_seq) %>% unique() %>% tibble::deframe()
    ),
    threshold = 1
  )@listData)
motif_hit_df <- motif_hit_df %>%
  mutate(
    seq_compare = mapply(compare_strings, hum_seq, lftover_seq),
    frac_matches_new = 1 - (str_count(seq_compare, "\\?") / nchar(seq_compare))
  ) %>%
  left_join(tibble(
    tfbs_nm.species = liftover_pwm_scores$sequence,
    lftover_pwm_score = liftover_pwm_scores$score
  )) %>%
  mutate(lftover_score_over_thresh = lftover_pwm_score >= trip_opt_thresh)

psl_df_max_ages_tmp <- psl_df_max_ages_tmp %>%
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

# detach("package:GenomeInfoDb", unload=TRUE)
# unloadNamespace("GenomeInfoDb")
# detach("package:IRanges", unload=TRUE)
# unloadNamespace("IRanges")
# detach("package:BiocGenerics", unload=TRUE)
# unloadNamespace("BiocGenerics")
# detach("package:Biostrings", unload=TRUE)
# unloadNamespace("Biostrings")

# join with PWM scores
# psl_df_max_ages_tmp_pwm_scores <- psl_df_max_ages_tmp %>%
#   left_join(lazy_dt(PWM_scores) %>% select(tfbs_nm, PWM_score)) %>% as_tibble()

## Summarize ages
cluster <- new_cluster(parallel::detectCores() - 1)
cluster_copy(cluster, c("get_adjusted_age", "type_species", "mammal_dates", "tree", "ages_df"))
cluster_library(cluster, c("dplyr", "ape"))

psl_df_max_ages <-
  psl_df_max_ages_tmp %>%
  # select(-frac_matches_new) %>%
  group_by(tfbs_nm) %>%
  partition(cluster) %>%
  summarize(
    n_species_cons_0 = dplyr::n_distinct(.data$species[frac_0]),
    n_species_cons_075 = dplyr::n_distinct(.data$species[frac_075]),
    n_species_cons_090 = dplyr::n_distinct(.data$species[frac_090]),
    n_species_cons_100 = dplyr::n_distinct(.data$species[frac_100]),
    n_species_cons_pwm = dplyr::n_distinct(.data$species[lftover_score_over_thresh]),

    # n_species_cons_075_TF_informed_age =
    #   dplyr::n_distinct(.data$species[frac_075 & (TF_age >= species_age)]),
    # n_species_cons_090_TF_informed_age =
    #   dplyr::n_distinct(.data$species[frac_090 & (TF_age >= species_age)]),
    # n_species_cons_100_TF_informed_age =
    #   dplyr::n_distinct(.data$species[frac_100 & (TF_age >= species_age)]),
    max_age_0 = max(.data$species_age[frac_0], na.rm = T),
    max_age_075 = max(.data$species_age[frac_075], na.rm = T),
    max_age_090 = max(.data$species_age[frac_090], na.rm = T),
    max_age_100 = max(.data$species_age[frac_100], na.rm = T),
    max_age_pwm = max(.data$species_age[lftover_score_over_thresh], na.rm = T),
    max_age_species_0 =
      paste(unique(.data$species[which(.data$species_age == max_age_0 & frac_0)]),
        collapse = ";"
      ),
    max_age_species_075 =
      paste(unique(.data$species[which(.data$species_age == max_age_075 & frac_075)]),
        collapse = ";"
      ),
    max_age_species_090 =
      paste(unique(.data$species[which(.data$species_age == max_age_090 & frac_090)]),
        collapse = ";"
      ),
    max_age_species_100 =
      paste(unique(.data$species[which(.data$species_age == max_age_100 & frac_100)]),
        collapse = ";"
      ),
    max_age_species_pwm =
      paste(
        unique(.data$species[which(.data$species_age == max_age_pwm &
          lftover_score_over_thresh)]),
        collapse = ";"
      ),
    species_0_all =
      paste(unique(.data$species[which(.data$species_age <= max_age_0 & frac_0)]),
        collapse = ";"
      ),
    species_075_all =
      paste(unique(.data$species[which(.data$species_age <= max_age_075 & frac_075)]),
        collapse = ";"
      ),
    species_090_all =
      paste(unique(.data$species[which(.data$species_age <= max_age_090 & frac_090)]),
        collapse = ";"
      ),
    species_100_all =
      paste(unique(.data$species[which(.data$species_age <= max_age_100 & frac_100)]),
        collapse = ";"
      ),
    species_pwm_all =
      paste(
        unique(.data$species[which(.data$species_age <= max_age_pwm &
          lftover_score_over_thresh)]),
        collapse = ";"
      ),
    # not using TF ages anymore
    # max_age_075_TF_informed_age = max(.data$species_age[frac_075 & (TF_age >= species_age)], na.rm = T),
    # max_age_090_TF_informed_age = max(.data$species_age[frac_090 & (TF_age >= species_age)], na.rm = T),
    # max_age_100_TF_informed_age = max(.data$species_age[frac_100 & (TF_age >= species_age)], na.rm = T),
    # max_age_species_075_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age == max_age_075_TF_informed_age & frac_075 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    # max_age_species_090_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age == max_age_090_TF_informed_age & frac_090 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    # max_age_species_100_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age == max_age_100_TF_informed_age & frac_100 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    # species_075_all_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age <= max_age_075_TF_informed_age & frac_075 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    # species_090_all_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age <= max_age_090_TF_informed_age & frac_090 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    # species_100_all_TF_informed_age =
    #   paste(
    #     unique(.data$species[which(.data$species_age <= max_age_100_TF_informed_age & frac_100 &
    #                                  (TF_age >= species_age))]),
    #     collapse = ";"
    #   ),
    adjusted_age_0 = unlist(mapply(get_adjusted_age, species_0_all, max_age_0)),
    adjusted_age_075 = unlist(mapply(get_adjusted_age, species_075_all, max_age_075)),
    adjusted_age_090 = unlist(mapply(get_adjusted_age, species_090_all, max_age_090)),
    adjusted_age_100 = unlist(mapply(get_adjusted_age, species_100_all, max_age_100)),
    adjusted_age_pwm = unlist(mapply(get_adjusted_age, species_pwm_all, max_age_pwm)),
    species_0_all_adjusted_age =
      paste(
        unique(.data$species[which(.data$species_age <= adjusted_age_0 & frac_0 &
          (TF_age >= species_age))]),
        collapse = ";"
      ),
    species_075_all_adjusted_age =
      paste(
        unique(.data$species[which(.data$species_age <= adjusted_age_075 & frac_075 &
          (TF_age >= species_age))]),
        collapse = ";"
      ),
    species_090_all_adjusted_age =
      paste(
        unique(.data$species[which(.data$species_age <= adjusted_age_090 & frac_090 &
          (TF_age >= species_age))]),
        collapse = ";"
      ),
    species_100_all_adjusted_age =
      paste(
        unique(.data$species[which(.data$species_age <= adjusted_age_100 & frac_100 &
          (TF_age >= species_age))]),
        collapse = ";"
      ),
    species_pwm_all_adjusted_age =
      paste(
        unique(.data$species[which(.data$species_age <= adjusted_age_pwm &
          lftover_score_over_thresh)]),
        collapse = ";"
      )
  ) %>%
  collect() %>%
  unique() %>% # select(-frac_075, -frac_090, -frac_100) %>%
  # as_tibble() %>%
  # fix infinite values
  mutate(across(everything(), ~ ifelse(is.infinite(.), NA, .)))

## write files
write_csv(psl_df_max_ages, file.path(psl_dir, paste0(TF, "_tfbs_liftover_max_ages.csv.gz")))
