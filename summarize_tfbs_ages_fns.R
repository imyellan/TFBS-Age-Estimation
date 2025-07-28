#functions for running summarize tfbs_ages

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
    if(grepl("ancestral", fa_path)){
      species <- gsub(".*_(fullTreeAnc[A-Za-z0-9]+)_ancestral.fa",
                      "\\1", basename(fa_path))
    } else{
      species <- gsub(".*_([A-Za-z]+_[A-Za-z]+).fa", "\\1", basename(fa_path))
    }
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

moods_replica <- function(moods_pwm, dna_seq, hum_length) {
  if (!is.na(dna_seq)) {
    if (!grepl("N|-", dna_seq) & nchar(dna_seq) == hum_length) {
      dna_seq_split <- strsplit(dna_seq, "")[[1]]
      # browser()
      sum(sapply(1:nchar(dna_seq), function(x) moods_pwm[dna_seq_split[x],x]))
    } else{return(NA)}
  } else{return(NA)}
}

traverse_ancestral <- function(df_filt, tfbs_nm, frac_thresh){
  # df_filt <- df %>% filter(tfbs_nm/ == !!tfbs_nm)
  ancest_age_df <- data.frame()
  ancest_specs <- ages_df$ancestral_species
  for (thresh in frac_thresh) {
    for (i in 1:nrow(df_filt)) {
      pass <- df_filt %>% filter(species == ancest_specs[i]) %>% 
        pull(.data[[thresh]])
      # print(paste(tfbs_nm, i))
      # if(tfbs_nm == "chr11_1806386_1806406" & i == 6){
      #   browser()
      # }
      # check if last species in df_filt has been reached
      # check if tfbs was aligned to species - if not pass will be NULL and length 0,
      if (length(pass) == 0 || is.na(pass)) {
        # if i == 1, then no alignment with any ancestral genome, so age is 0
        passing_clade <- ifelse(i > 1, ages_df$clade[(i - 1)], "Homo")
        passing_age <- ifelse(i > 1, ages_df$precomputed_age[(i - 1)], 0)
        break
      } else if (!pass & i > 1) {
        passing_clade <- ages_df$clade[(i - 1)]
        passing_age <- ages_df$precomputed_age[(i - 1)]
        break
      } else if (i == nrow(df_filt)) {
        passing_clade <- ifelse(pass, ages_df$clade[i], 
                                ages_df$clade[(i - 1)])
        passing_age <- ifelse(pass, 
                              ages_df$precomputed_age[i], 
                              ages_df$precomputed_age[(i - 1)])
      }
    }
    ancest_age_df <- rbind(ancest_age_df, 
                           data.frame(tfbs_nm = tfbs_nm, frac_thresh = thresh,
             ancestral_clade = passing_clade, ancestral_age = passing_age))
  }
  ancest_age_df
}