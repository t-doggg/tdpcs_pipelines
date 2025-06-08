#!/usr/bin/env Rscript

input_csv <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_05_tdpcs_long_single/03-Results/f_res_ncbi.csv"
output_csv <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_05_tdpcs_long_single/03-Results/s_res_ncbi.csv"
synonym_file <- "/home/diablo/publ_pipeline/others/bacteria_synonyms.txt"

# === Load Libraries ===
library(dplyr)
library(stringr)
library(readr)

# === Load CSV ===
data <- read_csv(input_csv, show_col_types = FALSE)

# === Load Synonym Mappings ===
synonym_lines <- readLines(synonym_file)
synonym_map <- list()

for (line in synonym_lines) {
  parts <- str_split(line, ",")[[1]]
  canonical <- str_trim(parts[1])
  synonyms <- str_trim(parts)
  for (syn in synonyms) {
    synonym_map[[tolower(syn)]] <- canonical
  }
}

# === Map Entries to Canonical Names using substring matching ===
data <- data %>%
  rename(bacteria = Stamm) %>%
  rowwise() %>%
  mutate(
    match_key = {
      matched <- NA_character_
      bac <- tolower(bacteria)
      if (!is.na(bac)) {
        for (syn in names(synonym_map)) {
          if (!is.na(syn) && nzchar(syn) && str_detect(bac, fixed(syn, ignore_case = TRUE))) {
            matched <- synonym_map[[syn]]
            break
          }
        }
      }
      if (!is.na(matched)) matched else bacteria
    }
  ) %>%
  ungroup()

# === Aggregate Counts ===
result <- data %>%
  group_by(match_key) %>%
  summarise(total_count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  rename(bacteria = match_key)

# === Write Output ===
write_csv(result, output_csv)
