#!/usr/bin/env Rscript

cat("Running script summa_bacteria_ncbi.R\n")

# === Command Line Argument Handling ===
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript summa_bacteria.R <input_csv> <output_csv> <synonym_file>")
}

input_csv <- args[1]
output_csv <- args[2]
synonym_file <- args[3]

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
