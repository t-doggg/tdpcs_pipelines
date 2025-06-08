#!/usr/bin/env Rscript

# === Load libraries ===
library(dplyr)
library(stringr)
library(readr)

# === Hardcoded file paths ===
input_file <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_16_tdpcs_long_single/02-Classification/out_ncbi.tsv"
output_file <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_16_tdpcs_long_single/02-Classification/out_ncbi_collapsed.tsv"
low_conf_file <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_16_tdpcs_long_single/02-Classification/out_ncbi_low_confidence.tsv"

# === Read the input file ===
data <- read_tsv(input_file, col_names = FALSE, show_col_types = FALSE)

# Set column names
colnames(data) <- c("read_id", "sacc", "stitle", "ssciname", "nident", "qlen", "qcovs", "bitscore", "pident")

# === Helper: Extract first 2 words of stitle ===
extract_first_two_words <- function(x) {
  words <- str_split_fixed(x, " ", 3)
  str_trim(paste(words[, 1], words[, 2]))
}

data <- data %>%
  mutate(base_class = extract_first_two_words(stitle))

# === Group by read_id and apply consensus logic ===
collapsed <- data %>%
  group_by(read_id) %>%
  group_modify(~{
    df <- .
    n <- nrow(df)
    
    if (n < 3) {
      return(df[0, ])  # low confidence for n == 1 or 2
    }
    
    counts <- table(df$base_class)
    max_count <- max(counts)
    
    if (max_count >= 2) {
      # Majority class exists
      majority_class <- names(counts)[counts == max_count][1]
      selected <- df %>% filter(base_class == majority_class)
      return(selected[1, ])
    } else {
      return(df[0, ])  # no majority match
    }
  }) %>%
  ungroup() %>%
  select(-base_class)


# === Identify low-confidence reads ===
collapsed_ids <- collapsed$read_id
low_conf <- data %>% filter(!(read_id %in% collapsed_ids))

# === Write outputs ===
write_tsv(collapsed, output_file)
write_tsv(low_conf, low_conf_file)
