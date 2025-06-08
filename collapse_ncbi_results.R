#!/usr/bin/env Rscript

# === Load libraries ===
library("dplyr")
library("stringr")
library("readr")

# === Parse command-line arguments ===
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript collapse_bacteria.R <input_tsv> <output_tsv>")
}

input_file <- args[1]
output_file <- args[2]
low_conf_file <- args[3]

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
    
    counts <- table(df$base_class)
    max_count <- max(counts)
    
    if ((n == 1) ||
        ((n == 2 || n == 3) && max_count >= 2) ||
        ((n == 4 || n == 5) && max_count >= 2)) {
      
      majority_class <- names(counts)[counts == max_count][1]
      selected <- df %>% filter(base_class == majority_class)
      return(selected[1, ])
    } else {
      return(df[0, ])  # no consensus
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
