#!/usr/bin/env Rscript

cat("Running script summa_bacteria_kraken2.R\n")

# === Command Line Argument Handling ===
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript summa_bacteria_kraken2.R <kraken2_report_csv> <output_csv>")
}

kraken2_csv <- args[1]
output_csv <- args[2]

# === Load Libraries ===
library(dplyr)
library(readr)
library(stringr)

# === Parameters ===
valid_ranks <- c("S", "S1")    # Species and strain level entries only

# === Load Kraken2 Report CSV ===
data <- read_delim(
  kraken2_csv,
  delim = ",",
  col_names = c("percent_reads", "num_reads", "num_reads_direct", "rank_code", "ncbi_taxid", "name"),
  trim_ws = TRUE,
  show_col_types = FALSE
)

# === Clean up and parse ===
data <- data %>%
  mutate(
    name = str_squish(name),
    num_reads = as.numeric(num_reads),
    rank_code = str_trim(rank_code),
    percent_reads = as.numeric(percent_reads)
  )

# === Filter and summarize ===
result <- data %>%
  filter(percent_reads > 0.4, rank_code %in% valid_ranks) %>%
  rename(bacteria = name, Count = num_reads) %>%
  group_by(bacteria) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Count))

# === Write Output CSV ===
write_csv(result, output_csv)
