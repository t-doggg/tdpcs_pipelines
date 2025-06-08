#!/usr/bin/env Rscript

# === Load Libraries ===
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# === Parse Arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript summarize_armfinder_hits.R <input_tsv> <output_summary_csv>")
}

input_file <- args[1]
output_file <- args[2]

# === Parameters ===
min_identity <- 90      # Minimum % identity
min_coverage <- 50      # Minimum % coverage

# === Read Input File ===
df <- read_tsv(input_file, show_col_types = FALSE, comment = "#")

# === Filter High-Confidence Hits ===
filtered <- df %>%
  filter(`% Identity to reference sequence` >= min_identity,
         `% Coverage of reference sequence` >= min_coverage)

# === Summarize by Class and Subclass ===
summary_table <- filtered %>%
  count(Class, Subclass, name = "count") %>%
  arrange(desc(count))

# === Write Output as CSV ===
write_csv(summary_table, output_file)

cat("✅ Resistance summary written to:", output_file, "\n")
