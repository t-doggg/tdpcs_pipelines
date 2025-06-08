#!/usr/bin/env Rscript

cat("Running script compare_kraken_ncbi.R\n")

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compare_bacteria.R <csv1_path> <csv2_path> <output_path>")
}

csv1_path <- args[1]
csv2_path <- args[2]
output_path <- args[3]

# Debug output
print(paste("CSV Path:", csv1_path))
print(paste("CSV Path:", csv2_path))
print(paste("Output Path:", output_path))

# --- Load required libraries ---
library(dplyr)
library(stringr)
library(readr)

# --- Read CSV files ---
csv1 <- read_csv(csv1_path, show_col_types = FALSE)
csv2 <- read_csv(csv2_path, show_col_types = FALSE)

# --- Standardize column names and clean entries ---
colnames(csv1) <- c("bacteria", "count_1")
colnames(csv2) <- c("bacteria", "count_2")

csv1 <- csv1 %>%
  mutate(bacteria = str_squish(bacteria))

csv2 <- csv2 %>%
  mutate(bacteria = str_squish(bacteria))

# --- Merge datasets on 'bacteria' ---
merged <- full_join(csv1, csv2, by = "bacteria")

# --- Replace NA counts with 0 ---
merged <- merged %>%
  mutate(
    count_1 = ifelse(is.na(count_1), 0, count_1),
    count_2 = ifelse(is.na(count_2), 0, count_2),
    confidence = ifelse(count_1 > 0 & count_2 > 0, 1, 0)
  )

# --- Write output ---
write_csv(merged, output_path)
cat("Comparison saved to:", output_path, "\n")
