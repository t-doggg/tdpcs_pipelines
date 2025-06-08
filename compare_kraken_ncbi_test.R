#!/usr/bin/env Rscript

cat("Running script compare_kraken_ncbi.R\n")

csv1_path <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_05_tdpcs_long_single/03-Results/s_res_kraken2.csv"
csv2_path <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_05_tdpcs_long_single/03-Results/s_res_ncbi.csv"
output_path <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_05_tdpcs_long_single/03-Results/s_res_compared_conficence.csv"

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
