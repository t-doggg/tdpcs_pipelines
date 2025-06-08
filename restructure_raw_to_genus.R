#!/usr/bin/env Rscript

cat("Running script restructure_raw_to_genus.R\n")

# Load required libraries
library(dplyr)
library(stringr)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check arguments
if (length(args) < 2) {
  stop("Bitte sowohl den Pfad zur Eingabedatei als auch den Pfad zur Ausgabedatei als Argumente angeben.")
}

input_file <- args[1]
output_file <- args[2]

# Debug output
print(paste("CSV Path:", input_file))
print(paste("Output Path:", output_file))

# Generate temporary file
temp_output_file <- tempfile(fileext = ".tsv")

# Function to process bacteria species counts
process_bacteria_counts_file <- function(input_file, temp_output_file) {
  
  # Read TSV with header
  data <- read.delim(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Check for enough columns
  if (ncol(data) < 3) {
    stop("Input file hat nicht mindestens drei Spalten.")
  }
  
  # Extract third column (classification) and clean brackets
  bacteria_column_clean <- str_replace_all(data[[3]], "\\[|\\]", "")
  
  # Extract first two words
  bacteria_names <- sapply(strsplit(bacteria_column_clean, " "), function(x) {
    if (length(x) >= 2) {
      paste(x[1:2], collapse = " ")
    } else {
      x[1]
    }
  })
  
  # Count occurrences
  bacteria_counts <- as.data.frame(table(bacteria_names))
  colnames(bacteria_counts) <- c("Bacteria", "Count")
  
  # Write temporary output
  write.table(bacteria_counts, file = temp_output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  print(paste("Bacteria counts saved to:", temp_output_file))
}

# Function to convert TSV to final CSV
convert_tsv_to_csv <- function(tsv_file, csv_file) {
  output_dir <- dirname(csv_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  data <- read.delim(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(data) <- c("Stamm", "Count")
  
  write.csv(data, file = csv_file, row.names = FALSE, quote = FALSE)
  print(paste("CSV file has been saved to", csv_file))
}

# Run script steps
process_bacteria_counts_file(input_file, temp_output_file)

# Convert TSV to CSV
csv_output_file <- sub("\\.tsv$", ".csv", output_file)
convert_tsv_to_csv(temp_output_file, csv_output_file)

# Delete temp file
file.remove(temp_output_file)
print("Temporary file has been deleted.")
