# === Required Libraries ===
library(dplyr)
library(readr)
library(stringr)

# === HARDCODED PATHS ===
kraken_output <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_13test_c_50_f_0/02-Classification/kraken2_output.txt"
kraken_report <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_13test_c_50_f_0/02-Classification/kraken2_report.txt"
fasta_input   <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_13test_c_50_f_0/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
output_ids    <- "/home/diablo/publ/6_TDPCS_ABGLEICH/b_13test_c_50_f_0/phylum_only_ids.txt"

# === PARAMETERS ===
rank_cutoff <- "C"  # Keep reads with terminal taxon only up to Phylum
rank_levels <- c("D", "P", "C", "O", "F", "G", "S")

# === Load Kraken2 Output ===
kraken_df <- read_tsv(kraken_output, col_names = FALSE, show_col_types = FALSE)
colnames(kraken_df) <- c("status", "read_id", "length", "mapped", "tax_path")

# === Classification Depth Function ===
extract_depth <- function(path) {
  taxid_parts <- unlist(strsplit(path, " "))
  taxids <- sapply(strsplit(taxid_parts, ":"), `[`, 1)
  sum(taxids != "0")
}

# === Extract Last (Terminal) TaxID Function ===
extract_terminal_taxid <- function(path) {
  taxid_parts <- str_split(path, " ", simplify = TRUE)
  taxids <- str_replace(taxid_parts, ":.*", "")  # remove after colon
  taxids <- str_trim(taxids)                     # remove leading/trailing spaces
  taxids <- taxids[taxids != "0"]                # filter non-zero
  if (length(taxids) == 0) return(NA)
  return(as.character(tail(taxids, 1)))
}

# === Apply to Kraken Output ===
kraken_df <- kraken_df %>%
  mutate(
    classification_depth = sapply(tax_path, extract_depth),
    terminal_taxid = sapply(tax_path, extract_terminal_taxid)
  )

# === Load Kraken Report (TaxID to Rank Mapping) ===
report_df <- read_tsv(kraken_report, col_names = FALSE, show_col_types = FALSE, comment = "#")
colnames(report_df) <- c("pct", "reads", "clade_reads", "rank", "taxid", "name")
taxid_to_rank <- report_df %>%
  select(taxid, rank) %>%
  mutate(taxid = as.character(taxid))

# === Join Taxonomic Rank Info ===
kraken_df <- kraken_df %>%
  left_join(taxid_to_rank, by = c("terminal_taxid" = "taxid")) %>%
  rename(terminal_rank = rank)

# === Define Allowed Ranks (up to specified level) ===
allowed_ranks <- rank_levels[1:which(rank_levels == rank_cutoff)]

# === Adaptive Sampling Filter ===
selected_reads <- kraken_df %>%
  filter(status == "U" | terminal_rank %in% allowed_ranks) %>%
  pull(read_id)

# === Output Selected Read IDs ===
write_lines(selected_reads, output_ids)

cat("✅ Adaptive sampling completed.\n")
cat("→ Rank threshold: <", rank_cutoff, ">\n")
cat("→ Total reads selected:", length(selected_reads), "\n")
cat("→ Output saved to:", output_ids, "\n")
