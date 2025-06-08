# random_sample_fasta.R
# Usage: Rscript random_sample_fasta.R input_fasta output_fasta fraction

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript random_sample_fasta.R input_fasta output_fasta fraction")
}

input_fasta <- args[1]
output_fasta <- args[2]
fraction <- as.numeric(args[3])

if (fraction == 0) {
  cat("Fraction is 0. Skipping sampling.\n")
  quit(status = 0)
}

# Read FASTA
lines <- readLines(input_fasta)
headers <- grep("^>", lines)
n_reads <- length(headers)
n_sample <- ceiling(n_reads * fraction)

set.seed(42) # for reproducibility
sampled <- sort(sample(headers, n_sample))

out_lines <- c()
for (h in sampled) {
  seq_start <- h
  seq_end <- ifelse(any(headers > h), min(headers[headers > h]) - 1, length(lines))
  out_lines <- c(out_lines, lines[seq_start:seq_end])
}

writeLines(out_lines, output_fasta)