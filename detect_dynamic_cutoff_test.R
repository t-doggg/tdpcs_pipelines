#!/usr/bin/env Rscript

# --- Parse command line arguments ---
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]
cutoff_output <- args[2]

# --- Function to extract read lengths from a FASTA file ---
get_fasta_read_lengths <- function(file_path) {
  lines <- readLines(file_path)
  read_lengths <- c()
  seq_accum <- ""
  for (line in lines) {
    if (startsWith(line, ">")) {
      if (nchar(seq_accum) > 0) {
        read_lengths <- c(read_lengths, nchar(seq_accum))
        seq_accum <- ""
      }
    } else {
      seq_accum <- paste0(seq_accum, line)
    }
  }
  if (nchar(seq_accum) > 0) {
    read_lengths <- c(read_lengths, nchar(seq_accum))
  }
  return(read_lengths)
}

# --- Get read lengths from FASTA ---
lengths_all <- get_fasta_read_lengths(fasta_file)
lengths <- lengths_all[lengths_all <= 10000]

if (length(lengths) == 0) {
  stop("No read lengths ≤ 10,000 bp found. Cannot compute histogram.")
}

# --- Generate histogram ---
bin_width <- 50
max_bin <- min(10000, ceiling(max(lengths) / bin_width) * bin_width)
breaks_seq <- seq(0, max_bin, by = bin_width)

hist_data <- hist(lengths, breaks = breaks_seq, plot = FALSE)
x_vals <- hist_data$mids
cum_counts <- cumsum(hist_data$counts)

# --- Compute slope ---
slope_vals <- diff(cum_counts)
x_slope <- x_vals[-1]

# --- Get cutoff from first point with max slope ---
max_slope_idx <- which(slope_vals == max(slope_vals))[1]
dynamic_cutoff <- x_slope[max_slope_idx]

# --- Save cutoff ---
write(dynamic_cutoff, file = cutoff_output)

# --- Save plot ---
png(filename = sub("\\.txt$", ".png", cutoff_output), width = 900, height = 600)
plot(x_vals, cum_counts, type = "l", col = "blue", lwd = 2,
     main = "Cumulative Histogram with Cutoff",
     xlab = "Read Length (bp)", ylab = "Cumulative Count")
abline(v = dynamic_cutoff, col = "red", lwd = 2, lty = 2)
text(dynamic_cutoff, max(cum_counts), labels = paste0("Cutoff: ", round(dynamic_cutoff)),
     pos = 4, col = "red", offset = 1)
dev.off()
