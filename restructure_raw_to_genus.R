## Anmerkung: TSV wird ab 0.04b jetzt im restructure_raw_to_genus.R-Skript zur CSV umgewandelt.
## Zuvor war es im Shell Code enthalten aber über R lässt sich es einfacher umsetzen
## Umwandlung ist ganz unten, temp file wird gelöscht aber TSV-File Output
## bleibt bestehen!

# Notwendige Bibliotheken laden
library(dplyr)
library(stringr)

# Holen der Kommandozeilenargumente
args <- commandArgs(trailingOnly = TRUE)

# Überprüfen, ob die erforderlichen Argumente bereitgestellt wurden
if (length(args) < 2) {
  stop("Bitte sowohl den Pfad zur Eingabedatei als auch den Pfad zur Ausgabedatei als Argumente angeben.")
}

# Zuweisung der Pfade für die Eingabe- und Ausgabedatei aus den Argumenten
input_file <- args[1]
output_file <- args[2]

# Debugging: Überprüfen, ob die Pfade korrekt übergeben wurden
print(paste("CSV Path:", input_file))
print(paste("FASTQ Path:", output_file))

# Generieren eines temporären Dateipfads
temp_output_file <- tempfile(fileext = ".tsv")

# Funktion zur Verarbeitung der Bakterienzählungen aus der Eingabedatei
process_bacteria_counts_file <- function(input_file, temp_output_file) {
  
  # Schritt 1: Einlesen der TSV-Datei ohne Header
  data <- read.delim(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Überprüfen, ob die Eingabedatei mindestens drei Spalten hat
  if (ncol(data) < 3) {
    stop("Input file hat nicht mindestens drei Spalten.")
  }
  
  # Schritt 2: Extrahieren der dritten Spalte (vermutlich Bakteriennamen)
  bacteria_column <- data[[3]]  # Extract the third column directly
  
  # Schritt 3: Entfernen von "[" und "]" aus den Bakteriennamen - macht BLAST machnal
  # weil die Namen si in der nt_prok hinterlegt sind
  bacteria_column_clean <- str_replace_all(bacteria_column, "\\[|\\]", "")
  
  # Schritt 4: Extrahieren der ersten zwei Wörter aus jedem Eintrag der bereinigten Bakterien-Spalte
  bacteria_names <- sapply(strsplit(bacteria_column_clean, " "), function(x) paste(x[1:2], collapse = " "))
  
  # Schritt 5: Zählen der Vorkommen jedes einzigartigen Bakteriennamens
  bacteria_counts <- bacteria_names %>%
    table() %>%
    as.data.frame()
  
  # Schritt 6: Umbenennen der Spalten für bessere Verständlichkeit
  colnames(bacteria_counts) <- c("Bacteria", "Count")
  
  # Schritt 7: Schreiben des Ergebnisses in die temporäre TSV-Datei
  write.table(bacteria_counts, file = temp_output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  print(paste("Bacteria counts have been saved to the temporary file:", temp_output_file))
}

# Funktion zur Verarbeitung der Gattungszählungen aus der temporären Datei
process_genus_counts <- function(temp_output_file, output_file) {
  # Schritt 1: Einlesen der TSV-Datei mit Bakteriennamen und Zählungen
  data <- read.delim(temp_output_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Schritt 2: Extrahieren des ersten Wortes (Genus) aus der "Bakterium"-Spalte
  genus_names <- sapply(strsplit(data$Bacteria, " "), function(x) x[1])
  
  # Schritt 3: Kombinieren von Genusnamen und der entsprechenden Zählungen
  genus_data <- data.frame(Genus = genus_names, Count = data$Count)
  
  # Step 4: Sum der Zählungen für jede einzigartige Genusnamen
  genus_counts <- genus_data %>%
    group_by(Genus) %>%
    summarise(Total_Count = sum(Count))  # Summing counts for the same genus
  
  # Schritt 5: Schreiben des Ergebnisses in die finale Ausgabedatei
  write.table(genus_counts, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  print(paste("Genus counts have been saved to", output_file))
}

# Funktion zur Konvertierung der finalen TSV-Datei in eine CSV-Datei mit modifizierten Spaltennamen
convert_tsv_to_csv <- function(tsv_file, csv_file) {
  # Einlesen der TSV-Datei
  data <- read.delim(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Umbenennen der Spalten
  colnames(data) <- c("Stamm", "Count")
  
  # Schreiben der Daten in eine CSV-Datei
  write.csv(data, file = csv_file, row.names = FALSE, quote = FALSE)
  
  print(paste("CSV file has been saved to", csv_file))
}


# Run das R-Skript
process_bacteria_counts_file(input_file, temp_output_file)   # Step 1: Process bacteria counts to temp file
process_genus_counts(temp_output_file, output_file)          # Step 2: Process genus counts to final output file

# Generiere CSV
# ab jetzt im R-Skript integriert weil hier über Shell doof
csv_output_file <- sub("\\.tsv$", ".csv", output_file)  # Replace ".tsv" with ".csv"
convert_tsv_to_csv(output_file, csv_output_file)        # Convert the TSV to CSV

# Löschen der temporären Datei
file.remove(temp_output_file)
print("Temporary file has been deleted.")
