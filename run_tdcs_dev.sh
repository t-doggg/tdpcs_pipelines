#!/bin/bash

# Versionsnummer definieren
VERSION="0.09a-dev"

INFQ="$2"
OUTDIR="$4"
REFHUM="$6"
DATABASE_PATH="$8"
THREADS="${10}"
COVERAGE="${12}"
FRACTION="${14}"

### Pipeline Code folgt hier:

## Definition der unterschiedlichen Ausgaben:
datetime.task(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] TASK:"; }
datetime.info(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] INFO:"; }
datetime.warning(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] WARNING:"; }
datetime.error(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] ERROR:"; }
datetime.skip(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] SKIP:"; }
datetime.done(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] DONE:"; }

echo "  _______ _____  _____   _____  _____ "
echo " |__   __|  __ \|  __ \ / ____|/ ____|"
echo "    | |  | |  | | |__) | |    | (___  "
echo "    | |  | |  | |  ___/| |     \___ \ "
echo "    | |  | |__| | |    | |____ ____) |"
echo "    |_|  |_____/|_|     \_____|_____/ "
echo ""
echo "    DEV-Mode - TDCS Pipeline v$VERSION"
echo "    $INFQ"
echo "    $OUTDIR"
echo "    $REFHUM"
echo "    $DATABASE_PATH"
echo "    $THREADS Threads"
echo "    $COVERAGE Coverage"
echo "    $FRACTION FRAKTION"
echo "    DEV-Mode - TDCS Pipeline v$VERSION"
echo ""

### All functions and variables are defined 
# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe folgender Task
copy_input_to_output() {
	# Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation und späteren Vollständigkeit der Daten
	mkdir -p "$OUTDIR/00-InputData"

	# Ausgabe des Startzeitpunkts in die Protokolldatei
	echo "$(datetime.task) Starte Analyse-Pipeline DEV-Mode"
	echo "$(datetime.task) Starte Analyse-Pipeline DEV-Mode" >> "$G_LOG_FILE"

	# Kopiere Input Fastq INFQ in OUTDIR/00-InputData
	echo "$(datetime.task) Kopiere die Input-Fastq Datei in das Ausgabeverzeichnis" >> "$G_LOG_FILE"
	cp "$INFQ" "$OUTDIR/00-InputData"
	echo "$(datetime.done) $INFQ in $OUTDIR/00-InputData erfolgreich kopiert" >> "$G_LOG_FILE"
	echo "$(datetime.done) ---------------------------------------------------------------" >> "$G_LOG_FILE"
}

# Funktion für das Alignment gegen das Human-Referenzgenom und Extraktion der unalignierten Reads
align_with_minimap2() {
	# Definiere den Pfad zur Protokolldatei des Mappings
	ALIGN_LOG="$OUTDIR/01-CleanedReads/align_to_human_ref.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$ALIGN_LOG")"

	# Überprüfe, ob alle erforderlichen Variablen gesetzt wurden
	if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REFHUM" ]; then
		echo "$(datetime.error) Erforderliche Variablen fehlen."
		exit 1
	fi

	# Ausgabe des Startzeitpunkts in die Protokolldatei
	echo "$(datetime.task) Starte Minimap2-Alignment gegen $REFHUM" >> "$ALIGN_LOG"
	echo "$(datetime.task) Starte Alignment gegen $REFHUM" >> "$G_LOG_FILE"

	# Alignment mit minimap2
	mkdir -p "$OUTDIR/01-CleanedReads"
	minimap2 -ax map-ont -t "$THREADS" "$REFHUM" "$OUTDIR/00-InputData"/*fastq -o "$OUTDIR/01-CleanedReads/minimap2_hostalignment.sam"

	# Samtools erstellt aus Sam-File eine Bam-File
	samtools view -bS "$OUTDIR/01-CleanedReads/minimap2_hostalignment.sam" > "$OUTDIR/01-CleanedReads/minimap2_hostalignment.bam"

	# Sortiere und indexiere die BAM-Datei
	samtools sort "$OUTDIR/01-CleanedReads/minimap2_hostalignment.bam" -o "$OUTDIR/01-CleanedReads/minimap2_hostalignment.sorted.bam"
	samtools index "$OUTDIR/01-CleanedReads/minimap2_hostalignment.sorted.bam"

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
	echo "$(datetime.done) Alignment zu $REFHUM abgeschlossen. Ausgabedaten in $OUTDIR/01-CleanedReads" | tee -a "$ALIGN_LOG"
	echo "$(datetime.done) Alignment gegen $REFHUM abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) ---------------------------------------------------------------" >> "$G_LOG_FILE"

	#### Extraktion der unaligned Reads aus der BAM-Datei (minimap2)
	UNALIGNED_LOG="$OUTDIR/01-CleanedReads/extract_unaligned_reads.log"
	mkdir -p "$(dirname "$UNALIGNED_LOG")"

	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$UNALIGNED_LOG"
	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$G_LOG_FILE"

	# Extrahiere primäre, nicht alignierte Reads
	aligned_bam="$OUTDIR/01-CleanedReads/minimap2_hostalignment.sorted.bam"
	minimap2_unaligned_fastq="$OUTDIR/01-CleanedReads/minimap2_hostalignment_unaligned.fastq"
	samtools view -f 4 -u "$aligned_bam" | samtools fastq - > "$minimap2_unaligned_fastq"

	if [ ! -s "$minimap2_unaligned_fastq" ]; then
		echo "$(datetime.error) WARNUNG: Die Datei $minimap2_unaligned_fastq ist leer!" | tee -a "$UNALIGNED_LOG"
	fi

	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen. Ausgabedaten in $minimap2_unaligned_fastq" | tee -a "$UNALIGNED_LOG"
	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen" >> "$G_LOG_FILE"
}

align_with_bwa() {
	UNALIGNED_LOG="$OUTDIR/01-CleanedReads/extract_unaligned_reads.log"
	mkdir -p "$(dirname "$UNALIGNED_LOG")"

	minimap2_unaligned_fastq="$OUTDIR/01-CleanedReads/minimap2_hostalignment_unaligned.fastq"

	# Alignment mit bwa
	/home/diablo/miniconda3/envs/nanophase-v0.2.2/bin/bwa mem -t "$THREADS" "$REFHUM" "$minimap2_unaligned_fastq" > "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sam"

	# Samtools erstellt aus Sam-File eine Bam-File
	samtools view -bS "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sam" > "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.bam"

	# Sortiere und indexiere die BWA BAM-Datei
	samtools sort "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.bam" -o "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"
	samtools index "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"

	# Extrahiere nicht alignierte Reads aus BWA-Mapping
	bwa_aligned_bam="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"
	bwa_minimap2_unaligned_fastq="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fastq"
	samtools view -f 4 -u "$bwa_aligned_bam" | samtools fastq - > "$bwa_minimap2_unaligned_fastq"

	if [ ! -s "$bwa_minimap2_unaligned_fastq" ]; then
		echo "$(datetime.error) WARNUNG: Die Datei $bwa_minimap2_unaligned_fastq ist leer!" | tee -a "$UNALIGNED_LOG"
	fi

	# Umwandlung von FASTQ in FASTA
	unaligned_fasta="${bwa_minimap2_unaligned_fastq%.fastq}.fa"
	awk 'NR%4==1 {print ">"substr($0,2)} NR%4==2 {print}' "$bwa_minimap2_unaligned_fastq" > "$unaligned_fasta"

	echo "$(datetime.task) Umgewandelte FASTQ-Datei in FASTA: $unaligned_fasta" >> "$UNALIGNED_LOG"
	echo "$(datetime.task) Umgewandelte FASTQ-Datei in FASTA: $unaligned_fasta" >> "$G_LOG_FILE"
}

# Funktion zum dynamischen Cutoff und Längenfilter
# R-Scripte sind noch Hard-Coded
# Cut-Off sollte noch im 01-Ordner passieren
dynamic_cutoff_and_filter() {
    # Paths
    unaligned_fasta="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fa"
    filtered_fasta="$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
    CUTOFF_FILE="$OUTDIR/02-Classification/dynamic_cutoff.txt"
    R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
    BLASTN_LOG="$OUTDIR/02-Classification/blastn.log"

    # Ensure output directory exists
    mkdir -p "$(dirname "$CUTOFF_FILE")"

    Rscript "$R_SCRIPT_DIR/detect_dynamic_cutoff.R" "$unaligned_fasta" "$CUTOFF_FILE"
    CUTOFF=$(cat "$CUTOFF_FILE" | cut -d'.' -f1)
    echo "$(datetime.info) Dynamischer Cutoff ermittelt: $CUTOFF bp" | tee -a "$BLASTN_LOG"

    awk -v minlen="$CUTOFF" '
        /^>/ {header=$0; next}
        length($0) >= minlen {print header; print $0}
    ' "$unaligned_fasta" > "$filtered_fasta"

    echo "$(datetime.done) Dynamischer Längenfilter angewendet. Neue Datei: $filtered_fasta" | tee -a "$BLASTN_LOG"
}

# Funktion für Kraken2-Classification
run_kraken2() {
	KRAKEN_LOG="$OUTDIR/02-Classification/kraken2.log"
	mkdir -p "$(dirname "$KRAKEN_LOG")"
	echo "$(datetime.task) Starte Kraken2-Classification" >> "$KRAKEN_LOG"
	echo "$(datetime.task) Starte Kraken2-Classification" >> "$G_LOG_FILE"
	kraken2 --db "/home/diablo/kraken2_db" --threads "$THREADS" --report "$OUTDIR/02-Classification/kraken2_report.txt" --output "$OUTDIR/02-Classification/kraken2_output.txt" "$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
	echo "$(datetime.done) Kraken2-Classification abgeschlossen." | tee -a "$KRAKEN_LOG"
	echo "$(datetime.done) Kraken2-Classification abgeschlossen." >> "$G_LOG_FILE"
	echo "$(datetime.done) ---------------------------------------------------------------" >> "$G_LOG_FILE"
}

# Funktion für AMRFinderPlus-Assembly und -Classification
# HIER IST NOCH FRAGLICH OB ASSEMBLY WIRKLICH BENÖTIGT WIRD
# <--------------------
# <--------------------
run_amrfinder() {
    # === Pfade und Variablen ===
    filtered_reads="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fastq"
    flye_outdir="$OUTDIR/01-CleanedReads/flye_assembly"
    contigs_fasta="$flye_outdir/assembly.fasta"
    amr_output="$OUTDIR/02-Classification/amrfinder_results.tsv"
    amr_log="$OUTDIR/02-Classification/amrfinder.log"
    amr_db="/home/diablo/miniconda3/envs/nanophase-v0.2.2/share/amrfinderplus/data/2023-11-15.1"

    # Sicherstellen, dass Output-Ordner existieren
    mkdir -p "$flye_outdir"
    mkdir -p "$(dirname "$amr_output")"

    echo "$(datetime.info) Starte Assembly mit Flye im Metagenom-Modus..." | tee -a "$amr_log"

    # === Assembly mit Flye: ONT-Metagenomdaten ===
    flye --nano-raw "$filtered_reads" --meta -t "$THREADS" -o "$flye_outdir" 2>&1 | tee -a "$amr_log"

    if [ ! -f "$contigs_fasta" ]; then
        echo "$(datetime.error) ❌ Flye-Assembly fehlgeschlagen oder keine assembly.fasta gefunden!" | tee -a "$amr_log"
        return 1
    fi

    echo "$(datetime.done) Assembly abgeschlossen. Gefundene Contigs: $contigs_fasta" | tee -a "$amr_log"
    echo "$(datetime.done) Starte AMRFinderPlus auf Datei: $contigs_fasta" >> "$G_LOG_FILE"
    echo "$(datetime.info) Starte AMRFinderPlus auf Datei: $contigs_fasta" | tee -a "$amr_log"

    # === AMRFinder ausführen ===
    amrfinder -n "$contigs_fasta" \
              -o "$amr_output" \
              --database "$amr_db" \
              2>&1 | tee -a "$amr_log"

    if [ "${PIPESTATUS[0]}" -eq 0 ]; then
        echo "$(datetime.done) AMRFinderPlus erfolgreich abgeschlossen: $amr_output" | tee -a "$amr_log"
    else
        echo "$(datetime.error) ❌ AMRFinderPlus fehlgeschlagen!" | tee -a "$amr_log"
    fi
}


# Funktion zum Zusammenfassen der AMR-Ergebnisse
summarize_armfinder_results() {
    R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
	amr_output="$OUTDIR/02-Classification/amrfinder_results.tsv"
    output_summary="$OUTDIR/03-Classification/amrfinder_results_summary.csv"

	# Prüfe, ob die Datei existiert und mehr als nur die Headerzeile enthält
	if [[ ! -s "$amr_output" ]] || [[ $(awk 'END{print NR}' "$amr_output") -le 1 ]]; then
		echo "[ERROR] armfinder CSV is missing or contains no entries: $amr_output"
		return 1
	fi

    echo "🔍 Summarizing AMR hits from $amr_output ..."
    Rscript "$R_SCRIPT_DIR/summarize_armfinder_hits.R" "$amr_output" "$output_summary"
}

# Funktion zum Extrahieren der "Unknown"- und "low-resolution-classified" Reads aus Kraken2 und Umwandlung in FASTA
extract_unknown_reads_to_fasta() {
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
	KRAKEN_OUTPUT="$OUTDIR/02-Classification/kraken2_output.txt"
	KRAKEN_REPORT="$OUTDIR/02-Classification/kraken2_report.txt"
	FILTERED_FASTA="$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
	UNKNOWN_LOW_RES_IDS="$OUTDIR/02-Classification/kraken2_unknown_ids.txt"
	UNKNOWN_LOW_RES_FASTA="$OUTDIR/02-Classification/kraken2_unknown_low_res.fasta"

	echo "$(datetime.task) Starte Extraktion der Unknown und low-resolution Reads"
	echo "$(datetime.task) Starte Extraktion der Unknown und low-resolution Reads" >> "$G_LOG_FILE"
	# R-Skript aufrufen, übergibt: kraken_output, fasta_input, fasta_output
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	Rscript "$R_SCRIPT_DIR/extract_unknown_and_shallow_reads.R" "$KRAKEN_OUTPUT" "$KRAKEN_REPORT" "$FILTERED_FASTA" "$UNKNOWN_LOW_RES_FASTA" "$UNKNOWN_LOW_RES_IDS"
}

# Funktion zum Erstellen einer zufälligen Stichprobe von 10% der Reads aus der unaligned FASTA-Datei
random_sample_fasta() {
    R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
    input_fasta="$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
    output_fasta="$OUTDIR/02-Classification/ncbi_random_sample.fasta"
    mkdir -p "$(dirname "$output_fasta")"
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
    Rscript "$R_SCRIPT_DIR/random_sample_fasta.R" "$input_fasta" "$output_fasta" "$FRACTION"
}

# Funktion zum Kombinieren von ncbi_random_sample und kraken2_unknown zu einer neuen FASTA-Datei
combine_random_and_unknown_fasta() {
	random_sample_fasta="$OUTDIR/02-Classification/ncbi_random_sample.fasta"
	kraken2_unknown_fasta="$OUTDIR/02-Classification/kraken2_unknown_low_res.fasta"
	random_sample_kranken2_unknown="$OUTDIR/02-Classification/kraken2_unknown_random_sample.fasta"

	# Prüfe, ob beide Dateien existieren
	if [ -f "$random_sample_fasta" ] && [ -f "$kraken2_unknown_fasta" ]; then
		cat "$random_sample_fasta" "$kraken2_unknown_fasta" > "$random_sample_kranken2_unknown"
		echo "$(datetime.done) Kombinierte FASTA-Datei erstellt: $random_sample_kranken2_unknown" | tee -a "$G_LOG_FILE"
	else
		echo "$(datetime.warning) Eine der Eingabedateien für das Kombinieren fehlt. Vorgang wird übersprungen." | tee -a "$G_LOG_FILE"
	fi
}

# Funktion für BLASTn-Classification nur auf Unknown-Reads
run_blastn_on_unknown_five_hitter() {
	BLASTN_LOG="$OUTDIR/02-Classification/blastn.log"
	mkdir -p "$(dirname "$BLASTN_LOG")"
	echo "$(datetime.task) Starte BLASTn nur für Unknown-Reads" >> "$BLASTN_LOG"
	echo "$(datetime.task) Starte BLASTn nur für Unknown-Reads" >> "$G_LOG_FILE"

	input_ncbi_fasta="$OUTDIR/02-Classification/kraken2_unknown_random_sample.fasta"
	alt_input_ncbi_fasta="$OUTDIR/02-Classification/kraken2_unknown_low_res.fasta"
	ncbi_unknowns="$OUTDIR/02-Classification/out_ncbi.tsv"

	if [ -f "$input_ncbi_fasta" ]; then
		echo "$(datetime.info) Verwende $input_ncbi_fasta als Input für BLASTn." | tee -a "$BLASTN_LOG"
		blast_input="$input_ncbi_fasta"
	else
		echo "$(datetime.warning) $input_ncbi_fasta nicht gefunden, verwende stattdessen $alt_input_ncbi_fasta." | tee -a "$BLASTN_LOG"
		echo "$(datetime.warning) $input_ncbi_fasta nicht gefunden, verwende stattdessen $alt_input_ncbi_fasta." >> "$G_LOG_FILE"
		blast_input="$alt_input_ncbi_fasta"
	fi

	blastn -task megablast -db "$DATABASE_PATH" -num_threads "$THREADS" -outfmt '6 qseqid sacc stitle ssciname nident qlen qcovs bitscore pident' -max_target_seqs 5 -max_hsps 1 -query "$blast_input" > "$ncbi_unknowns"
	echo "$(datetime.done) BLASTn für Unknown-Reads abgeschlossen." | tee -a "$BLASTN_LOG"
	echo "$(datetime.done) BLASTn für Unknown-Reads abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) ---------------------------------------------------------------" >> "$G_LOG_FILE"
}

# Funktion zum Collapsen der NCBI BLASTn Ergebnisse mit collapse_ncbi.R
collapse_ncbi_five_hitter() {
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
	ncbi_unknowns="$OUTDIR/02-Classification/out_ncbi.tsv"
	ncbi_collapsed="$OUTDIR/02-Classification/out_ncbi_collapsed.tsv"
	ncbi_low_confidence="$OUTDIR/02-Classification/out_ncbi_low_confidence.tsv"
	COLLAPSE_LOG="$OUTDIR/02-Classification/collapse_ncbi.log"
	mkdir -p "$(dirname "$COLLAPSE_LOG")"
	echo "$(datetime.task) Starte collapse_ncbi_results.R für $ncbi_unknowns" >> "$COLLAPSE_LOG"
	# Umgebungsvariable R_LIBS_USER setzen
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	Rscript "$R_SCRIPT_DIR/collapse_ncbi_results.R" "$ncbi_unknowns" "$ncbi_collapsed" "$ncbi_low_confidence"
	echo "$(datetime.done) collapse_ncbi_results.R abgeschlossen: $ncbi_collapsed" >> "$COLLAPSE_LOG"
}

# Funktion zum Filtern der NCBI-Ergebnisse nach Coverage und Länge 300 bp
filter_ncbi_by_coverage() {
	ncbi_collapsed="$OUTDIR/02-Classification/out_ncbi_collapsed.tsv"
	f_output_file="$OUTDIR/02-Classification/f_out_ncbi.tsv"
	# $6 = qlen, $7 = qcovs
	awk -F'\t' -v cov="$COVERAGE" '$7 > cov && $6 >= 300' "$ncbi_collapsed" > "$f_output_file"
	echo "Filtered rows saved to $f_output_file"
}

restructure_ncbi_output() {
	# Umgebungsvariable R_LIBS_USER setzen
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"

	# Definiere den Pfad zur Protokolldatei für BLASTn
	RESTRUC_LOG="$OUTDIR/03-Results/restructured.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$RESTRUC_LOG")"

	# Verzeichnis, in dem sich die TSV-RAW-Datei befindet
	NCBI_RAW="$OUTDIR/02-Classification/f_out_ncbi.tsv"

	# Verzeichnis, in dem sich die TSV-MODIFIED-Datei befinden wird
	NCBI_RESTRUCTURED="$OUTDIR/02-Classification/f_res_ncbi.tsv"

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	echo "$(datetime.task) Starte Restructuring von $NCBI_RAW" >> "$RESTRUC_LOG"
	Rscript "$R_SCRIPT_DIR/restructure_raw_to_genus.R" "$NCBI_RAW" "$NCBI_RESTRUCTURED"
	echo "$(datetime.done) Restructuring abgeschlossen: $NCBI_RESTRUCTURED" >> "$RESTRUC_LOG"
}

summarize_ncbi_bacteria() {
	# Script: Bacteria summarize with synonyms - NCBI

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	# Umgebungsvariable R_LIBS_USER setzen
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	# NCBI_RESTRUCTURED muss hier noch zuerst als CSV formatiert werden (macht bisher plot_Script)
	# Daher: Konvertiere TSV zu CSV, falls noch nicht vorhanden
	NCBI_RESTRUCTURED="$OUTDIR/02-Classification/f_res_ncbi.tsv"
	NCBI_RESTRUCTURED_CSV="$OUTDIR/02-Classification/f_res_ncbi.csv"

	if [ -f "$NCBI_RESTRUCTURED" ]; then
		# Konvertiere TSV zu CSV, falls CSV nicht existiert
		if [ ! -f "$NCBI_RESTRUCTURED_CSV" ]; then
			awk 'BEGIN{OFS=","} {gsub("\t",","); print}' "$NCBI_RESTRUCTURED" > "$NCBI_RESTRUCTURED_CSV"
		fi
	fi

	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden soll
	NCBI_SUMMARIZED="$OUTDIR/03-Results/s_res_ncbi.csv"

	# Verzeichnis, in dem sich die TXT-Datei mit den Synonymen befinden wird
	BACTERIA_SYNO="/home/diablo/publ_pipeline/others/bacteria_synonyms.txt"

	# Führe das Summarize-Skript aus, falls CSV existiert
	if [ -f "$NCBI_RESTRUCTURED_CSV" ]; then
		Rscript "$R_SCRIPT_DIR/summa_bacteria_ncbi.R" "$NCBI_RESTRUCTURED_CSV" "$NCBI_SUMMARIZED" "$BACTERIA_SYNO"
	fi
}

# Funktion zum Zusammenfassen der Kraken2-Ergebnisse
summarize_kraken2_bacteria() {
	# Script: Bacteria summarize with synonyms - KRAKEN2

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	# Umgebungsvariable R_LIBS_USER setzen
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	# KRAKEN2_RAW_OUTPUT muss hier noch zuerst als CSV formatiert werden
	KRAKEN2_RAW_OUTPUT="$OUTDIR/02-Classification/kraken2_report.txt"
	KRAKEN2_RESTRUCTURED_CSV="$OUTDIR/02-Classification/kraken2_report.csv"
	if [ -f "$KRAKEN2_RAW_OUTPUT" ]; then
		# Konvertiere TSV zu CSV, falls CSV nicht existiert
		if [ ! -f "$KRAKEN2_RESTRUCTURED_CSV" ]; then
			awk 'BEGIN{OFS=","} {gsub("\t",","); print}' "$KRAKEN2_RAW_OUTPUT" > "$KRAKEN2_RESTRUCTURED_CSV"
		fi
	fi

	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden soll
	KRAKEN2_SUMMARIZED="$OUTDIR/03-Results/s_res_kraken2.csv"

	# Führe das Summarize-Skript aus, falls CSV existiert
	if [ -f "$KRAKEN2_RESTRUCTURED_CSV" ]; then
		Rscript "$R_SCRIPT_DIR/summa_bacteria_kraken2.R" "$KRAKEN2_RESTRUCTURED_CSV" "$KRAKEN2_SUMMARIZED"
	fi
}

compare_ncbi_kraken() {
	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden wird
	NCBI_SUMMARIZED="$OUTDIR/03-Results/s_res_ncbi.csv"
	
	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden wird
	KRAKEN_SUMMARIZED="$OUTDIR/03-Results/s_res_kraken2.csv"

	# Verzeichnis, in dem die compared-Datei angelegt werden soll
	COMPARED_NCBI_KRAKEN="$OUTDIR/03-Results/s_res_compared_conficence.csv"

	# Führe den Vergleich aus
	# Umgebungsvariable R_LIBS_USER setzen
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
	echo "$(datetime.task) Starte Vergleich NCBI/Kraken2" >> "$G_LOG_FILE"
	echo "$(datetime.task) Starte Vergleich NCBI/Kraken2"
	# Führe das Summarize-Skript aus, falls CSV existiert
	if [ -f "$NCBI_SUMMARIZED" ] && [ -f "$KRAKEN_SUMMARIZED" ]; then
		Rscript "$R_SCRIPT_DIR/compare_kraken_ncbi.R" "$NCBI_SUMMARIZED" "$KRAKEN_SUMMARIZED" "$COMPARED_NCBI_KRAKEN"
	fi
}

# Shiny App zur Datenauswertung
shiny_gui() {

	# Pfad zur CSV-Datei
	CSV_FILE="$OUTDIR/03-Results/s_res_compared_conficence.csv"
	AMR_FILE="$OUTDIR/03-Classification/amrfinder_results_summary.tsv"

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	export R_LIBS_USER="/home/diablo/R/x86_64-pc-linux-gnu-library/4.4/"
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	# Definiere den Pfad zur Protokolldatei für BLASTn
	RSCRIPT_LOG="$OUTDIR/03-Results/r_plots.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$RSCRIPT_LOG")"

	# Ausgabe des Startzeitpunkts in die Protokolldatei für BLASTn
	echo "$(datetime.task) Starte R-Skript" >> "$RSCRIPT_LOG"

	# Übergabe der Werte an das R-Skript.
	# Setze den Pfad zur Ausgabedatei als Umgebungsvariable
	# Hier steht noch falsch TSV_PATH, sollte aber CSV sein
	export FASTQ_PATH="$INFQ"

	# R-Skript aufrufen
	Rscript "$R_SCRIPT_DIR/plot_script.R" "$CSV_FILE" "$FASTQ_PATH" "$AMR_FILE"

	# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
	sleep 20

	# Öffne den Webbrowser mit der Shiny-App
	xdg-open "http://127.0.0.1:4010"

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für BLASTn
	echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$RSCRIPT_LOG"
	echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$G_LOG_FILE"

}

#------------------------------------------------------------------------------------------------------------------------------------- 
echo "#-------------------------------------------------------------------------------------------------------------------------------------"
echo "$(datetime.done) Pipeline in DEV-MODE BETA wird gestartet"
#------------------------------------------------------------------------------------------------------------------------------------- 

# Prüfen ob OUTDIR bereits vorhanden ist auf dem angegebenen Pfad
if [ -d "$OUTDIR" ]; then
    echo "$(datetime.warning) Der Ordner $OUTDIR ist bereits vorhanden. Möchten Sie fortfahren? (Ja/Nein)"
    read -r response
    case $response in
        Ja|ja|J|j|Yes|yes|Y|y )
            echo "$(datetime.info) Starte Analyse-Pipeline an letzter Stelle" 
            ;;
        *)
            echo "$(datetime.warning) Pipeline wird beendet."
            exit 1
            ;;
    esac
else
    # Erstellen des OUTDIR-Verzeichnisses, falls es nicht existiert
    echo "$(datetime.info) Erstelle das Verzeichnis $OUTDIR"
    mkdir -p "$OUTDIR"
fi

# Erstellen einer globalen Log-File
G_LOG_FILE="$OUTDIR/global_file.log"
	
# Erstelle das Verzeichnis, falls es nicht vorhanden ist
mkdir -p "$(dirname "$G_LOG_FILE")"	
				
# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe folgender Task
if [ -d "$OUTDIR" ]; then	
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.warning) DAS AUSGABEVERZEICHNIS IST BEREITS VORHANDEN. UEBERSPRINGE DIE SCHON BEARBEITETEN DATEIEN" >> "$G_LOG_FILE"
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"	
else				
	# Erstellen des OUTPUT Verzeichnisses
	mkdir $OUTDIR

fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Überprüfen, ob die Datei $INFQ existiert
if [ ! -f "$INFQ" ]; then
	echo "$(datetime.error) Die Datei $INFQ existiert nicht. Das Skript wird abgebrochen."
	exit 1
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation und späteren Vollständigkeit der Daten
if [ -d "$OUTDIR/01-CleanedReads" ]; then
	echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen"
else
	copy_input_to_output
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Wenn der Ordner $OUTDIR/02-Classification bereits existiert, dann überspringe folgender Task, sonst Funktion aufrufen
if [ -d "$OUTDIR/02-Classification" ]; then
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.warning) Alignments gegen Human Referenzgenom wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Alignments gegen Human Referenzgenom wird uebersprungen" 
else
	align_with_minimap2
	align_with_bwa
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Wenn der Ordner $OUTDIR/03-Results bereits existiert, dann überspringe
if [ -d "$OUTDIR/03-Results" ]; then
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Klassifizierung mit BLASTn und Kraken2 wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Klassifizierung mit BLASTn und Kraken2 wird uebersprungen"
else
	dynamic_cutoff_and_filter
	run_kraken2
	run_amrfinder # NEW
	summarize_armfinder_results # NEW
	extract_unknown_reads_to_fasta
	random_sample_fasta # NEW
	combine_random_and_unknown_fasta # NEW
	run_blastn_on_unknown_five_hitter
	collapse_ncbi_five_hitter
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Filter die NCBI-Ergebnisse nach Coverage, wenn die Datei existiert
ncbi_unbinned_dir="$OUTDIR/02-Classification/out_ncbi.tsv"
if [ -f "$ncbi_unbinned_dir" ]; then
	filter_ncbi_by_coverage
else
	echo "$(datetime.warning) BLASTn Output-File $ncbi_unbinned_dir nicht vorhanden."
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Prüfe, ob die RAW-Datei existiert, und führe die Funktion aus, falls ja
NCBI_RAW="$OUTDIR/02-Classification/f_out_ncbi.tsv"
if [ -f "$NCBI_RAW" ]; then
	restructure_ncbi_output
else
	RESTRUC_LOG="$OUTDIR/03-Results/restructured.log"
	mkdir -p "$(dirname "$RESTRUC_LOG")"
	echo "$(datetime.warning) $NCBI_RAW existiert nicht. Restructuring wird übersprungen." >> "$RESTRUC_LOG"
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Check if input file exists, then run function, else error message
NCBI_RESTRUCTURED="$OUTDIR/02-Classification/f_res_ncbi.csv"
if [ -f "$NCBI_RESTRUCTURED" ]; then
	summarize_ncbi_bacteria
else
	echo "$(datetime.error) Input file $NCBI_RESTRUCTURED does not exist. Skipping NCBI bacteria summarization."
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Check if input file exists, then run function, else error message
KRAKEN2_RESTRUCTURED="$OUTDIR/02-Classification/kraken2_output.txt"
if [ -f "$KRAKEN2_RESTRUCTURED" ]; then
	summarize_kraken2_bacteria
else
	echo "$(datetime.error) Input file $KRAKEN2_RESTRUCTURED does not exist. Skipping KRAKEN2 bacteria summarization."
fi

#------------------------------------------------------------------------------------------------------------------------------------

# Führe den Vergleich nur aus, wenn beide Dateien existieren
NCBI_SUMMARIZED="$OUTDIR/03-Results/s_res_ncbi.csv"
KRAKEN_SUMMARIZED="$OUTDIR/03-Results/s_res_kraken2.csv"
if [ -f "$NCBI_SUMMARIZED" ] && [ -f "$KRAKEN_SUMMARIZED" ]; then
	compare_ncbi_kraken
else
	echo "$(datetime.warning) Vergleich NCBI/Kraken2 wird übersprungen, da mindestens eine Datei fehlt."
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Prüfe, ob die Eingabedatei für die Shiny-App existiert, bevor shiny_gui aufgerufen wird
CSV_FILE="$OUTDIR/03-Results/s_res_compared_conficence.csv"
if [ -f "$CSV_FILE" ]; then
	# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
	sleep 20
	# Aufruf der Funktion Shiny GUI
	shiny_gui
else
	echo "$(datetime.error) Die Eingabedatei $CSV_FILE für die Shiny-App existiert nicht. Shiny-App wird nicht gestartet."
fi

#-------------------------------------------------------------------------------------------------------------------------------------
# Vergleich NCBI/Kraken2, falls die Dateien existieren
CSV_FILE="$OUTDIR/03-Results/s_res_compared_conficence.csv"
if [ -f "$CSV_FILE" ]; then
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"	
	echo "$(datetime.done) Pipeline in Long-Mode Single wurde erfolgreich abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) Pipeline in Long-Mode Single wurde erfolgreich abgeschlossen"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
else
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"	
	echo "$(datetime.done) Pipeline in Long-Mode Single fehlerhaft abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) Pipeline in Long-Mode Single fehlerhaft abgeschlossen"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
fi
#-------------------------------------------------------------------------------------------------------------------------------------

# Beende den lokalen Shiny-Server nach 10 Sekunden
sleep 10
# Finde und beende den Rscript-Prozess, der die Shiny-App auf Port 4010 startet
pkill -f "Rscript.*plot_script.R"