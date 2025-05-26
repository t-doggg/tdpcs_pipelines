#!/bin/bash

# Versionsnummer definieren
VERSION="0.09a-dev"

INFQ="$2"
OUTDIR="$4"
REFHUM="$6"
DATABASE_PATH="$8"
THREADS="${10}"

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

echo "Version: $VERSION"
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
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
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

	# hostdepletion with minimap2
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
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

	#### Extraktion der unaligned Reads aus der BAM-Datei
	UNALIGNED_LOG="$OUTDIR/01-CleanedReads/extract_unaligned_reads.log"
	mkdir -p "$(dirname "$UNALIGNED_LOG")"

	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$UNALIGNED_LOG"
	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$G_LOG_FILE"

	# Extrahiere primäre, nicht alignierte Reads
	aligned_bam="$OUTDIR/01-CleanedReads/minimap2_hostalignment.sorted.bam"
	unaligned_fastq="$OUTDIR/01-CleanedReads/minimap2_hostalignment_unaligned.fastq"
	samtools view -f 4 -u "$aligned_bam" | samtools fastq - > "$unaligned_fastq"

	if [ ! -s "$unaligned_fastq" ]; then
		echo "$(datetime.error) WARNUNG: Die Datei $unaligned_fastq ist leer!" | tee -a "$UNALIGNED_LOG"
	fi

	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen. Ausgabedaten in $unaligned_fastq" | tee -a "$UNALIGNED_LOG"
	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen" >> "$G_LOG_FILE"
}

align_with_bwa() {
	UNALIGNED_LOG="$OUTDIR/01-CleanedReads/extract_unaligned_reads.log"
	mkdir -p "$(dirname "$UNALIGNED_LOG")"

	unaligned_fastq="$OUTDIR/01-CleanedReads/minimap2_hostalignment_unaligned.fastq"

	# hostdepletion with bwa
	#/home/diablo/miniconda3/envs/nanophase-v0.2.2/bin/bwa index "$REFHUM"
	/home/diablo/miniconda3/envs/nanophase-v0.2.2/bin/bwa mem -t "$THREADS" "$REFHUM" "$unaligned_fastq" > "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sam"

	# Samtools erstellt aus Sam-File eine Bam-File
	samtools view -bS "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sam" > "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.bam"

	# Sortiere und indexiere die BWA BAM-Datei
	samtools sort "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.bam" -o "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"
	samtools index "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"

	# Extrahiere nicht alignierte Reads aus BWA-Mapping
	bwa_aligned_bam="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment.sorted.bam"
	bwa_unaligned_fastq="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fastq"
	samtools view -f 4 -u "$bwa_aligned_bam" | samtools fastq - > "$bwa_unaligned_fastq"

	if [ ! -s "$bwa_unaligned_fastq" ]; then
		echo "$(datetime.error) WARNUNG: Die Datei $bwa_unaligned_fastq ist leer!" | tee -a "$UNALIGNED_LOG"
	fi

	# Umwandlung von FASTQ in FASTA
	unaligned_fasta="${bwa_unaligned_fastq%.fastq}.fa"
	awk 'NR%4==1 {print ">"substr($0,2)} NR%4==2 {print}' "$bwa_unaligned_fastq" > "$unaligned_fasta"

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
	kraken2 --db "/home/diablo/kraken2_db" --threads "$THREADS" --report "$OUTDIR/02-Classification/kraken2_report.txt" --output "$OUTDIR/02-Classification/kraken2_output.txt" "$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fastq"
	echo "$(datetime.done) Kraken2-Classification abgeschlossen." | tee -a "$KRAKEN_LOG"
	echo "$(datetime.done) Kraken2-Classification abgeschlossen." >> "$G_LOG_FILE"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
}

# Funktion zum Extrahieren der "Unknown"-Reads aus Kraken2 und Umwandlung in FASTA
extract_unknown_reads_to_fasta() {
	KRAKEN_OUTPUT="$OUTDIR/02-Classification/kraken2_output.txt"
	INPUT_FASTQ="$OUTDIR/01-CleanedReads/bwa_minimap2_host_alignment_unaligned.fastq"
	UNKNOWN_IDS="$OUTDIR/02-Classification/kraken2_unknown_ids.txt"
	UNKNOWN_FASTQ="$OUTDIR/02-Classification/kraken2_unknown.fastq"
	UNKNOWN_FASTA="$OUTDIR/02-Classification/kraken2_unknown.fasta"

	# Extrahiere IDs der "U"-Reads (unclassified/Unknown) aus Kraken2
	awk '$1=="U"{print $2}' "$KRAKEN_OUTPUT" > "$UNKNOWN_IDS"

	# Extrahiere die entsprechenden Reads aus FASTQ
	awk 'NR==FNR{ids[$1]; next} /^@/ {h=substr($0,2); p=(h in ids)} p' "$UNKNOWN_IDS" "$INPUT_FASTQ" > "$UNKNOWN_FASTQ"

	# Wandle FASTQ in FASTA um
	awk 'NR%4==1{printf(">%s\n",substr($0,2))} NR%4==2{print}' "$UNKNOWN_FASTQ" > "$UNKNOWN_FASTA"
}

# Funktion für BLASTn-Classification nur auf Unknown-Reads
run_blastn_on_unknown_five_hitter() {
	BLASTN_LOG="$OUTDIR/02-Classification/blastn.log"
	mkdir -p "$(dirname "$BLASTN_LOG")"
	echo "$(datetime.task) Starte BLASTn nur für Unknown-Reads" >> "$BLASTN_LOG"
	echo "$(datetime.task) Starte BLASTn nur für Unknown-Reads" >> "$G_LOG_FILE"
	unknown_fasta="$OUTDIR/02-Classification/kraken2_unknown.fasta"
	ncbi_unbinned_dir="$OUTDIR/02-Classification/out_ncbi.tsv"
	blastn -task megablast -db "$DATABASE_PATH" -num_threads "$THREADS" -outfmt '6 qseqid sacc stitle ssciname nident qlen qcovs bitscore pident' -max_target_seqs 5 -max_hsps 1 -query "$unknown_fasta" > "$ncbi_unbinned_dir"
	echo "$(datetime.done) BLASTn für Unknown-Reads abgeschlossen." | tee -a "$BLASTN_LOG"
	echo "$(datetime.done) BLASTn für Unknown-Reads abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
}

# Funktion zum Collapsen der NCBI BLASTn Ergebnisse mit collapse_ncbi.R
collapse_ncbi_five_hitter() {
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"
	ncbi_unbinned_dir="$OUTDIR/02-Classification/out_ncbi.tsv"
	ncbi_collapsed_dir="$OUTDIR/02-Classification/out_ncbi_collapsed.tsv"
	COLLAPSE_LOG="$OUTDIR/02-Classification/collapse_ncbi.log"
	mkdir -p "$(dirname "$COLLAPSE_LOG")"
	echo "$(datetime.task) Starte collapse_ncbi.R für $ncbi_unbinned_dir" >> "$COLLAPSE_LOG"
	Rscript "$R_SCRIPT_DIR/collapse_ncbi.R" "$ncbi_unbinned_dir" "$ncbi_collapsed_dir"
	echo "$(datetime.done) collapse_ncbi.R abgeschlossen: $ncbi_collapsed_dir" >> "$COLLAPSE_LOG"
}

# Funktion zum Filtern der NCBI-Ergebnisse nach Coverage
filter_ncbi_by_coverage() {
	ncbi_collapsed_dir="$OUTDIR/02-Classification/out_ncbi.tsv"
	f_output_file="$OUTDIR/02-Classification/f_out_ncbi.tsv"
	awk -F'\t' '$7 > 75' "$ncbi_collapsed_dir" > "$f_output_file"
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
	NCBI_RESTRUCTURED="$OUTDIR/03-Results/f_res_ncbi.tsv"

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	echo "$(datetime.task) Starte Restructuring von $NCBI_RAW" >> "$RESTRUC_LOG"
	Rscript "$R_SCRIPT_DIR/restructure_raw_to_genus.R" "$NCBI_RAW" "$NCBI_RESTRUCTURED"
	echo "$(datetime.done) Restructuring abgeschlossen: $NCBI_RESTRUCTURED" >> "$RESTRUC_LOG"
}

summarize_ncbi_bacteria() {
	# Script: Bacteria summarize with synonyms - NCBI

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	# NCBI_RESTRUCTURED muss hier noch zuerst als CSV formatiert werden (macht bisher plot_Script)
	# Daher: Konvertiere TSV zu CSV, falls noch nicht vorhanden
	NCBI_RESTRUCTURED="$OUTDIR/03-Results/f_res_ncbi.tsv"
	NCBI_RESTRUCTURED_CSV="$OUTDIR/03-Results/f_res_ncbi.csv"

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
		Rscript "$R_SCRIPT_DIR/summa_bacteria.R" "$NCBI_RESTRUCTURED_CSV" "$NCBI_SUMMARIZED" "$BACTERIA_SYNO"
	fi
}

summarize_kraken2_bacteria() {
	# Script: Bacteria summarize with synonyms - KRAKEN2

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	R_SCRIPT_DIR="/home/diablo/publ_pipeline/R"

	# KRAKEN2_RESTRUCTURED muss hier noch zuerst als CSV formatiert werden
	KRAKEN2_RESTRUCTURED="$OUTDIR/02-Classification/kraken2_output.txt"
	KRAKEN2_RESTRUCTURED_CSV="$OUTDIR/02-Classification/kraken2_output.csv"
	if [ -f "$KRAKEN2_RESTRUCTURED" ]; then
		# Konvertiere TSV zu CSV, falls CSV nicht existiert
		if [ ! -f "$KRAKEN2_RESTRUCTURED_CSV" ]; then
			awk 'BEGIN{OFS=","} {gsub("\t",","); print}' "$KRAKEN2_RESTRUCTURED" > "$KRAKEN2_RESTRUCTURED_CSV"
		fi
	fi

	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden soll
	KRAKEN2_SUMMARIZED="$OUTDIR/03-Results/s_res_kraken2.csv"

	# Verzeichnis, in dem sich die TXT-Datei mit den Synonymen befinden wird
	BACTERIA_SYNO="/home/diablo/publ_pipeline/others/bacteria_synonyms.txt"

	# Führe das Summarize-Skript aus, falls CSV existiert
	if [ -f "$KRAKEN2_RESTRUCTURED_CSV" ]; then
		Rscript "$R_SCRIPT_DIR/summa_bacteria.R" "$KRAKEN2_RESTRUCTURED_CSV" "$KRAKEN2_SUMMARIZED" "$BACTERIA_SYNO"
	fi
}

compare_ncbi_kraken() {
	# Verzeichnis, in dem sich die csv-summarized-filtered-Datei befinden wird
	KRAKEN_SUMMARIZED="$OUTDIR/03-Results/s_res_kraken.csv"

	# Verzeichnis, in dem die compared-Datei angelegt werden soll
	COMPARED_NCBI_KRAKEN="$OUTDIR/03-Results/s_res_compared_conficence.csv"

	# Führe den Vergleich aus
	Rscript "$R_SCRIPT_DIR/compare_kraken_ncbi.R" "$NCBI_SUMMARIZED" "$KRAKEN_SUMMARIZED" "$COMPARED_NCBI_KRAKEN"
}

# Shiny App zur Datenauswertung
shiny_gui() {

	# Pfad zur CSV-Datei
	CSV_FILE="$OUTDIR/03-Results/f_res_ncbi.tsv"

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
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
	export TSV_PATH="$CSV_FILE"
	export FASTQ_PATH="$INFQ"

	# R-Skript aufrufen
	Rscript "$R_SCRIPT_DIR/plot_script.R" "$TSV_PATH" "$FASTQ_PATH" &

	# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
	sleep 2

	# Öffne den Webbrowser mit der Shiny-App
	xdg-open "http://127.0.0.1:4010"

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für BLASTn
	echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$RSCRIPT_LOG"
	echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$G_LOG_FILE"

}

#------------------------------------------------------------------------------------------------------------------------------------- 
echo "#-------------------------------------------------------------------------------------------------------------------------------------"
echo "$(datetime.done) Pipeline in Long-Mode Single wird gestartet"
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
	echo "$(datetime.warning) Alignment gegen Human Referenzgenom wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Alignment gegen Human Referenzgenom wird uebersprungen" 
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
	extract_unknown_reads_to_fasta
	run_blastn_on_unknown_five_hitter
	collapse_ncbi_five_hitter
fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Filter die NCBI-Ergebnisse nach Coverage, wenn die Datei existiert
ncbi_unbinned_dir="$OUTDIR/02-Classification/out_ncbi.tsv"
if [ -f "$ncbi_unbinned_dir" ]; then
	filter_ncbi_by_coverage
else
	echo "$(datetime.warning) Input file $ncbi_unbinned_dir does not exist. Skipping filtering step."
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
NCBI_RESTRUCTURED="$OUTDIR/03-Results/f_res_ncbi.tsv"
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
if [ -f "$NCBI_SUMMARIZED" ] && [ -f "$OUTDIR/03-Results/s_res_kraken.csv" ]; then
	compare_ncbi_kraken
else
	echo "$(datetime.warning) Vergleich NCBI/Kraken2 wird übersprungen, da mindestens eine Datei fehlt." >> "$G_LOG_FILE"
fi


#-------------------------------------------------------------------------------------------------------------------------------------

# Prüfe, ob die Eingabedatei für die Shiny-App existiert, bevor shiny_gui aufgerufen wird
CSV_FILE="$OUTDIR/03-Results/f_res_ncbi.tsv"
if [ -f "$CSV_FILE" ]; then
	# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
	sleep 20
	# Aufruf der Funktion Shiny GUI
	shiny_gui
else
	echo "$(datetime.error) Die Eingabedatei $CSV_FILE für die Shiny-App existiert nicht. Shiny-App wird nicht gestartet."
fi

#-------------------------------------------------------------------------------------------------------------------------------------

echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"	
echo "$(datetime.done) Pipeline in Long-Mode Single wurde erfolgreich abgeschlossen" >> "$G_LOG_FILE"
echo "$(datetime.done) Pipeline in Long-Mode Single wurde erfolgreich abgeschlossen"
echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

#-------------------------------------------------------------------------------------------------------------------------------------

# Beende den lokalen Shiny-Server nach 10 Sekunden
sleep 10
# Finde und beende den Rscript-Prozess, der die Shiny-App auf Port 4010 startet
pkill -f "Rscript.*plot_script.R"
# Beende den lokalen Shiny-Server
