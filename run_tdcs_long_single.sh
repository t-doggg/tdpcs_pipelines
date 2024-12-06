#!/bin/bash

# Versionsnummer definieren
VERSION="0.05a"

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

# Start der Pipeline im Long Run.                  

echo ""
echo "$(datetime.done) Long Mode" 

# Prüfen ob OUTDIR bereits vorhanden ist auf dem angegebenen Pfad
if [ -d "$OUTDIR" ]; then
    echo "$(datetime.warning) Der Ordner $OUTDIR ist bereits vorhanden. Möchten Sie fortfahren? (Ja/Nein)"
    read -r response
    case $response in
        Ja|ja|J|j|Yes|yes|Y|y )
            echo "$(datetime.info) Starte Analyse-Pipeline LOOP-RUN (BETA-MODE)"
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

# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe folgender Task
if [ -d "$OUTDIR/01-CleanedReads" ]; then
	echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" 
else	    	
	#### Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation und späteren Vollständigkeit der Daten
	mkdir $OUTDIR/00-InputData

	# Ausgabe des Startzeitpunkts in die Protokolldatei
	echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" 
	echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" >> "$G_LOG_FILE"

	# Kopiere Input Fastq INFQ in OUTDIR/00-InputData
	echo "$(datetime.task) Kopiere die Input-Fastq Datei in das Ausgabeverzeichnis" >> "$G_LOG_FILE"
	cp "$INFQ" "$OUTDIR/00-InputData"
	echo "$(datetime.done) $INFQ in $OUTDIR/00-InputData erfolgreich kopiert" >> "$G_LOG_FILE"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
fi
#-------------------------------------------------------------------------------------------------------------------------------------

# Wenn der Ordner $OUTDIR/02-Classification bereits existiert, dann überspringe folgender Task
if [ -d "$OUTDIR/02-Classification" ]; then
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.warning) Alignment gegen Human Referenzgenom wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Alignment gegen Human Referenzgenom wird uebersprungen" 
else
	#### Alignment gegen Human Referenzgenom
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

	# align to ref
	mkdir -p "$OUTDIR/01-CleanedReads"
	minimap2 -ax map-ont -t "$THREADS" "$REFHUM" "$OUTDIR/00-InputData"/*fastq -o "$OUTDIR/01-CleanedReads/align_to_ref.sam"

	# Samtools erstellt aus Sam-File eine Bam-File
	samtools view -bS $OUTDIR/01-CleanedReads/align_to_ref.sam > $OUTDIR/01-CleanedReads/align_to_ref.bam

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
	echo "$(datetime.done) Alignment zu $REFHUM abgeschlossen. Ausgabedaten in $OUTDIR/01-CleanedReads" | tee -a "$ALIGN_LOG"
	echo "$(datetime.done) Alignment gegen $REFHUM abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	#-------------------------------------------------------------------------------------------------------------------------------------
	#### Extraktion der unaligned Reads aus der Bam-Datei
	# Definiere den Pfad zur Protokolldatei
	UNALIGNED_LOG="$OUTDIR/01-CleanedReads/extract_unaligned_reads.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$UNALIGNED_LOG")"

	# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
	if [ -z "$OUTDIR" ]; then
		echo "$(datetime.error) Erforderliche Variablen fehlen."
		exit 1
	fi

	# Ausgabe des Startzeitpunkts in die Protokolldatei
	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$UNALIGNED_LOG"
	echo "$(datetime.task) Starte Extrahieren nicht alignierter Reads" >> "$G_LOG_FILE"

	# Funktion zum Extrahieren nicht alignierter Reads
	extract_unaligned_reads() {
		aligned_bam="$1"
		unaligned_fastq="${aligned_bam%.bam}_unaligned.fastq"

		# Extrahiere nicht alignierte Reads und schreibe sie in eine FASTQ-Datei
		samtools view -f 4 -u "$aligned_bam" | samtools fastq - > "$unaligned_fastq"

		echo "$unaligned_fastq"
	}

	# Aufruf der Funktion zum Extrahieren nicht alignierter Reads
	aligned_bam="$OUTDIR/01-CleanedReads/align_to_ref.bam"
	unaligned_fastq=$(extract_unaligned_reads "$aligned_bam")

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen. Ausgabedaten in $OUTDIR/03-UnalignedReads" | tee -a "$UNALIGNED_LOG"
	echo "$(datetime.done) Extrahieren nicht alignierter Reads abgeschlossen" >> "$G_LOG_FILE"

	# Funktion zur Umwandlung von FASTQ in FASTA - FA wird für BLAST Suche benötigt
	convert_fastq_to_fasta() {
		unaligned_fastq="$1"
		unaligned_fasta="${unaligned_fastq%.fastq}.fa"

		# FASTQ zu FASTA konvertieren: Überspringe jede zweite Zeile nach dem Header und Quality-Werte
		awk 'NR%4==1 {print ">"substr($0,2)} NR%4==2 {print}' "$unaligned_fastq" > "$unaligned_fasta"

		echo "$unaligned_fasta"
	}

	# Aufruf der Funktion zur Umwandlung von FASTQ in FASTA
	unaligned_fasta=$(convert_fastq_to_fasta "$unaligned_fastq")

	# Ausgabe des Ergebnisses in die Protokolldatei
	echo "$(datetime.task) Umgewandelte FASTQ-Datei in FASTA: $unaligned_fasta" >> "$UNALIGNED_LOG"
	echo "$(datetime.task) Umgewandelte FASTQ-Datei in FASTA: $unaligned_fasta" >> "$G_LOG_FILE"

	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.done)  ___ ___  ____  ____   ____  ___ ___   ____  ____       ___     ___   ____     ___ "
	echo "$(datetime.done) |   |   ||    ||    \ |    ||   |   | /    ||    \     |   \   /   \ |    \   /  _]"
	echo "$(datetime.done) | _   _ | |  | |  _  | |  | | _   _ ||  o  ||  o  )    |    \ |     ||  _  | /  [_ "
	echo "$(datetime.done) |  \_/  | |  | |  |  | |  | |  \_/  ||     ||   _/     |  D  ||  O  ||  |  ||    _]"
	echo "$(datetime.done) |   |   | |  | |  |  | |  | |   |   ||  _  ||  |       |     ||     ||  |  ||   [_ "
	echo "$(datetime.done) |   |   | |  | |  |  | |  | |   |   ||  |  ||  |       |     ||     ||  |  ||     |"
	echo "$(datetime.done) |___|___||____||__|__||____||___|___||__|__||__|       |_____| \___/ |__|__||_____|"

fi
#-------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------

# Wenn der Ordner $OUTDIR/03-Results bereits existiert, dann überspringe
if [ -d "$OUTDIR/03-Results" ]; then
	echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Klassifizierung mit BLASTn wird uebersprungen" >> "$G_LOG_FILE"
	echo "$(datetime.skip) Klassifizierung mit BLASTn wird uebersprungen" 


else

	# Definiere den Pfad zur Protokolldatei für BLASTn
	BLASTN_LOG="$OUTDIR/02-Classification/blastn.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$BLASTN_LOG")"

	# Ausgabe des Startzeitpunkts in die Protokolldatei für BLASTn
	echo "$(datetime.task) Starte BLASTn"
	echo "$(datetime.task) Starte BLASTn" >> "$BLASTN_LOG"
	echo "$(datetime.task) Starte BLASTn" >> "$G_LOG_FILE"
	echo "$(datetime.task) Suche nach Fasta in $OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa" >> "$BLASTN_LOG"

	# Funktion zum Filtern der Reads nach Länge (>= 400 bp)
	filter_reads_by_length() {
	    input_fasta="$1"
	    output_fasta="$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"
	    min_length=400

	    # Filter: nur Reads >= 400 bp in neue Datei schreiben
	    awk '/^>/ {header=$0; next} length($0) >= min_length {print header; print $0}' min_length="$min_length" "$input_fasta" > "$output_fasta"

	    echo "$(datetime.done) Längenfilter angewendet. Neue Datei: $output_fasta" | tee -a "$BLASTN_LOG"
	    echo "$output_fasta"
	}

	# Funktion zur Ausführung von BLASTn aufrufen
	run_blastn() {

	    # Anwendung des Längenfilters auf die ursprüngliche Datei
	    filtered_fasta=$(filter_reads_by_length "$OUTDIR/01-CleanedReads/align_to_ref_unaligned.fa")

		filtered_fasta="$OUTDIR/01-CleanedReads/align_to_ref_unaligned.fa"

	    # Ausgabe-TSV-Dateipfad festlegen
	    ncbi_unbinned_dir="$OUTDIR/02-Classification/out_ncbi.tsv"

	    # Befehl für Blastn vorbereiten
	    blastn -task megablast -db "$DATABASE_PATH" -num_threads "$THREADS" -outfmt '6 qseqid sacc stitle ssciname nident qlen' -max_target_seqs 1 -max_hsps 1 -query "$filtered_fasta" > "$ncbi_unbinned_dir"

	    # Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für BLASTn
	    echo "$(datetime.done) BLASTn abgeschlossen." | tee -a "$BLASTN_LOG"
	    echo "$(datetime.done) BLASTn abgeschlossen" >> "$G_LOG_FILE"
	    echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

	    # Korrekter Pfad zur TSV-Datei zurückgeben
	    echo "$ncbi_unbinned_dir"
	}

	# Aufruf der Funktion
	run_blastn "$OUTDIR/01-CleanedReads/align_to_ref_unaligned_filtered.fa"

fi

#-------------------------------------------------------------------------------------------------------------------------------------

# Ausgabe von BLAST neu anordnen
# Verzeichnis, in dem sich die TSV-RAW-Datei befindet
NCBI_RAW="$OUTDIR/02-Classification/out_ncbi.tsv"

# Verzeichnis, in dem sich die TSV-MODIFIED-Datei befinden wird
NCBI_RESTRUCTURED="$OUTDIR/03-Results/res_ncbi.tsv"

# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
R_SCRIPT_DIR="/home/drk/tdcs/R"

# Rufe das R-Skript auf und übergebe den Pfad zur Ausgabedatei
Rscript "$R_SCRIPT_DIR/restructure_raw_to_genus.R" "$NCBI_RAW" "$NCBI_RESTRUCTURED"


#-------------------------------------------------------------------------------------------------------------------------------------
# Starte Shiny App zur Datenauswertung

# Shiny App zur Datenauswertung
shiny_gui() {
	# Verzeichnis, in dem sich die CSV-Datei befindet
	TSV_DIR="$OUTDIR/03-Results"

	# Pfad zur CSV-Datei
	TSV_FILE="$TSV_DIR/res_ncbi.csv"

	# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
	R_SCRIPT_DIR="/home/drk/tdcs/R"

	# Definiere den Pfad zur Protokolldatei für BLASTn
	RSCRIPT_LOG="$OUTDIR/03-Results/r_plots.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$RSCRIPT_LOG")"

	# Ausgabe des Startzeitpunkts in die Protokolldatei für BLASTn
	echo "$(datetime.task) Starte R-Skript" >> "$RSCRIPT_LOG"

	# Übergabe der Werte an das R-Skript.
	# Setze den Pfad zur Ausgabedatei als Umgebungsvariable
	export TSV_PATH="$TSV_FILE"
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

# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
sleep 20 

# Aufruf der Funktion Shiny GUI
shiny_gui




#-------------------------------------------------------------------------------------------------------------------------------------
