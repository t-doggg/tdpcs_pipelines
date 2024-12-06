#!/bin/bash

# Versionsnummer definieren
VERSION="0.05a"

INFQ="$2"
OUTDIR="$4"
REFHUM="$6"
DATABASE_PATH="$8"
THREADS="${10}"

### Pipeline Code

## Definition der Ausgaben:
datetime.task(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] TASK:"; }
datetime.info(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] INFO:"; }
datetime.warning(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] WARNING:"; }
datetime.error(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] ERROR:"; }
datetime.skip(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] SKIP:"; }
datetime.done(){ Date=$(date "+%Y-%m-%d %H:%M:%S"); echo "[$Date] DONE:"; }

# Verarbeite die Eingangsdaten
# Zum Beispiel:
#echo "INFQ: 			$INFQ"
#echo "OUTDIR: 		$OUTDIR"
#echo "REFHUM: 		$REFHUM"
#echo "DATABASE_PATH: 		$DATABASE_PATH"
#echo "THREADS: 		$THREADS"


echo "  _____ _______ _____  _____     __ "
echo " |  __ \__   __/ ____|/ ____|   /_ |"
echo " | |  | | | | | |    | (_____   _| |"
echo " | |  | | | | | |     \___ \ \ / / |"
echo " | |__| | | | | |____ ____) \ V /| |"
echo " |_____/  |_|  \_____|_____/ \_/ |_|"


# Start der Pipeline im Single Run.                  

echo ""
echo "$(datetime.done) High Throughput Mode" 

        echo "$(datetime.info) Starte Analyse-Pipeline SINGLE-RUN (BETA-MODE) bei letzter Analyse-Position"
		 # Erstellen einer globalen Log-File
		G_LOG_FILE="$OUTDIR/global_file.log"
			
		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$G_LOG_FILE")"	
            	     
		# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe
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

		# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe
		if [ -d "$OUTDIR/01-CleanedReads" ]; then
			echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" >> "$G_LOG_FILE"
			echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" 
		else	    	
			#### Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation
			mkdir $OUTDIR/00-InputData

			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" 
			echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" >> "$G_LOG_FILE"

			# Kopiere Input Fastq in OutDIR
			echo "$(datetime.task) Kopiere die Input-Fastq Datei in das Ausgabeverzeichnis" >> "$G_LOG_FILE"
			cp "$INFQ" "$OUTDIR/00-InputData"
			echo "$(datetime.done) $INFQ in $OUTDIR/00-InputData erfolgreich kopiert" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		fi
		#-------------------------------------------------------------------------------------------------------------------------------------

		# Wenn der Ordner $OUTDIR/02-De_Novo_Assembly bereits existiert, dann überspringe
		if [ -d "$OUTDIR/02-De_Novo_Assembly" ]; then
			echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
			echo "$(datetime.warning) Alignment gegen Human Referenzgenom wird uebersprungen" >> "$G_LOG_FILE"
			echo "$(datetime.skip) Alignment gegen Human Referenzgenom wird uebersprungen" 
		else
			#### Alignment gegen Human Referenzgenom
			# Definiere den Pfad zur Protokolldatei
			ALIGN_LOG="$OUTDIR/01-CleanedReads/align_to_human_ref.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$ALIGN_LOG")"

			# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
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

		# Wenn der Ordner $OUTDIR/03-Bins bereits existiert, dann überspringe
		if [ -d "$OUTDIR/03-Bins" ]; then
			echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
			echo "$(datetime.skip) DeNovo Assembly wird uebersprungen" >> "$G_LOG_FILE"
				echo "$(datetime.skip) DeNovo Assembly wird uebersprungen" 
		else
				### DeNovo Assembly
				# Definiere den Pfad zur Protokolldatei
			FLYE_LOG="$OUTDIR/02-De_Novo_Assembly/de_novo_assembly.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$ALIGN_LOG")"

			# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
			if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REFHUM" ]; then
				echo "$(datetime.error) Erforderliche Variablen fehlen."
				exit 1
			fi
			

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$FLYE_LOG")"

			# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
			if [ -z "$OUTDIR" ] || [ -z "$THREADS" ]; then
				echo "$(datetime.error) Erforderliche Variablen fehlen."
				exit 1
			fi

			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte De-Novo-Assembly" >> "$FLYE_LOG"
			echo "$(datetime.task) Starte De-Novo-Assembly" >> "$G_LOG_FILE"

			# Assembly mit Flye durchführen
			mkdir -p "$OUTDIR/02-De_Novo_Assembly"
			flye --meta --nano-hq "$OUTDIR/01-CleanedReads/align_to_ref_unaligned.fastq" -t "$THREADS" -o "$OUTDIR/02-De_Novo_Assembly/Flye"

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
			echo "$(datetime.done) De-Novo-Assembly abgeschlossen. Ausgabedaten in $OUTDIR/02-De_Novo_Assembly/Flye" | tee -a "$FLYE_LOG"
			echo "$(datetime.done) De-Novo-Assembly abgeschlossen. Ausgabedaten in $OUTDIR/02-De_Novo_Assembly/Flye" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
			echo "$(datetime.done)  _____  _      __ __    ___      ___     ___   ____     ___ "
			echo "$(datetime.done) |     || |    |  |  |  /  _]    |   \   /   \ |    \   /  _]"
			echo "$(datetime.done) |   __|| |    |  |  | /  [_     |    \ |     ||  _  | /  [_ "
			echo "$(datetime.done) |  |_  | |___ |  ~  ||    _]    |  D  ||  O  ||  |  ||    _]"
			echo "$(datetime.done) |   _] |     ||___, ||   [_     |     ||     ||  |  ||   [_ "
			echo "$(datetime.done) |  |   |     ||     ||     |    |     ||     ||  |  ||     |"
			echo "$(datetime.done) |__|   |_____||____/ |_____|    |_____| \___/ |__|__||_____|"
		fi
		#-------------------------------------------------------------------------------------------------------------------------------------

		# Wenn der Ordner $OUTDIR/04-Classify bereits existiert, dann überspringe
		if [ -d "$OUTDIR/04-Classify" ]; then
			echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
			echo "$(datetime.skip) Binning mit Metawrap & Binning-Verfeinerung mit Metawrap wird uebersprungen" >> "$G_LOG_FILE"
				echo "$(datetime.skip) Binning mit Metawrap & Binning-Verfeinerung mit Metawrap wird uebersprungen" 
		else
				### Binning mit Metawrap
			# Definiere den Pfad zur Protokolldatei
			BINNING_LOG="$OUTDIR/03-Bins/binning.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$BINNING_LOG")"

			# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
			if [ -z "$OUTDIR" ] || [ -z "$THREADS" ]; then
				echo "$(datetime.error) Erforderliche Variablen fehlen."
				exit 1
			fi

			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte Binning mit Metawrap" >> "$BINNING_LOG"
			echo "$(datetime.task) Starte Binning mit Metawrap" >> "$G_LOG_FILE"

			# Metawrap Binning
			mkdir -p "$OUTDIR/03-Bins"
			mkdir -p "$OUTDIR/03-Bins/Binning"

			# Binning mit Metawrap
			metawrap binning -t "$THREADS" -a "$OUTDIR/02-De_Novo_Assembly/Flye/assembly.fasta" -o "$OUTDIR/03-Bins/Binning" --metabat2 --maxbin2 --single-end "$OUTDIR/00-InputData/"*fastq

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
			echo "$(datetime.done) Binning mit Metawrap abgeschlossen. Ausgabedaten in $OUTDIR/03-Bins/Binning" | tee -a "$BINNING_LOG"
			echo "$(datetime.done) Binning mit Metawrap abgeschlossen" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

			
			
			
			### Binning-Verfeinerung mit Metawrap

			# Definiere den Pfad zur Protokolldatei
			BIN_REFINEMENT_LOG="$OUTDIR/03-Bins/bin_refinement.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$BIN_REFINEMENT_LOG")"

			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte Binning-Verfeinerung mit Metawrap" >> "$BIN_REFINEMENT_LOG"
			echo "$(datetime.task) Starte Binning-Verfeinerung mit Metawrap" >> "$G_LOG_FILE"

			# Metawrap Binning-Verfeinerung
			mkdir -p "$OUTDIR/03-Bins/Bin_Refinement"

			# Binning-Verfeinerung mit Metawrap
			metawrap bin_refinement -t "$THREADS" -A "$OUTDIR/03-Bins/Binning/metabat2_bins" -B "$OUTDIR/03-Bins/Binning/maxbin2_bins" -o "$OUTDIR/03-Bins/Bin_Refinement" --skip-checkm

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
			echo "$(datetime.done) Binning-Verfeinerung mit Metawrap abgeschlossen. Ausgabedaten in $OUTDIR/03-Bins/Bin_Refinement" | tee -a "$BIN_REFINEMENT_LOG"
			echo "$(datetime.done) Binning-Verfeinerung mit Metawrap abgeschlossen" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

			# Definiere den Pfad zur Protokolldatei für die Split-Funktion
			SPLIT_LOG="$OUTDIR/03-Bins/split_fasta.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$SPLIT_LOG")"

			# Ausgabe des Startzeitpunkts in die Protokolldatei für die Split-Funktion
			echo "$(datetime.task) Starte Split-Funktion" >> "$SPLIT_LOG"
			echo "$(datetime.task) Starte Split-Funktion" >> "$G_LOG_FILE"

			# Funktion zur Aufteilung der unbinned.fa-Datei aufrufen, wenn vorhanden
			split_fasta() {
				refined_bins_dir="$1"
				unbinned_fasta_path=$(find "$refined_bins_dir" -type f -name "bin.unbinned.fa" | head -n 1)

				if [ -z "$unbinned_fasta_path" ]; then
				echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen."
				echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen." >> "$SPLIT_LOG"
				echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen." >> "$G_LOG_FILE"
				return 1
				fi

				contig_count=0
				while IFS= read -r line; do
				if [[ $line == ">"* ]]; then
					# Neue Contig gefunden
					contig_count=$((contig_count+1))
					contig_name="${line:1}"
					printf '%s\n' "$line" > "${unbinned_fasta_path%.fa}.contig${contig_count}.fa"
				else
					# Sequenz der aktuellen Contig speichern
					printf '%s\n' "$line" >> "${unbinned_fasta_path%.fa}.contig${contig_count}.fa"
				fi
				done < "$unbinned_fasta_path"

				echo "$(datetime.done) Die Datei $unbinned_fasta_path wurde erfolgreich aufgespalten."
			}

			# Aufruf der Funktion
			split_fasta "$OUTDIR/03-Bins/Bin_Refinement/work_files/binsM"

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für die Split-Funktion
			echo "$(datetime.done) Split-Funktion abgeschlossen." | tee -a "$SPLIT_LOG"
			echo "$(datetime.done) Split-Funktion abgeschlossen" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------" >> "$G_LOG_FILE"
			echo "$(datetime.done)  ___ ___    ___ ______   ____  __    __  ____    ____  ____       ___     ___   ____     ___ "
			echo "$(datetime.done) |   |   |  /  _]      | /    ||  |__|  ||    \  /    ||    \     |   \   /   \ |    \   /  _]"
			echo "$(datetime.done) | _   _ | /  [_|      ||  o  ||  |  |  ||  D  )|  o  ||  o  )    |    \ |     ||  _  | /  [_ "
			echo "$(datetime.done) |  \_/  ||    _]_|  |_||     ||  |  |  ||    / |     ||   _/     |  D  ||  O  ||  |  ||    _]"
			echo "$(datetime.done) |   |   ||   [_  |  |  |  _  ||        ||    \ |  _  ||  |       |     ||     ||  |  ||   [_ "
			echo "$(datetime.done) |   |   ||     | |  |  |  |  | \      / |  .  \|  |  ||  |       |     ||     ||  |  ||     |"
			echo "$(datetime.done) |___|___||_____| |__|  |__|__|  \_/\_/  |__|\_||__|__||__|       |_____| \___/ |__|__||_____|"
		fi
		#-------------------------------------------------------------------------------------------------------------------------------------

		### Klassifizierung mit GTDB-Tk
		# Prüfe, ob der Ausgabeordner bereits existiert
		if [ -d "$OUTDIR/04-Classify" ]; then
			# Wenn der Ausgabeordner existiert, entferne ihn und vermerke es im Terminal und Log
			echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/04-Classify" >> "$G_LOG_FILE"
			echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/04-Classify"
			echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/05-Results" >> "$G_LOG_FILE"
			echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/05-Results"
			rm -r "$OUTDIR/04-Classify"
			rm -r "$OUTDIR/05-Results"
		fi

		classify_gtdbtk(){
			# Definiere den Pfad zur Protokolldatei
			CLASSIFY_LOG="$OUTDIR/04-Classify/classify_with_gtdbtk.log"

			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$CLASSIFY_LOG")"

			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte Klassifizierung mit GTDB-Tk" >> "$CLASSIFY_LOG"
			echo "$(datetime.task) Starte Klassifizierung mit GTDB-Tk" >> "$G_LOG_FILE"

			# classify mit GTDB-Tk
			mkdir -p "$OUTDIR/04-Classify"
			mkdir -p "$OUTDIR/04-Classify/GTDBTK"
			gtdbtk classify_wf --genome_dir "$OUTDIR/03-Bins/Bin_Refinement/work_files/binsM" --out_dir "$OUTDIR/04-Classify/GTDBTK" --cpus "$THREADS" --extension fa --skip_ani_screen --min_perc_aa 0.5

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
			echo "$(datetime.done) Klassifizierung mit GTDB-Tk abgeschlossen. Ausgabedaten in $OUTDIR/04-Classify/GTDBTK" >> "$CLASSIFY_LOG"
			echo "$(datetime.done) Klassifizierung mit GTDB-Tk abgeschlossen" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		}

		# Aufruf der Funktion zum Klassifizieren
		classify_gtdbtk

		# Funktion zur Verarbeitung des GTDB-Tk-Ausgangs aufrufen
		process_gtdbtk_output() {
			# Ausgabe des Startzeitpunkts in die Protokolldatei
			echo "$(datetime.task) Starte Aufbereitung der GTDB-Tk Ergebnisse" >> "$CLASSIFY_LOG"
			echo "$(datetime.task) Starte Aufbereitung der GTDB-Tk Ergebnisse" >> "$G_LOG_FILE"

			gtdbtk_output_dir="$1"

			# Ein leeres Array erstellen, um die Anzahl jedes Stammes zu zählen
			declare -A stamm_count

			# Pfad zur TSV-Datei unter Verwendung von gtdbtk_output_dir
			tsv_paths=("$gtdbtk_output_dir"/gtdbtk*.summary.tsv)

			# Erstelle das Verzeichnis für die Ergebnisdatei, falls es nicht vorhanden ist
			results_dir="$OUTDIR/05-Results"
			mkdir -p "$results_dir"

			# CSV-Datei erstellen
			results_csv_path="$results_dir/ergebnisse.csv"
			# Header in die CSV-Datei schreiben
			echo "Stamm,Count" > "$results_csv_path"

			# Für jede TSV-Datei im Verzeichnis iterieren
			for tsv_path in "${tsv_paths[@]}"; do
				# Überprüfe, ob die TSV-Datei vorhanden ist
				if [ -f "$tsv_path" ]; then
					# CSV-Datei bearbeiten
					awk 'NR>1 { print $2 }' FS='\t' OFS=',' "$tsv_path" | sort | uniq -c | awk '{ print $2","$1 }' >>"$results_csv_path"

					echo "Verarbeitung von $tsv_path abgeschlossen."
				else
					echo "Die Datei $tsv_path existiert nicht."
				fi
			done

			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
			echo "$(datetime.done) Aufbereitung GTDB-Tk abgeschlossen. Ausgabedaten in $OUTDIR/05-Results/ergebnisse.csv" | tee -a "$CLASSIFY_LOG"
			echo "$(datetime.done) Aufbereitung GTDB-Tk abgeschlossen" >> "$G_LOG_FILE"
			echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		}

		# Aufruf der Funktion zur Verarbeitung des GTDB-Tk-Ausgangs
		process_gtdbtk_output "$OUTDIR/04-Classify/GTDBTK"
		#-------------------------------------------------------------------------------------------------------------------------------------
		
		# Extrahieren der Stammnamen aus GTDBTK Results
		extract_output() {
			# Eingabedatei
			EINGABEDATEI="$OUTDIR/05-Results/ergebnisse.csv"
			# Ausgabedatei
			AUSGABEDATEI="$OUTDIR/05-Results/extrakt.csv"
			# Temporäre Datei
			TEMP="$OUTDIR/05-Results/temp.csv"
			# Trennzeichen
			TRENNZEICHEN=";"
			# Leere Ausgabedatei erstellen
			touch $AUSGABEDATEI
			# Spalten in der Ergebnisdatei
			SPALTEN="1,2"
			# Header in die Ausgabedatei schreiben
			echo "Stamm,Count" > $AUSGABEDATEI
			# Ergebnisdatei zeilenweise verarbeiten
			tail -n +2 "$EINGABEDATEI" | while IFS=, read -r STAMM COUNT; do
				# Stammname extrahieren
				if [[ $STAMM == "Unclassified" ]]; then
					# "Unclassified" unverändert übernehmen
					NAME=$STAMM
				elif [[ $STAMM =~ ^$TRENNZEICHEN ]]; then
					# Kein Stammname vorhanden, Trennzeichen nach vorne schieben
					NAME="${STAMM:1}"
				else
					# Letzten Namen extrahieren
					NAME=${STAMM##*;}
				fi

				# Leerzeichen im Namen entfernen
				NAME=${NAME//[[:space:]]/}

				# Stammname und Zählwert in temporäre Datei schreiben
				echo "$NAME,$COUNT" >> $AUSGABEDATEI
			done
		}

		# Aufruf der Funktion zur Extraktion des GTDB-Tk-Output
		extract_output
		#-------------------------------------------------------------------------------------------------------------------------------------

		# Shiny App zur Datenauswertung
		shiny_gui() {
			# Verzeichnis, in dem sich die CSV-Datei befindet
			CSV_DIR="$OUTDIR/05-Results"
			# Pfad zur CSV-Datei
			CSV_FILE="$CSV_DIR/extrakt.csv"
			# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
			R_SCRIPT_DIR="/home/drk/atk_shell_based/R"
			# Definiere den Pfad zur Protokolldatei für BLASTn
			RSCRIPT_LOG="$OUTDIR/05-Results/r_plots.log"
			# Erstelle das Verzeichnis, falls es nicht vorhanden ist
			mkdir -p "$(dirname "$RSCRIPT_LOG")"
			# Ausgabe des Startzeitpunkts in die Protokolldatei für BLASTn
			echo "$(datetime.task) Starte R-Skript" >> "$RSCRIPT_LOG"
			# Übergabe der Werte an das R-Skript.
			# Setze den Pfad zur Ausgabedatei als Umgebungsvariable
			export CSV_PATH="$CSV_FILE"
			# Wechsle in das Verzeichnis des R-Skripts
			cd "$R_SCRIPT_DIR" || exit
			# Rufe dein R-Skript auf und übergebe den Pfad zur Ausgabedatei
			Rscript "$R_SCRIPT_DIR/plot_script.R" "$CSV_PATH" "$INFQ" &
			# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
			sleep 2
			# Öffne den Webbrowser mit der Shiny-App
			xdg-open "http://127.0.0.1:4010"
			# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für BLASTn
			echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$RSCRIPT_LOG"
			echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$G_LOG_FILE"
			sleep 120
		}

		# Aufruf der Funktion Shiny GUI
		shiny_gui
		#-------------------------------------------------------------------------------------------------------------------------------------

		;;

    *)

            echo "$(datetime.warning) Pipeline wird beendet."
            exit 1
            ;;
    esac

else	
	echo "$(datetime.info) Starte Analyse-Pipeline SINGLE-RUN (BETA-MODE) bei letzter Analyse-Position"
		# Erstellen einer globalen Log-File
	G_LOG_FILE="$OUTDIR/global_file.log"
		
	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$G_LOG_FILE")"	
					
	# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe
	if [ -d "$OUTDIR" ]; then	
		echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.warning) DAS AUSGABEVERZEICHNIS IST BEREITS VORHANDEN. UEBERSPRINGE DIE SCHON BEARBEITETEN DATEIEN" >> "$G_LOG_FILE"
		echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"	
	else				
		# Erstellen des OUTPUT Verzeichnisses
		mkdir $OUTDIR

	fi

	#-------------------------------------------------------------------------------------------------------------------------------------

	# Wenn der Ordner $OUTDIR/01-CleanedReads bereits existiert, dann überspringe
	if [ -d "$OUTDIR/01-CleanedReads" ]; then
		echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" >> "$G_LOG_FILE"
		echo "$(datetime.skip) Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation wird uebersprungen" 
	else	    	
			#### Kopieren der Input Datei in das Output Verzeichnis zur Dokumentation
		mkdir $OUTDIR/00-InputData

		# Ausgabe des Startzeitpunkts in die Protokolldatei
		echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" 
		echo "$(datetime.task) Starte Analyse-Pipeline BETA-Mode" >> "$G_LOG_FILE"

		# Kopiere Input Fastq in OutDIR
		echo "$(datetime.task) Kopiere die Input-Fastq Datei in das Ausgabeverzeichnis" >> "$G_LOG_FILE"
		cp "$INFQ" "$OUTDIR/00-InputData"
		echo "$(datetime.done) $INFQ in $OUTDIR/00-InputData erfolgreich kopiert" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	fi
	#-------------------------------------------------------------------------------------------------------------------------------------

	# Wenn der Ordner $OUTDIR/02-De_Novo_Assembly bereits existiert, dann überspringe
	if [ -d "$OUTDIR/02-De_Novo_Assembly" ]; then
		echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.warning) Alignment gegen Human Referenzgenom wird uebersprungen" >> "$G_LOG_FILE"
		echo "$(datetime.skip) Alignment gegen Human Referenzgenom wird uebersprungen" 
	else
		#### Alignment gegen Human Referenzgenom
		# Definiere den Pfad zur Protokolldatei
		ALIGN_LOG="$OUTDIR/01-CleanedReads/align_to_human_ref.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$ALIGN_LOG")"

		# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
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

	# Wenn der Ordner $OUTDIR/03-Bins bereits existiert, dann überspringe
	if [ -d "$OUTDIR/03-Bins" ]; then
		echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.skip) DeNovo Assembly wird uebersprungen" >> "$G_LOG_FILE"
			echo "$(datetime.skip) DeNovo Assembly wird uebersprungen" 
	else
			### DeNovo Assembly
			# Definiere den Pfad zur Protokolldatei
		FLYE_LOG="$OUTDIR/02-De_Novo_Assembly/de_novo_assembly.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$ALIGN_LOG")"

		# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
		if [ -z "$OUTDIR" ] || [ -z "$THREADS" ] || [ -z "$REFHUM" ]; then
			echo "$(datetime.error) Erforderliche Variablen fehlen."
			exit 1
		fi
		

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$FLYE_LOG")"

		# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
		if [ -z "$OUTDIR" ] || [ -z "$THREADS" ]; then
			echo "$(datetime.error) Erforderliche Variablen fehlen."
			exit 1
		fi

		# Ausgabe des Startzeitpunkts in die Protokolldatei
		echo "$(datetime.task) Starte De-Novo-Assembly" >> "$FLYE_LOG"
		echo "$(datetime.task) Starte De-Novo-Assembly" >> "$G_LOG_FILE"

		# Assembly mit Flye durchführen
		mkdir -p "$OUTDIR/02-De_Novo_Assembly"
		flye --meta --nano-hq "$OUTDIR/01-CleanedReads/align_to_ref_unaligned.fastq" -t "$THREADS" -o "$OUTDIR/02-De_Novo_Assembly/Flye"

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
		echo "$(datetime.done) De-Novo-Assembly abgeschlossen. Ausgabedaten in $OUTDIR/02-De_Novo_Assembly/Flye" | tee -a "$FLYE_LOG"
		echo "$(datetime.done) De-Novo-Assembly abgeschlossen. Ausgabedaten in $OUTDIR/02-De_Novo_Assembly/Flye" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.done)  _____  _      __ __    ___      ___     ___   ____     ___ "
		echo "$(datetime.done) |     || |    |  |  |  /  _]    |   \   /   \ |    \   /  _]"
		echo "$(datetime.done) |   __|| |    |  |  | /  [_     |    \ |     ||  _  | /  [_ "
		echo "$(datetime.done) |  |_  | |___ |  ~  ||    _]    |  D  ||  O  ||  |  ||    _]"
		echo "$(datetime.done) |   _] |     ||___, ||   [_     |     ||     ||  |  ||   [_ "
		echo "$(datetime.done) |  |   |     ||     ||     |    |     ||     ||  |  ||     |"
		echo "$(datetime.done) |__|   |_____||____/ |_____|    |_____| \___/ |__|__||_____|"
	fi

	#-------------------------------------------------------------------------------------------------------------------------------------

	# Wenn der Ordner $OUTDIR/04-Classify bereits existiert, dann überspringe
	if [ -d "$OUTDIR/04-Classify" ]; then
		echo "$(datetime.warning) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.skip) Binning mit Metawrap & Binning-Verfeinerung mit Metawrap wird uebersprungen" >> "$G_LOG_FILE"
			echo "$(datetime.skip) Binning mit Metawrap & Binning-Verfeinerung mit Metawrap wird uebersprungen" 
	else
			### Binning mit Metawrap
		# Definiere den Pfad zur Protokolldatei
		BINNING_LOG="$OUTDIR/03-Bins/binning.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$BINNING_LOG")"

		# Überprüfe, ob alle erforderlichen Variablen gesetzt sind
		if [ -z "$OUTDIR" ] || [ -z "$THREADS" ]; then
			echo "$(datetime.error) Erforderliche Variablen fehlen."
			exit 1
		fi

		# Ausgabe des Startzeitpunkts in die Protokolldatei
		echo "$(datetime.task) Starte Binning mit Metawrap" >> "$BINNING_LOG"
		echo "$(datetime.task) Starte Binning mit Metawrap" >> "$G_LOG_FILE"

		# Metawrap Binning
		mkdir -p "$OUTDIR/03-Bins"
		mkdir -p "$OUTDIR/03-Bins/Binning"

		# Binning mit Metawrap
		metawrap binning -t "$THREADS" -a "$OUTDIR/02-De_Novo_Assembly/Flye/assembly.fasta" -o "$OUTDIR/03-Bins/Binning" --metabat2 --maxbin2 --single-end "$OUTDIR/00-InputData/"*fastq

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
		echo "$(datetime.done) Binning mit Metawrap abgeschlossen. Ausgabedaten in $OUTDIR/03-Bins/Binning" | tee -a "$BINNING_LOG"
		echo "$(datetime.done) Binning mit Metawrap abgeschlossen" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

		
		
		
		### Binning-Verfeinerung mit Metawrap

		# Definiere den Pfad zur Protokolldatei
		BIN_REFINEMENT_LOG="$OUTDIR/03-Bins/bin_refinement.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$BIN_REFINEMENT_LOG")"

		# Ausgabe des Startzeitpunkts in die Protokolldatei
		echo "$(datetime.task) Starte Binning-Verfeinerung mit Metawrap" >> "$BIN_REFINEMENT_LOG"
		echo "$(datetime.task) Starte Binning-Verfeinerung mit Metawrap" >> "$G_LOG_FILE"

		# Metawrap Binning-Verfeinerung
		mkdir -p "$OUTDIR/03-Bins/Bin_Refinement"

		# Binning-Verfeinerung mit Metawrap
		metawrap bin_refinement -t "$THREADS" -A "$OUTDIR/03-Bins/Binning/metabat2_bins" -B "$OUTDIR/03-Bins/Binning/maxbin2_bins" -o "$OUTDIR/03-Bins/Bin_Refinement" --skip-checkm

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
		echo "$(datetime.done) Binning-Verfeinerung mit Metawrap abgeschlossen. Ausgabedaten in $OUTDIR/03-Bins/Bin_Refinement" | tee -a "$BIN_REFINEMENT_LOG"
		echo "$(datetime.done) Binning-Verfeinerung mit Metawrap abgeschlossen" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"

		# Definiere den Pfad zur Protokolldatei für die Split-Funktion
		SPLIT_LOG="$OUTDIR/03-Bins/split_fasta.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$SPLIT_LOG")"

		# Ausgabe des Startzeitpunkts in die Protokolldatei für die Split-Funktion
		echo "$(datetime.task) Starte Split-Funktion" >> "$SPLIT_LOG"
		echo "$(datetime.task) Starte Split-Funktion" >> "$G_LOG_FILE"

		# Funktion zur Aufteilung der unbinned.fa-Datei aufrufen, wenn vorhanden
		split_fasta() {
			refined_bins_dir="$1"
			unbinned_fasta_path=$(find "$refined_bins_dir" -type f -name "bin.unbinned.fa" | head -n 1)

			if [ -z "$unbinned_fasta_path" ]; then
			echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen."
			echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen." >> "$SPLIT_LOG"
			echo "$(datetime.warning) Die Datei bin.unbinned.fa wurde nicht gefunden. Die Split-Funktion wird übersprungen." >> "$G_LOG_FILE"
			return 1
			fi

			contig_count=0
			while IFS= read -r line; do
			if [[ $line == ">"* ]]; then
				# Neue Contig gefunden
				contig_count=$((contig_count+1))
				contig_name="${line:1}"
				printf '%s\n' "$line" > "${unbinned_fasta_path%.fa}.contig${contig_count}.fa"
			else
				# Sequenz der aktuellen Contig speichern
				printf '%s\n' "$line" >> "${unbinned_fasta_path%.fa}.contig${contig_count}.fa"
			fi
			done < "$unbinned_fasta_path"

			echo "$(datetime.done) Die Datei $unbinned_fasta_path wurde erfolgreich aufgespalten."
		}

		# Aufruf der Funktion
		split_fasta "$OUTDIR/03-Bins/Bin_Refinement/work_files/binsM"

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für die Split-Funktion
		echo "$(datetime.done) Split-Funktion abgeschlossen." | tee -a "$SPLIT_LOG"
		echo "$(datetime.done) Split-Funktion abgeschlossen" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------" >> "$G_LOG_FILE"
		echo "$(datetime.done)  ___ ___    ___ ______   ____  __    __  ____    ____  ____       ___     ___   ____     ___ "
		echo "$(datetime.done) |   |   |  /  _]      | /    ||  |__|  ||    \  /    ||    \     |   \   /   \ |    \   /  _]"
		echo "$(datetime.done) | _   _ | /  [_|      ||  o  ||  |  |  ||  D  )|  o  ||  o  )    |    \ |     ||  _  | /  [_ "
		echo "$(datetime.done) |  \_/  ||    _]_|  |_||     ||  |  |  ||    / |     ||   _/     |  D  ||  O  ||  |  ||    _]"
		echo "$(datetime.done) |   |   ||   [_  |  |  |  _  ||        ||    \ |  _  ||  |       |     ||     ||  |  ||   [_ "
		echo "$(datetime.done) |   |   ||     | |  |  |  |  | \      / |  .  \|  |  ||  |       |     ||     ||  |  ||     |"
		echo "$(datetime.done) |___|___||_____| |__|  |__|__|  \_/\_/  |__|\_||__|__||__|       |_____| \___/ |__|__||_____|"
	fi

	#-------------------------------------------------------------------------------------------------------------------------------------

	### Klassifizierung mit GTDB-Tk
	# Prüfe, ob der Ausgabeordner bereits existiert
	if [ -d "$OUTDIR/04-Classify" ]; then
		# Wenn der Ausgabeordner existiert, entferne ihn und vermerke es im Terminal und Log
		echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/04-Classify" >> "$G_LOG_FILE"
		echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/04-Classify"
		echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/05-Results" >> "$G_LOG_FILE"
		echo "$(datetime.warning) Entferne vorhandenen Ausgabeordner $OUTDIR/05-Results"
		rm -r "$OUTDIR/04-Classify"
		rm -r "$OUTDIR/05-Results"
	fi

	# Definiere den Pfad zur Protokolldatei
	CLASSIFY_LOG="$OUTDIR/04-Classify/classify_with_gtdbtk.log"

	# Erstelle das Verzeichnis, falls es nicht vorhanden ist
	mkdir -p "$(dirname "$CLASSIFY_LOG")"

	# Ausgabe des Startzeitpunkts in die Protokolldatei
	echo "$(datetime.task) Starte Klassifizierung mit GTDB-Tk" >> "$CLASSIFY_LOG"
	echo "$(datetime.task) Starte Klassifizierung mit GTDB-Tk" >> "$G_LOG_FILE"

	# classify mit GTDB-Tk
	mkdir -p "$OUTDIR/04-Classify"
	mkdir -p "$OUTDIR/04-Classify/GTDBTK"
	gtdbtk classify_wf --genome_dir "$OUTDIR/03-Bins/Bin_Refinement/work_files/binsM" --out_dir "$OUTDIR/04-Classify/GTDBTK" --cpus "$THREADS" --extension fa --skip_ani_screen --min_perc_aa 0.5

	# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
	echo "$(datetime.done) Klassifizierung mit GTDB-Tk abgeschlossen. Ausgabedaten in $OUTDIR/04-Classify/GTDBTK" >> "$CLASSIFY_LOG"
	echo "$(datetime.done) Klassifizierung mit GTDB-Tk abgeschlossen" >> "$G_LOG_FILE"
	echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"


	# Funktion zur Verarbeitung des GTDB-Tk-Ausgangs aufrufen
	process_gtdbtk_output() {
		# Ausgabe des Startzeitpunkts in die Protokolldatei
		echo "$(datetime.task) Starte Aufbereitung der GTDB-Tk Ergebnisse" >> "$CLASSIFY_LOG"
		echo "$(datetime.task) Starte Aufbereitung der GTDB-Tk Ergebnisse" >> "$G_LOG_FILE"

		gtdbtk_output_dir="$1"

		# Ein leeres Array erstellen, um die Anzahl jedes Stammes zu zählen
		declare -A stamm_count

		# Pfad zur TSV-Datei unter Verwendung von gtdbtk_output_dir
		tsv_paths=("$gtdbtk_output_dir"/gtdbtk*.summary.tsv)

		# Erstelle das Verzeichnis für die Ergebnisdatei, falls es nicht vorhanden ist
		results_dir="$OUTDIR/05-Results"
		mkdir -p "$results_dir"

		# CSV-Datei erstellen
		results_csv_path="$results_dir/ergebnisse.csv"
		# Header in die CSV-Datei schreiben
		echo "Stamm,Count" > "$results_csv_path"

		# Für jede TSV-Datei im Verzeichnis iterieren
		for tsv_path in "${tsv_paths[@]}"; do
			# Überprüfe, ob die TSV-Datei vorhanden ist
			if [ -f "$tsv_path" ]; then
				# CSV-Datei bearbeiten
				awk 'NR>1 { print $2 }' FS='\t' OFS=',' "$tsv_path" | sort | uniq -c | awk '{ print $2","$1 }' >>"$results_csv_path"

				echo "Verarbeitung von $tsv_path abgeschlossen."
			else
				echo "Die Datei $tsv_path existiert nicht."
			fi
		done

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben
		echo "$(datetime.done) Aufbereitung GTDB-Tk abgeschlossen. Ausgabedaten in $OUTDIR/05-Results/ergebnisse.csv" | tee -a "$CLASSIFY_LOG"
		echo "$(datetime.done) Aufbereitung GTDB-Tk abgeschlossen" >> "$G_LOG_FILE"
		echo "$(datetime.done) ------------------------------------------------------------------------------------------------------------------------" >> "$G_LOG_FILE"
	}

	# Aufruf der Funktion zur Verarbeitung des GTDB-Tk-Ausgangs
	process_gtdbtk_output "$OUTDIR/04-Classify/GTDBTK"

	#-------------------------------------------------------------------------------------------------------------------------------------
	### Extrahieren der Stammnamen aus GTDBTK Results
	# Eingabedatei
	EINGABEDATEI="$OUTDIR/05-Results/ergebnisse.csv"

	# Ausgabedatei
	AUSGABEDATEI="$OUTDIR/05-Results/extrakt.csv"

	# Temporäre Datei
	TEMP="$OUTDIR/05-Results/temp.csv"

	# Trennzeichen
	TRENNZEICHEN=";"

	# Leere Ausgabedatei erstellen
	touch $AUSGABEDATEI

	# Spalten in der Ergebnisdatei
	SPALTEN="1,2"

	# Header in die Ausgabedatei schreiben
	echo "Stamm,Count" > $AUSGABEDATEI

	# Ergebnisdatei zeilenweise verarbeiten
	tail -n +2 "$EINGABEDATEI" | while IFS=, read -r STAMM COUNT; do
		# Stammname extrahieren
		if [[ $STAMM == "Unclassified" ]]; then
			# "Unclassified" unverändert übernehmen
			NAME=$STAMM
		elif [[ $STAMM =~ ^$TRENNZEICHEN ]]; then
			# Kein Stammname vorhanden, Trennzeichen nach vorne schieben
			NAME="${STAMM:1}"
		else
			# Letzten Namen extrahieren
			NAME=${STAMM##*;}
		fi

		# Leerzeichen im Namen entfernen
		NAME=${NAME//[[:space:]]/}

		# Stammname und Zählwert in temporäre Datei schreiben
		echo "$NAME,$COUNT" >> $AUSGABEDATEI
	done


	#-------------------------------------------------------------------------------------------------------------------------------------
		# Starte Shiny App zur Datenauswertung
		# Verzeichnis, in dem sich die CSV-Datei befindet
		CSV_DIR="$OUTDIR/05-Results"

		# Pfad zur CSV-Datei
		CSV_FILE="$CSV_DIR/extrakt.csv"

		# Pfad zum Verzeichnis des R-Skripts (Übergeordnetes Verzeichnis)
		R_SCRIPT_DIR="/home/drk/atk_shell_based/R"

		# Definiere den Pfad zur Protokolldatei für BLASTn
		RSCRIPT_LOG="$OUTDIR/05-Results/r_plots.log"

		# Erstelle das Verzeichnis, falls es nicht vorhanden ist
		mkdir -p "$(dirname "$RSCRIPT_LOG")"

		# Ausgabe des Startzeitpunkts in die Protokolldatei für BLASTn
		echo "$(datetime.task) Starte R-Skript" >> "$RSCRIPT_LOG"

		# Übergabe der Werte an das R-Skript.
		# Setze den Pfad zur Ausgabedatei als Umgebungsvariable
		export CSV_PATH="$CSV_FILE"

		# Wechsle in das Verzeichnis des R-Skripts
		cd "$R_SCRIPT_DIR" || exit

		# Rufe dein R-Skript auf und übergebe den Pfad zur Ausgabedatei
		Rscript "$R_SCRIPT_DIR/plot_script.R" "$CSV_PATH" "$INFQ" &

		# Warte eine kurze Zeit, um sicherzustellen, dass die Shiny-App gestartet wurde
		sleep 2

		# Öffne den Webbrowser mit der Shiny-App
		xdg-open "http://127.0.0.1:4010"

		# Ausgabe, dass der Vorgang abgeschlossen ist und Protokoll in die Datei schreiben für BLASTn
		echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$RSCRIPT_LOG"
		echo "$(datetime.done) Übergabe an R-Skript abgeschlossen." | tee -a "$G_LOG_FILE"

		sleep 120
    
