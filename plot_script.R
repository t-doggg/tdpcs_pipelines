# TDPCS-Analysis - Enhanced Version

# Lade die benötigten Pakete
library("shiny")
library("RColorBrewer")
library("shinythemes")
library("plotly")
library("DT")

# Lies die Argumente aus
args <- commandArgs(trailingOnly = TRUE)

# Argument 1: Pfad zur CSV-Datei
csv_path <- args[1]

# Argument 2: Pfad zur FASTQ-Datei
fastq_path <- args[2]

# Debugging: Überprüfe, ob die Pfade korrekt übergeben wurden
print(paste("CSV Path:", csv_path))
print(paste("FASTQ Path:", fastq_path))

# CSV-Datei einlesen
data <- read.csv(csv_path)

# Funktion zum Lesen der FASTQ-Datei
read_fastq <- function(file_path) {
  con <- file(file_path, "r")
  reads <- list()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (substr(line, 1, 1) == "@") {
      seq_id <- line
      seq <- readLines(con, n = 1, warn = FALSE)
      plus <- readLines(con, n = 1, warn = FALSE)
      qscore <- readLines(con, n = 1, warn = FALSE)
      reads[[length(reads) + 1]] <- list(seq_id = seq_id, seq = seq, qscore = qscore)
    }
  }
  close(con)
  return(reads)
}

# Lese die FASTQ-Daten ein
fastq_data <- read_fastq(fastq_path)

# Extrahiere die Read-Größen und Q-Scores
read_sizes <- sapply(fastq_data, function(read) nchar(read$seq))
q_scores <- unlist(lapply(fastq_data, function(read) utf8ToInt(read$qscore) - 33))

# Berechne zusätzliche statistische Werte
read_median <- median(data$Count)
read_mean <- round(mean(data$Count), 2)
total_reads <- sum(data$Count)

# UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  titlePanel("TDPCS-Analysis - Enhanced Version"),
  sidebarLayout(
    sidebarPanel(
      h4("Analyse-Einstellungen"),
      radioButtons("plot_type", "Plot-Typ:",
                   choices = c("Balkendiagramm" = "bar", "Kuchen-Chart" = "pie", "Histogramm" = "histogram"),
                   selected = "bar"),
      selectInput("sort_by", "Sortieren nach:", choices = c("Stamm", "Count")),
      sliderInput("count_threshold", "Count-Schwellenwert (%):", min = 0.1, max = 5, value = 0.1, step = 0.05, post = "%"),
      checkboxInput("remove_human", "Human entfernen", value = FALSE),
      selectInput("filter_stem", "Nach Stamm filtern:", choices = c("Alle", unique(data$Stamm))),
      fileInput("file", "Daten exportieren", accept = c('.csv')),
      actionButton("export_button", "Exportieren"),
      radioButtons("x_axis_scale", "X-Achsen-Skalierung:",
                   choices = c("1kb" = "1kb", "2kb" = "2kb", "5kb" = "5kb", "Bis zum größten Read" = "max"),
                   selected = "max"),
      tags$h4("Statistische Daten"),
      textOutput("read_median"),
      textOutput("read_mean"),
      textOutput("total_reads"),
      hr(),
      actionButton("refresh_button", "Daten aktualisieren")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotlyOutput("plot", height = "600px")),
        tabPanel("Data View", DTOutput("data_view")),
        tabPanel("Histogramm Read-Größe", plotlyOutput("original_read_distribution", height = "600px")),
        tabPanel("Q-Score-Verteilung", plotlyOutput("qscore_distribution", height = "600px")),
        tabPanel("Zusammenfassung", verbatimTextOutput("summary_output"))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Berechne die statistischen Werte basierend auf read_sizes
  read_median <- median(read_sizes, na.rm = TRUE)
  read_mean <- round(mean(read_sizes, na.rm = TRUE), 2)
  total_reads <- length(read_sizes)
  
  # Setze die statistischen Werte in die UI-Elemente
  output$read_median <- renderText({ paste("Read Median:", read_median) })
  output$read_mean <- renderText({ paste("Read Durchschnitt:", read_mean) })
  output$total_reads <- renderText({ paste("Gesamtanzahl der Reads:", total_reads) })
  
  filtered_data <- reactive({
    threshold_value <- (input$count_threshold / 100) * total_reads
    filtered <- data[data$Count > threshold_value, ]
    if (input$remove_human) {
      filtered <- filtered[filtered$Stamm != "Human", ]
    }
    if (input$filter_stem != "Alle") {
      filtered <- filtered[grepl(input$filter_stem, filtered$Stamm), ]
    }
    return(filtered)
  })
  
  output$plot <- renderPlotly({
    filtered_data <- filtered_data()
    if (input$plot_type == "bar") {
      if (input$sort_by == "Stamm") {
        filtered_data <- filtered_data[order(filtered_data$Stamm), ]
      } else {
        filtered_data <- filtered_data[order(filtered_data$Count), ]
      }
      
      plot_ly(
        x = ~filtered_data$Stamm,
        y = ~filtered_data$Count,
        type = "bar",
        marker = list(color = brewer.pal(n = length(filtered_data$Stamm), name = "Set3"))
      ) %>% layout(title = "Balkendiagramm", xaxis = list(title = "Stamm"), yaxis = list(title = "Count"))
    } else if (input$plot_type == "pie") {
      agg_data <- aggregate(Count ~ Stamm, data = filtered_data, FUN = sum)
      filtered_data <- agg_data[agg_data$Count > (input$count_threshold / 100) * total_reads, ]
      filtered_data <- filtered_data[order(filtered_data$Count, decreasing = TRUE), ]
      
      plot_ly(
        labels = ~filtered_data$Stamm,
        values = ~filtered_data$Count,
        type = "pie",
        marker = list(colors = brewer.pal(n = length(filtered_data$Stamm), name = "Set3"))
      ) %>% layout(title = "Kuchen-Chart")
    } else if (input$plot_type == "histogram") {
      plot_ly(
        x = ~filtered_data$Count,
        type = "histogram",
        marker = list(color = "skyblue")
      ) %>% layout(title = "Histogramm der Counts", xaxis = list(title = "Count"), yaxis = list(title = "Frequency"))
    }
  })
  
  output$original_read_distribution <- renderPlotly({
    xlim <- switch(input$x_axis_scale,
                   "1kb" = 1000,
                   "2kb" = 2000,
                   "5kb" = 5000,
                   "max" = max(read_sizes, na.rm = TRUE))
    
    plot_data <- read_sizes[read_sizes <= xlim]
    plot_ly(
      x = ~plot_data,
      type = "histogram",
      marker = list(color = "lightgreen")
    ) %>% layout(title = "Histogramm der Original-Read-Verteilung", xaxis = list(title = "Read Größe"), yaxis = list(title = "Häufigkeit"))
  })
  
  output$qscore_distribution <- renderPlotly({
    plot_ly(
      x = ~q_scores,
      type = "histogram",
      marker = list(color = "salmon")
    ) %>% layout(title = "Histogramm der Q-Score-Verteilung", xaxis = list(title = "Q-Score"), yaxis = list(title = "Häufigkeit"))
  })
  
  output$summary_output <- renderPrint({
    summary(filtered_data())
  })
  
  observeEvent(input$export_button, {
    if (!is.null(input$file)) {
      export_path <- paste0(dirname(input$file$datapath), "/", sub("\\..*", "", input$file$name))
      
      # Exportieren der Daten
      write.csv(filtered_data(), paste0(export_path, ".csv"))
    }
  })
  
  observeEvent(input$refresh_button, {
    session$reload()
  })
}

# Shiny-App starten
shinyApp(
  ui = ui,
  server = server,
  options = list(port = 4010)
)