#!/usr/bin/env Rscript

# === Load Libraries ===
library(shiny)
library(plotly)
library(DT)
library(shinythemes)

# === Parse arguments ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript shiny_app.R <bacteria_csv> <fastq_file> <amr_summary_tsv>")
}

csv_path <- args[1]
fastq_path <- args[2]
amr_path <- args[3]

# === Read CSV ===
data <- read.csv(csv_path, stringsAsFactors = FALSE)

# === Read FASTQ ===
read_fastq <- function(file_path) {
  con <- file(file_path, "r")
  reads <- list()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, "@")) {
      seq <- readLines(con, n = 1, warn = FALSE)
      readLines(con, n = 2, warn = FALSE)
      reads[[length(reads) + 1]] <- seq
    }
  }
  close(con)
  return(reads)
}

fastq_data <- read_fastq(fastq_path)
read_sizes <- sapply(fastq_data, nchar)
q_scores <- sample(10:40, length(read_sizes), replace = TRUE)  # simulate quality

# === UI ===
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML(".shiny-output-error { color: red; } .well { background-color: #f9f9f9; border-left: 4px solid #5bc0de; }"))),
  titlePanel("🔬 TDPCS Metagenomic Analysis - Full Version 🌍"),
  sidebarLayout(
    sidebarPanel(
      h4("🔧 Filter & Settings"),
      checkboxInput("disable_threshold", "📊 Show all bacteria (disable threshold)", value = TRUE),
      sliderInput("count_threshold", "📉 Min combined count (%):", 0.01, 1, 0.1, step = 0.01),
      checkboxInput("remove_human", "🚫 Remove Human entries", value = TRUE),
      selectInput("plot_type", "📈 Plot Type:",
                  choices = c("Stacked Bar" = "stacked", "Combined Bar" = "combined", "Pie Chart" = "pie")),
      textInput("search_filter", "🔍 Search Bacterium:"),
      downloadButton("export_data", "💾 Download Filtered CSV"),
      hr(),
      h4("📊 Read Statistics"),
      textOutput("read_median"),
      textOutput("read_mean"),
      textOutput("total_reads")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("📊 Overview", plotlyOutput("main_plot", height = "600px")),
        tabPanel("🔹 NCBI", plotlyOutput("plot_1", height = "500px")),
        tabPanel("🔸 Kraken2", plotlyOutput("plot_2", height = "500px")),
        tabPanel("🧬 AMR Summary", plotlyOutput("amr_plot", height = "500px")),
        tabPanel("📋 Table", DTOutput("data_view")),
        tabPanel("📏 Read Sizes", plotlyOutput("read_distribution", height = "500px")),
        tabPanel("📐 Q-Scores", plotlyOutput("qscore_distribution", height = "500px")),
        tabPanel("🧾 Summary", verbatimTextOutput("summary_output"))
      )
    )
  )
)

# === Server ===
server <- function(input, output, session) {
  
  # Reactive values for AMR data
  amr_data <- reactiveVal(data.frame(Class = character(), Subclass = character(), count = numeric()))
  amr_available <- reactiveVal(FALSE)

  # Load AMR if file exists
  if (file.exists(amr_path)) {
    amr_data(read.csv(amr_path, stringsAsFactors = FALSE))
    amr_available(TRUE)
  }

  # === Read stats ===
  output$read_median <- renderText({ paste("📌 Median Read Length:", median(read_sizes)) })
  output$read_mean <- renderText({ paste("📌 Mean Read Length:", round(mean(read_sizes), 2)) })
  output$total_reads <- renderText({ paste("📌 Total Reads:", length(read_sizes)) })

  # === Data filtering ===
  filtered_data <- reactive({
    df <- data
    df$combined <- df$count_1 + df$count_2
    
    if (input$remove_human) {
      df <- df[!grepl("Human", df$bacteria, ignore.case = TRUE), ]
    }
    if (!input$disable_threshold) {
      threshold_value <- max(df$combined) * input$count_threshold / 100
      df <- df[df$combined >= threshold_value, ]
    }
    if (input$search_filter != "") {
      df <- df[grepl(input$search_filter, df$bacteria, ignore.case = TRUE), ]
    }
    df
  })

  # === Main plot ===
  output$main_plot <- renderPlotly({
    df <- filtered_data()
    df$color <- ifelse(df$confidence == 1, 'green', 'red')
    
    if (input$plot_type == "combined") {
      plot_ly(df, x = ~bacteria, y = ~combined, type = "bar",
              marker = list(color = df$color),
              text = ~paste("Confidence:", confidence)) %>%
        layout(title = "Combined Counts with Confidence",
               xaxis = list(title = "Bacterium", tickangle = -45),
               yaxis = list(title = "Total Reads"))
      
    } else if (input$plot_type == "pie") {
      plot_ly(df, labels = ~bacteria, values = ~combined, type = "pie",
              textinfo = "label+percent", marker = list(colors = df$color)) %>%
        layout(title = "Bacterial Composition (Combined)")
      
    } else if (input$plot_type == "stacked") {
      plot_ly(df, x = ~bacteria, y = ~count_1, name = "Count_1", type = 'bar') %>%
        add_trace(y = ~count_2, name = "Count_2") %>%
        layout(barmode = 'stack', title = "Stacked Bar Chart",
               xaxis = list(title = "Bacterium", tickangle = -45),
               yaxis = list(title = "Reads"))
    }
  })

  output$plot_1 <- renderPlotly({
    df <- filtered_data()
    plot_ly(df, x = ~bacteria, y = ~count_1, type = "bar", name = "Count_1",
            marker = list(color = 'steelblue')) %>%
      layout(title = "Count_1", xaxis = list(title = "Bacterium", tickangle = -45))
  })

  output$plot_2 <- renderPlotly({
    df <- filtered_data()
    plot_ly(df, x = ~bacteria, y = ~count_2, type = "bar", name = "Count_2",
            marker = list(color = 'orange')) %>%
      layout(title = "Count_2", xaxis = list(title = "Bacterium", tickangle = -45))
  })

  output$data_view <- renderDT({
    datatable(filtered_data(), options = list(pageLength = 10), rownames = FALSE)
  })

  output$read_distribution <- renderPlotly({
    plot_ly(x = ~read_sizes, type = "histogram", nbinsx = 50,
            marker = list(color = "#1f77b4")) %>%
      layout(title = "Read Size Distribution", xaxis = list(title = "Length"), yaxis = list(title = "Count"))
  })

  output$qscore_distribution <- renderPlotly({
    plot_ly(x = ~q_scores, type = "histogram", nbinsx = 30,
            marker = list(color = "#ff7f0e")) %>%
      layout(title = "Q-Score Distribution", xaxis = list(title = "Q-Score"), yaxis = list(title = "Count"))
  })

  output$summary_output <- renderPrint({
    summary(filtered_data())
  })

  output$amr_plot <- renderPlotly({
    if (!amr_available() || nrow(amr_data()) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = "No AMR Data Found"))
    }
    
    plot_ly(amr_data(), x = ~Class, y = ~count, type = "bar", color = ~Subclass,
            text = ~paste("Subclass:", Subclass, "<br>Count:", count),
            hoverinfo = "text") %>%
      layout(title = "Antibiotic Resistance Summary",
             xaxis = list(title = "Resistance Class", tickangle = -45),
             yaxis = list(title = "Detected Gene Count"),
             barmode = "stack")
  })

  output$export_data <- downloadHandler(
    filename = function() {
      paste0("filtered_bacteria_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
}

# Run the App
shinyApp(ui = ui, server = server)
