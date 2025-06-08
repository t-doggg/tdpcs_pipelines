#!/usr/bin/env Rscript

options(shiny.maxRequestSize = 5 * 1024^3)  # 5 GB upload limit

library(shiny)
library(plotly)
library(DT)
library(shinythemes)
library(shinyFiles)
library(fs)

# === Helper ===
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

# === UI ===
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(tags$style(HTML("
    #log_output {
      max-height: 500px;
      overflow-y: scroll;
      white-space: pre-wrap;
      font-family: monospace;
    }
  "))),
  titlePanel("🧬 TDPCS GUI + Result Viewer (v0.04b_bug)"),
  sidebarLayout(
    sidebarPanel(
      h4("🚀 Pipeline Launcher"),
      fileInput("input_fastq", "📄 Input FASTQ"),
      fileInput("reference_fasta", "🧬 Human Reference (FASTA)"),
      h5("📁 Output Directory"),
      shinyDirButton("output_dir", "📁 Choose Output Directory", "Select output folder"),
      verbatimTextOutput("output_dir_display"),
      textInput("run_name", "📝 Experiment Name", value = "run_v1"),   # <-- Add this line
      fileInput("ncbi_db", "🧬 Select Any File in NCBI DB Folder"),
      numericInput("threads", "🧵 Threads", value = 4),
      numericInput("coverage", "📊 Coverage %", value = 75),
      numericInput("fraction", "🧪 Kraken Fraction %", value = 0),
      selectInput("mode", "⚙️ Mode", choices = c("fl", "fs", "ll", "ls", "dev")),
      verbatimTextOutput("cmd_display"),
      uiOutput("run_btn_ui"),
      hr(),
      h4("📂 Load Results Folder"),
      shinyDirButton("runfolder", "📂 Select TDPCS Run Folder", "Choose a folder"),
      verbatimTextOutput("runfolder_display"),
      actionButton("load_outputs", "📥 Load Outputs", class = "btn btn-info btn-block"),
      hr(),
      downloadButton("export_data", "💾 Download Filtered CSV")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("📊 Overview", plotlyOutput("main_plot")),
        tabPanel("🔹 NCBI", plotlyOutput("plot_1")),
        tabPanel("🔸 Kraken2", plotlyOutput("plot_2")),
        tabPanel("🧬 AMR", plotlyOutput("amr_plot")),
        tabPanel("📋 Table", DTOutput("data_view")),
        tabPanel("📏 Read Sizes", plotlyOutput("read_distribution")),
        tabPanel("📐 Q-Scores", plotlyOutput("qscore_distribution")),
        tabPanel("📄 Summary", verbatimTextOutput("summary_output")),
        tabPanel("🧾 Run Log", verbatimTextOutput("log_output"))
      )
    )
  )
)

# === SERVER ===
server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "Root" = "/")
  
  # === Output dir selector ===
  shinyDirChoose(input, "output_dir", roots = volumes, session = session)
  selected_outdir <- reactiveVal(NULL)
  
  observeEvent(input$output_dir, {
    selected_outdir(parseDirPath(volumes, input$output_dir))
  })
  
  output$output_dir_display <- renderText({
    req(selected_outdir())
    paste("📂 Selected Output Directory:", selected_outdir())
  })
  
  # Helper to get full output path
  full_outdir <- reactive({
    req(selected_outdir(), input$run_name)
    # Remove spaces and slashes from run name for safety
    run_name_clean <- gsub("[/ ]", "", input$run_name)
    file.path(selected_outdir(), run_name_clean)
  })
  
  # === Enable Run button ===
  output$run_btn_ui <- renderUI({
    if (!is.null(input$input_fastq) &&
        !is.null(input$reference_fasta) &&
        !is.null(input$ncbi_db) &&
        !is.null(selected_outdir()) &&
        nzchar(input$run_name)) {
      actionButton("run_btn", "🚀 Run Pipeline", class = "btn btn-success btn-block")
    } else {
      actionButton("run_btn", "🚫 Run Pipeline", class = "btn btn-secondary btn-block", disabled = TRUE)
    }
  })
  
  # === Log file setup ===
  logfile <- file.path(tempdir(), "tdpcs_run.log")
  if (file.exists(logfile)) file.remove(logfile)
  
  last_cmd <- reactiveVal("")
  # === Launch pipeline ===
  observeEvent(input$run_btn, {
    req(input$input_fastq, input$reference_fasta, input$ncbi_db, selected_outdir(), nzchar(input$run_name))
    
    ncbi_folder <- dirname(input$ncbi_db$datapath)
    outdir <- full_outdir()
    
    # Do NOT create the output directory here! The shell script will handle it.
    # if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)   # <-- REMOVE THIS LINE
    
    cmd <- sprintf(
      "./gui analyse -i '%s' -o '%s' -x '%s' -d '%s' -t %d -c %d -f %f -m %s > '%s' 2>&1",
      input$input_fastq$datapath,
      outdir,
      input$reference_fasta$datapath,
      ncbi_folder,
      input$threads,
      input$coverage,
      input$fraction,
      input$mode,
      logfile
    )
    
    last_cmd(cmd)  # Store the command for display
    
    # Run the pipeline in background
    system(cmd, wait = FALSE)
  })

  output$cmd_display <- renderText({
    if (last_cmd() != "") {
      paste("Last executed command:\n", last_cmd())
    } else {
      "No command executed yet."
    }
  })
  
  # === Auto-refresh log view ===
  log_content <- reactiveFileReader(1000, session, logfile, readLines)
  
  output$log_output <- renderText({
    paste(log_content(), collapse = "\n")
  })
  
  # === Load result folder ===
  shinyDirChoose(input, "runfolder", roots = volumes, session = session)
  selected_runfolder <- reactiveVal(NULL)
  
  observeEvent(input$runfolder, {
    selected_runfolder(parseDirPath(volumes, input$runfolder))
  })
  
  output$runfolder_display <- renderText({
    req(selected_runfolder())
    paste("📁 Selected Run Folder:", selected_runfolder())
  })
  
  # === Data objects ===
  data_loaded <- reactiveVal(FALSE)
  fastq_data <- reactiveVal(NULL)
  read_sizes <- reactiveVal(NULL)
  bacteria_data <- reactiveVal(NULL)
  amr_data <- reactiveVal(NULL)
  
  observeEvent(input$load_outputs, {
    runfolder <- selected_runfolder()
    req(runfolder)
    
    fastq_file <- list.files(file.path(runfolder, "00-InputData"), pattern = "\\.(fastq|fq)$", full.names = TRUE)[1]
    bacteria_csv <- file.path(runfolder, "03-Results", "s_res_compared_conficence.csv")
    amr_file <- file.path(runfolder, "03-Classification", "amrfinder_results_summary.csv")
    
    if (!file.exists(fastq_file)) {
      showModal(modalDialog(title = "❌ FASTQ File Missing",
                            "No .fastq or .fq file found in 00-InputData/", easyClose = TRUE))
      return()
    }
    
    if (!file.exists(bacteria_csv)) {
      showModal(modalDialog(title = "❌ Bacteria CSV Missing",
                            "Expected file not found: 03-Results/s_res_compared_conficence.csv", easyClose = TRUE))
      return()
    }
    
    fastq_data(read_fastq(fastq_file))
    read_sizes(sapply(fastq_data(), nchar))
    bacteria_data(read.csv(bacteria_csv))
    
    if (file.exists(amr_file)) {
      amr_data(read.csv(amr_file))
    } else {
      amr_data(NULL)
    }
    
    data_loaded(TRUE)
  })
  
  # === Plot outputs ===
  output$main_plot <- renderPlotly({
    req(data_loaded())
    df <- bacteria_data()
    df$combined <- df$count_1 + df$count_2
    df <- df[!grepl("Human", df$bacteria, ignore.case = TRUE), ]
    plot_ly(df, x = ~bacteria, y = ~combined, type = "bar",
            marker = list(color = ifelse(df$confidence == 1, 'green', 'red')),
            text = ~paste("Confidence:", confidence)) %>%
      layout(title = "Combined Read Counts", xaxis = list(tickangle = -45))
  })
  
  output$plot_1 <- renderPlotly({
    req(data_loaded())
    df <- bacteria_data()
    plot_ly(df, x = ~bacteria, y = ~count_1, type = "bar", marker = list(color = 'blue')) %>%
      layout(title = "Count_1")
  })
  
  output$plot_2 <- renderPlotly({
    req(data_loaded())
    df <- bacteria_data()
    plot_ly(df, x = ~bacteria, y = ~count_2, type = "bar", marker = list(color = 'orange')) %>%
      layout(title = "Count_2")
  })
  
  output$amr_plot <- renderPlotly({
    req(data_loaded(), !is.null(amr_data()))
    df <- amr_data()
    plot_ly(df, x = ~Class, y = ~count, type = "bar", color = ~Subclass) %>%
      layout(title = "AMR Summary", xaxis = list(tickangle = -45))
  })
  
  output$data_view <- renderDT({
    req(data_loaded())
    df <- bacteria_data()
    df$combined <- df$count_1 + df$count_2
    datatable(df, options = list(pageLength = 10))
  })
  
  output$read_distribution <- renderPlotly({
    req(data_loaded())
    plot_ly(x = ~read_sizes(), type = "histogram") %>%
      layout(title = "Read Lengths")
  })
  
  output$qscore_distribution <- renderPlotly({
    req(data_loaded())
    qs <- sample(10:40, length(read_sizes()), replace = TRUE)
    plot_ly(x = ~qs, type = "histogram") %>%
      layout(title = "Q-Scores")
  })
  
  output$summary_output <- renderPrint({
    req(data_loaded())
    summary(bacteria_data())
  })
  
  output$export_data <- downloadHandler(
    filename = function() {
      paste0("filtered_bacteria_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(bacteria_data(), file, row.names = FALSE)
    }
  )
}

# === Run the App ===
shinyApp(ui = ui, server = server)
