# MCI Cardiomyopathy Resource Browser
# Resource version: v0.2-resource-reframe-corrected-source

library(shiny)
library(DT)
library(ggplot2)

DATA_PATH <- "../results/resource_tables/MCI_BROWSER_READY_DISPLAY_TABLE_v0_2.csv"

if (!file.exists(DATA_PATH)) {
  stop(paste("Could not find browser table at:", DATA_PATH))
}

mci <- read.csv(DATA_PATH, stringsAsFactors = FALSE, check.names = FALSE)

numeric_cols <- c("MCI", "Adj_MCI", "MCI_CI95_lower", "MCI_CI95_upper", "sigma_disease_to_GTEx_ratio", "DCM_MCI_if_available")

for (col in numeric_cols) {
  if (col %in% names(mci)) {
    mci[[col]] <- suppressWarnings(as.numeric(mci[[col]]))
  }
}

ui <- fluidPage(
  titlePanel("MCI Cardiomyopathy Resource Browser"),
  tags$p("A browser-ready resource for evaluating transcriptomic concordance of ClinVar-annotated cardiomyopathy genes across public human heart datasets."),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("gene", "Search gene", choices = sort(unique(mci$gene_symbol)), selected = sort(unique(mci$gene_symbol))[1]),
      selectInput("tier", "Filter by MCI tier", choices = c("ALL", sort(unique(mci$MCI_tier))), selected = "ALL"),
      selectInput("gtex_flag", "Filter by GTEx low-confidence flag", choices = c("ALL", "TRUE", "FALSE"), selected = "ALL"),
      checkboxInput("only_dcm", "Show only genes with DCM generalization entry", value = FALSE),
      tags$hr(),
      downloadButton("download_filtered", "Download filtered table"),
      tags$hr(),
      tags$p(tags$b("Interpretation boundary:")),
      tags$p("MCI is an evidence-auditing resource. It does not prove causality, clinical actionability, or therapeutic validity.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Gene view", br(), h3(textOutput("gene_title")), tableOutput("gene_summary"), br(), h4("Resource interpretation"), verbatimTextOutput("gene_interpretation"), br(), h4("MCI and GTEx context"), plotOutput("gene_plot", height = "300px")),
        tabPanel("Resource table", br(), DTOutput("resource_table")),
        tabPanel("Tier distribution", br(), plotOutput("tier_plot", height = "400px")),
        tabPanel("MCI vs GTEx ratio", br(), plotOutput("mci_gtex_plot", height = "450px")),
        tabPanel("About", br(), h3("About this resource"),
                 tags$p("The Molecular Concordance Index combines direction agreement, effect-size consistency, and statistical reproducibility across disease cohorts."),
                 tags$p("GTEx baseline benchmarking compares disease-associated variability against normal left-ventricle expression variability."),
                 tags$p("A GTEx-low-confidence flag means disease-associated variability does not exceed GTEx baseline variability in the current bulk data. It does not mean the gene is biologically irrelevant."),
                 tags$p("This browser is designed for database/resource manuscript framing and gene-level evidence review."))
      )
    )
  )
)

server <- function(input, output, session) {

  filtered_data <- reactive({
    df <- mci
    if (!is.null(input$tier) && input$tier != "ALL") { df <- df[df$MCI_tier == input$tier, ] }
    if (!is.null(input$gtex_flag) && input$gtex_flag != "ALL") {
      target_flag <- ifelse(input$gtex_flag == "TRUE", "TRUE", "FALSE")
      df <- df[toupper(as.character(df$GTEx_low_confidence_flag)) == target_flag, ]
    }
    if (isTRUE(input$only_dcm)) {
      df <- df[toupper(as.character(df$has_DCM_generalization_entry)) == "TRUE", ]
    }
    df
  })

  selected_gene_data <- reactive({
    df <- mci[mci$gene_symbol == input$gene, ]
    if (nrow(df) == 0) return(NULL)
    df[1, ]
  })

  output$gene_title <- renderText({ paste("Gene:", input$gene) })

  output$gene_summary <- renderTable({
    row <- selected_gene_data()
    if (is.null(row)) return(NULL)
    keep <- c("gene_symbol", "disease_group", "stratum", "MCI", "Adj_MCI", "MCI_tier", "bootstrap_majority_tier", "MCI_CI95_lower", "MCI_CI95_upper", "sigma_disease_to_GTEx_ratio", "GTEx_low_confidence_flag", "has_DCM_generalization_entry", "DCM_MCI_if_available", "DCM_tier_if_available")
    keep <- keep[keep %in% names(row)]
    out <- data.frame(field = keep, value = as.character(unlist(row[keep])), stringsAsFactors = FALSE)
    out
  }, striped = TRUE, bordered = TRUE)

  output$gene_interpretation <- renderText({
    row <- selected_gene_data()
    if (is.null(row)) return("No gene selected.")
    if ("resource_interpretation" %in% names(row)) return(row$resource_interpretation)
    "No interpretation field available."
  })

  output$gene_plot <- renderPlot({
    row <- selected_gene_data()
    if (is.null(row)) return(NULL)
    vals <- data.frame(metric = c("MCI", "Adj_MCI", "DCM MCI"), value = c(row$MCI, row$Adj_MCI, row$DCM_MCI_if_available))
    vals <- vals[!is.na(vals$value), ]
    ggplot(vals, aes(x = metric, y = value)) +
      geom_col() +
      geom_hline(yintercept = 0.70, linetype = 'dashed') +
      geom_hline(yintercept = 0.45, linetype = 'dashed') +
      ylim(0, max(1.2, max(vals$value, na.rm = TRUE))) +
      labs(x = "", y = "Score", title = paste("MCI profile for", row$gene_symbol)) +
      theme_minimal(base_size = 13)
  })

  output$resource_table <- renderDT({
    datatable(filtered_data(), options = list(pageLength = 15, scrollX = TRUE, autoWidth = TRUE), rownames = FALSE, filter = 'top')
  })

  output$tier_plot <- renderPlot({
    df <- filtered_data()
    counts <- as.data.frame(table(df$MCI_tier), stringsAsFactors = FALSE)
    names(counts) <- c("MCI_tier", "n_genes")
    ggplot(counts, aes(x = MCI_tier, y = n_genes)) +
      geom_col() +
      labs(x = "MCI tier", y = "Number of genes", title = "MCI tier distribution") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  output$mci_gtex_plot <- renderPlot({
    df <- filtered_data()
    df <- df[!is.na(df$MCI) & !is.na(df$sigma_disease_to_GTEx_ratio), ]
    ggplot(df, aes(x = MCI, y = sigma_disease_to_GTEx_ratio)) +
      geom_point() +
      geom_hline(yintercept = 1.0, linetype = 'dashed') +
      geom_vline(xintercept = 0.70, linetype = 'dashed') +
      geom_vline(xintercept = 0.45, linetype = 'dashed') +
      labs(x = "MCI", y = "sigma_disease / sigma_GTEx", title = "MCI versus GTEx baseline variability ratio") +
      theme_minimal(base_size = 13)
  })

  output$download_filtered <- downloadHandler(
    filename = function() { paste0('MCI_filtered_resource_table_', Sys.Date(), '.csv') },
    content = function(file) { write.csv(filtered_data(), file, row.names = FALSE) }
  )
}

shinyApp(ui = ui, server = server)
