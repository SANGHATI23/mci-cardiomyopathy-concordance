
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(scales)

mci <- read_csv("data/mci_concordance_table.csv", show_col_types = FALSE)
de  <- read_csv("data/mci_per_cohort_de_table.csv", show_col_types = FALSE)

mci <- mci %>%
  mutate(
    gene_symbol = toupper(trimws(gene_symbol)),
    MCI_display = ifelse(is.na(MCI), NA, round(MCI, 3)),
    Adj_MCI_display = ifelse(is.na(Adj_MCI), NA, round(Adj_MCI, 3)),
    CI_display = paste0(
      ifelse(is.na(MCI_CI95_lower), "NA", round(MCI_CI95_lower, 3)),
      " - ",
      ifelse(is.na(MCI_CI95_upper), "NA", round(MCI_CI95_upper, 3))
    ),
    GTEx_status_display = case_when(
      isTRUE(GTEx_low_confidence_flag) ~ "LOW CONFIDENCE: disease variability within GTEx baseline",
      GTEx_low_confidence_flag == TRUE ~ "LOW CONFIDENCE: disease variability within GTEx baseline",
      GTEx_low_confidence_flag == FALSE ~ "Disease variability exceeds GTEx baseline",
      TRUE ~ as.character(GTEx_adjustment_status)
    )
  )

de <- de %>%
  mutate(
    gene_symbol = toupper(trimws(gene_symbol)),
    direction = ifelse(log2FC >= 0, "Up", "Down")
  )

gene_choices <- sort(unique(mci$gene_symbol))

ui <- dashboardPage(
  dashboardHeader(title = "MCI Concordance Browser"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gene Query", tabName = "gene", icon = icon("search")),
      menuItem("Stratum Browser", tabName = "stratum", icon = icon("chart-box")),
      menuItem("Full Table", tabName = "table", icon = icon("table")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    ),
    selectizeInput(
      inputId = "gene",
      label = "Search gene",
      choices = gene_choices,
      selected = ifelse("MYH7" %in% gene_choices, "MYH7", gene_choices[1]),
      options = list(placeholder = "Type gene symbol")
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f7f9fb; }
        .small-box { border-radius: 10px; }
        .box { border-radius: 10px; }
        .main-header .logo { font-weight: bold; }
      "))
    ),
    
    tabItems(
      tabItem(
        tabName = "gene",
        fluidRow(
          valueBoxOutput("mciBox", width = 3),
          valueBoxOutput("adjBox", width = 3),
          valueBoxOutput("tierBox", width = 3),
          valueBoxOutput("gtexBox", width = 3)
        ),
        fluidRow(
          box(width = 12, title = "Gene-level concordance summary", status = "primary", solidHeader = TRUE,
              DTOutput("geneSummary"))
        ),
        fluidRow(
          box(width = 7, title = "Visual concordance profile: per-cohort log2FC", status = "primary", solidHeader = TRUE,
              plotlyOutput("forestPlot", height = "450px")),
          box(width = 5, title = "GTEx baseline overlay", status = "warning", solidHeader = TRUE,
              plotlyOutput("gtexPlot", height = "450px"))
        ),
        fluidRow(
          box(width = 12, title = "Per-cohort DE statistics for selected gene", status = "primary", solidHeader = TRUE,
              DTOutput("deTable"))
        )
      ),
      
      tabItem(
        tabName = "stratum",
        fluidRow(
          box(width = 6, title = "MCI distribution by stratum", status = "primary", solidHeader = TRUE,
              plotlyOutput("stratumPlot", height = "500px")),
          box(width = 6, title = "GTEx ratio by stratum", status = "warning", solidHeader = TRUE,
              plotlyOutput("ratioStratumPlot", height = "500px"))
        ),
        fluidRow(
          box(width = 12, title = "Stratum-level table", status = "primary", solidHeader = TRUE,
              DTOutput("stratumTable"))
        )
      ),
      
      tabItem(
        tabName = "table",
        fluidRow(
          box(width = 12, title = "Download complete concordance table", status = "primary", solidHeader = TRUE,
              downloadButton("downloadData", "Download CSV"),
              br(), br(),
              DTOutput("fullTable"))
        )
      ),
      
      tabItem(
        tabName = "about",
        fluidRow(
          box(width = 12, title = "About this browser", status = "primary", solidHeader = TRUE,
              HTML("
              <p><b>Molecular Concordance Index Browser</b></p>
              <p>This Shiny browser provides gene-level query access to Molecular Concordance Index results for ClinVar-annotated cardiomyopathy genes.</p>
              <p>The app displays primary MCI, GTEx-adjusted MCI, bootstrap 95% confidence intervals, tier assignments, GTEx low-confidence flags, per-cohort differential-expression statistics, stratum-level distributions, and downloadable concordance data.</p>
              <p>Interpretation: MCI measures transcript-level reproducibility across disease cohorts. GTEx adjustment adds a normal left-ventricle baseline comparison. Low-confidence GTEx flags do not mean a gene is biologically irrelevant; they mean disease-associated variability did not exceed normal GTEx left-ventricle variability in this bulk-tissue benchmark.</p>
              ")
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  selected_mci <- reactive({
    mci %>% filter(gene_symbol == input$gene)
  })
  
  selected_de <- reactive({
    de %>% filter(gene_symbol == input$gene)
  })
  
  output$mciBox <- renderValueBox({
    x <- selected_mci()
    valueBox(
      value = ifelse(nrow(x) == 0 || is.na(x$MCI_display[1]), "NA", x$MCI_display[1]),
      subtitle = "Primary MCI",
      icon = icon("dna"),
      color = "blue"
    )
  })
  
  output$adjBox <- renderValueBox({
    x <- selected_mci()
    valueBox(
      value = ifelse(nrow(x) == 0 || is.na(x$Adj_MCI_display[1]), "NA", x$Adj_MCI_display[1]),
      subtitle = "GTEx-adjusted MCI",
      icon = icon("chart-line"),
      color = "purple"
    )
  })
  
  output$tierBox <- renderValueBox({
    x <- selected_mci()
    tier <- ifelse(nrow(x) == 0, "NA", as.character(x$MCI_tier[1]))
    color <- ifelse(tier == "HIGH", "green", ifelse(tier == "MODERATE", "yellow", "red"))
    valueBox(
      value = tier,
      subtitle = "MCI Tier",
      icon = icon("layer-group"),
      color = color
    )
  })
  
  output$gtexBox <- renderValueBox({
    x <- selected_mci()
    flag <- ifelse(nrow(x) == 0, NA, x$GTEx_low_confidence_flag[1])
    status <- ifelse(nrow(x) == 0, NA, as.character(x$GTEx_adjustment_status[1]))
    
    label <- ifelse(
      isTRUE(flag) || flag == TRUE,
      "GTEx BASELINE",
      ifelse(grepl("EXCEEDS", status, ignore.case = TRUE), "EXCEEDS GTEx", "GTEx CHECK")
    )
    
    # Important: LOW CONF is not colored red because it is a separate evidence axis,
    # not the same thing as an UNSTABLE MCI tier.
    color <- ifelse(
      isTRUE(flag) || flag == TRUE,
      "yellow",
      ifelse(grepl("EXCEEDS", status, ignore.case = TRUE), "purple", "aqua")
    )
    
    valueBox(
      value = label,
      subtitle = "GTEx baseline status",
      icon = icon("scale-balanced"),
      color = color
    )
  })
  
  output$geneSummary <- renderDT({
    x <- selected_mci()
    keep <- intersect(c(
      "gene_symbol", "disease_group", "stratum", "MCI", "Adj_MCI",
      "MCI_CI95_lower", "MCI_CI95_upper", "MCI_tier",
      "bootstrap_majority_tier", "bootstrap_prob_HIGH",
      "bootstrap_prob_MODERATE", "bootstrap_prob_UNSTABLE",
      "sigma_disease", "sigma_GTEx", "sigma_disease_to_GTEx_ratio",
      "GTEx_low_confidence_flag", "GTEx_adjustment_status"
    ), names(x))
    
    out <- x[, keep, drop = FALSE]
    
    # Display formatting only. Raw CSV stays full precision.
    numeric_cols <- intersect(c(
      "MCI", "Adj_MCI", "MCI_CI95_lower", "MCI_CI95_upper",
      "bootstrap_prob_HIGH", "bootstrap_prob_MODERATE", "bootstrap_prob_UNSTABLE",
      "sigma_disease", "sigma_GTEx", "sigma_disease_to_GTEx_ratio"
    ), names(out))
    
    for (cc in numeric_cols) {
      out[[cc]] <- ifelse(is.na(out[[cc]]), NA, sprintf("%.3f", as.numeric(out[[cc]])))
    }
    
    datatable(out, options = list(scrollX = TRUE, pageLength = 5))
  })
  
  output$deTable <- renderDT({
    x <- selected_de()
    keep <- intersect(c("gene_symbol", "cohort", "log2FC", "SE", "p_value", "FDR", "direction", "stratum", "matched_in_DE"), names(x))
    datatable(x[, keep, drop = FALSE], options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$forestPlot <- renderPlotly({
    x <- selected_de()
    validate(need(nrow(x) > 0, "No per-cohort DE rows available for this gene."))
    
    if (!("SE" %in% names(x))) {
      x$SE <- NA_real_
    }
    
    x <- x %>%
      mutate(
        lower = ifelse(is.na(SE), log2FC, log2FC - 1.96 * SE),
        upper = ifelse(is.na(SE), log2FC, log2FC + 1.96 * SE),
        cohort = factor(cohort, levels = unique(cohort))
      )
    
    p <- ggplot(x, aes(y = cohort, x = log2FC, text = paste0(
      "Gene: ", gene_symbol,
      "<br>Cohort: ", cohort,
      "<br>log2FC: ", round(log2FC, 3),
      "<br>FDR: ", signif(FDR, 3),
      "<br>Direction: ", direction
    ))) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
      geom_point(aes(shape = direction), size = 3) +
      labs(x = "log2FC HCM vs Control", y = "Cohort") +
      theme_minimal(base_size = 13)
    
    ggplotly(p, tooltip = "text")
  })
  
  output$gtexPlot <- renderPlotly({
    x <- selected_mci()
    validate(need(nrow(x) > 0, "No MCI row available for this gene."))
    
    plot_df <- data.frame(
      metric = c("sigma_disease", "sigma_GTEx"),
      value = c(x$sigma_disease[1], x$sigma_GTEx[1])
    )
    
    p <- ggplot(plot_df, aes(x = metric, y = value, text = paste0(metric, ": ", round(value, 4)))) +
      geom_col() +
      labs(x = "", y = "Sigma") +
      theme_minimal(base_size = 13)
    
    ratio <- ifelse("sigma_disease_to_GTEx_ratio" %in% names(x), x$sigma_disease_to_GTEx_ratio[1], NA)
    p <- p + ggtitle(paste0("Disease/GTEx ratio = ", ifelse(is.na(ratio), "NA", round(ratio, 3))))
    
    ggplotly(p, tooltip = "text")
  })
  
  output$stratumPlot <- renderPlotly({
    x <- mci %>% filter(!is.na(stratum), !is.na(MCI))
    p <- ggplot(x, aes(x = stratum, y = MCI, text = paste0(
      "Gene: ", gene_symbol,
      "<br>Stratum: ", stratum,
      "<br>MCI: ", round(MCI, 3),
      "<br>Tier: ", MCI_tier
    ))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.8) +
      geom_hline(yintercept = 0.45, linetype = "dotted") +
      geom_hline(yintercept = 0.70, linetype = "dashed") +
      coord_flip() +
      labs(x = "Stratum", y = "MCI") +
      theme_minimal(base_size = 13)
    ggplotly(p, tooltip = "text")
  })
  
  output$ratioStratumPlot <- renderPlotly({
    x <- mci %>% filter(!is.na(stratum), !is.na(sigma_disease_to_GTEx_ratio))
    p <- ggplot(x, aes(x = stratum, y = sigma_disease_to_GTEx_ratio, text = paste0(
      "Gene: ", gene_symbol,
      "<br>Ratio: ", round(sigma_disease_to_GTEx_ratio, 3),
      "<br>GTEx flag: ", GTEx_low_confidence_flag
    ))) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.15, alpha = 0.8) +
      geom_hline(yintercept = 1.0, linetype = "dashed") +
      coord_flip() +
      labs(x = "Stratum", y = "sigma_disease / sigma_GTEx") +
      theme_minimal(base_size = 13)
    ggplotly(p, tooltip = "text")
  })
  
  output$stratumTable <- renderDT({
    x <- mci %>%
      filter(!is.na(stratum)) %>%
      group_by(stratum) %>%
      summarise(
        n_genes = n(),
        median_MCI = median(MCI, na.rm = TRUE),
        mean_MCI = mean(MCI, na.rm = TRUE),
        median_Adj_MCI = median(Adj_MCI, na.rm = TRUE),
        n_HIGH = sum(MCI_tier == "HIGH", na.rm = TRUE),
        n_MODERATE = sum(MCI_tier == "MODERATE", na.rm = TRUE),
        n_UNSTABLE = sum(MCI_tier == "UNSTABLE", na.rm = TRUE),
        n_GTEx_low_confidence = sum(GTEx_low_confidence_flag == TRUE, na.rm = TRUE),
        .groups = "drop"
      )
    datatable(x, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$fullTable <- renderDT({
    datatable(mci, options = list(scrollX = TRUE, pageLength = 15))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("MCI_concordance_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(mci, file)
    }
  )
}

shinyApp(ui = ui, server = server)
