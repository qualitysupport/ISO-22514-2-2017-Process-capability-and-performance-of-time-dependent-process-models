# ISO 22514-2:2017 — Process Capability & Performance Calculator
#
# Implements the procedure described in ISO 22514-2:2017(E):
#   Statistical methods in process management — Capability and
#   performance — Part 2: Process capability and performance of
#   time-dependent process models.
#
# Covers:
#   Clause 5     - Time-dependent distribution models (A1, A2, B, C1-C4, D)
#   Clause 6.1   - Overview of performance/capability index formulas
#                  (Pp / Ppk, Cp / Cpk, one- and two-sided)
#   Clause 6.1.2 - Location calculation methods (Ml, l = 1-4)
#   Clause 6.1.3 - Dispersion calculation methods (Md, d = 1-5, incl.
#                  c4 / d2 bias-correction constants)
#   Clause 6.1.4 - Calculation of X0.135% and X99.865% (empirical
#                  percentile, fitted normal, or min/max)
#   Clause 6.2   - One-sided specification limits
#   Clause 6.3   - Guidance on which calculation methods suit which
#                  time-dependent model
#   Clause 7     - Reporting of process performance/capability indices
#
# This app implements the standard's calculation FORMULAS in original
# code; it does not reproduce any text or figures from the standard
# itself.
#
# Required packages: shiny, ggplot2, DT
#
# Author: Dan Lay Jr. | Calibration Support LLC
# www.calibrationsupport.com | linkedin.com/in/dlayjr
#
# Disclaimer: This tool implements calculation methods described in
# ISO 22514-2:2017 but is not affiliated with, endorsed by, or reviewed
# by ISO. It does not reproduce the standard itself; consult the official
# ISO 22514-2:2017 document (available for purchase from ISO or your
# national standards body) for authoritative guidance. Provided as-is,
# for educational and reference use.
########################################################################

library(shiny)
library(ggplot2)
library(DT)

# ---------------------------------------------------------------------------
# Reference constants (standard SPC bias-correction constants for subgroup
# statistics, used by dispersion methods d = 3 and d = 4)
# ---------------------------------------------------------------------------
c4_fun <- function(n) {
  if (n < 2) return(NA_real_)
  sqrt(2 / (n - 1)) * exp(lgamma(n / 2) - lgamma((n - 1) / 2))
}

d2_table <- c(
  `2` = 1.128, `3` = 1.693, `4` = 2.059, `5` = 2.326, `6` = 2.534,
  `7` = 2.704, `8` = 2.847, `9` = 2.970, `10` = 3.078, `11` = 3.173,
  `12` = 3.258, `13` = 3.336, `14` = 3.407, `15` = 3.472, `16` = 3.532,
  `17` = 3.588, `18` = 3.640, `19` = 3.689, `20` = 3.735, `21` = 3.778,
  `22` = 3.819, `23` = 3.858, `24` = 3.895, `25` = 3.931
)
d2_fun <- function(n) {
  key <- as.character(n)
  if (n >= 2 && n <= 25 && key %in% names(d2_table)) return(unname(d2_table[key]))
  NA_real_
}

# ---------------------------------------------------------------------------
# Model metadata (own wording — short labels + example use case only)
# ---------------------------------------------------------------------------
model_choices <- c(
  "A1 - constant location & dispersion, normal" = "A1",
  "A2 - constant location & dispersion, non-normal unimodal" = "A2",
  "B  - constant location, varying dispersion" = "B",
  "C1 - varying location (normal), constant dispersion" = "C1",
  "C2 - varying location (non-normal), constant dispersion" = "C2",
  "C3 - trending location, constant dispersion" = "C3",
  "C4 - shifting/batch location, constant dispersion" = "C4",
  "D  - both location & dispersion vary" = "D"
)

model_notes <- list(
  A1 = "Location and dispersion both stable; short-term and overall distributions are both normal. This is the classical 'process in statistical control' case.",
  A2 = "Location and dispersion both stable, but the underlying distribution is naturally not normal (e.g. a physically bounded characteristic).",
  B  = "Average level stays put, but spread changes over time (systematically or randomly) — e.g. gradual tool wear affecting spread across spindles.",
  C1 = "Spread stays put, but the average drifts randomly and is itself normally distributed from interval to interval.",
  C2 = "Spread stays put, but the average drifts randomly with a non-normal pattern.",
  C3 = "Spread stays put, but the average follows a systematic function of time (trend, cyclical wear, etc.).",
  C4 = "Spread stays put, but the average jumps systematically and randomly (e.g. tool changes, batch/lot changes).",
  D  = "Both average and spread change over time (systematic and/or random) — the general, least-controlled case."
)

# Soft guidance (paraphrased from the standard's explicit notes in 6.1.4):
# dispersion methods 2-4 estimate only within-subgroup spread and assume
# normality, so they are most defensible for model A1.
dispersion_guidance <- function(model, disp_method) {
  msgs <- character(0)
  if (disp_method %in% c(2, 3, 4) && model != "A1") {
    msgs <- c(msgs, "This dispersion method estimates only the within-subgroup spread and neglects between-subgroup variation - it is best suited to model A1 (stable process). For other models, consider dispersion method 1 or 5.")
  }
  if (disp_method %in% c(2, 3, 4, 5) && model %in% c("A2", "B", "C2", "C3", "C4", "D")) {
    msgs <- c(msgs, "This dispersion method assumes an underlying normal distribution. Because the selected model implies a non-normal or unstable resulting distribution, the estimate may be biased. Dispersion method 1 (quantile-based) is more general.")
  }
  msgs
}

# ---------------------------------------------------------------------------
# Core calculation engine
# ---------------------------------------------------------------------------

make_subgroups <- function(x, n) {
  if (n <= 1) return(as.list(x))
  k <- floor(length(x) / n)
  if (k < 1) return(list())
  lapply(seq_len(k), function(i) x[((i - 1) * n + 1):(i * n)])
}

calc_location <- function(x, subgroups, method) {
  switch(as.character(method),
         "1" = mean(x, na.rm = TRUE),
         "2" = median(x, na.rm = TRUE),
         "3" = mean(sapply(subgroups, mean, na.rm = TRUE)),
         "4" = mean(sapply(subgroups, median, na.rm = TRUE)),
         stop("Unknown location method")
  )
}

calc_dispersion <- function(x, subgroups, method, n, pctl_method = "empirical") {
  N <- length(x)
  
  if (method == 1) {
    # Quantile-based: X0.135% and X99.865% from the resulting distribution
    q <- get_quantiles(x, pctl_method)
    return(list(sigma = NA_real_, q_low = q$q_low, q_high = q$q_high, kind = "quantile"))
  }
  
  if (method == 2) {
    s2 <- sapply(subgroups, function(g) var(g, na.rm = TRUE))
    sigma <- sqrt(mean(s2, na.rm = TRUE))
  } else if (method == 3) {
    s <- sapply(subgroups, function(g) sd(g, na.rm = TRUE))
    sigma <- mean(s, na.rm = TRUE) / c4_fun(n)
  } else if (method == 4) {
    R <- sapply(subgroups, function(g) diff(range(g, na.rm = TRUE)))
    sigma <- mean(R, na.rm = TRUE) / d2_fun(n)
  } else if (method == 5) {
    sigma <- sd(x, na.rm = TRUE)
  } else {
    stop("Unknown dispersion method")
  }
  list(sigma = sigma, kind = "sigma")
}

get_quantiles <- function(x, pctl_method) {
  N <- length(x)
  if (pctl_method == "normal_fit") {
    mu <- mean(x, na.rm = TRUE); s <- sd(x, na.rm = TRUE)
    list(q_low = qnorm(0.00135, mu, s), q_high = qnorm(0.99865, mu, s))
  } else if (pctl_method == "minmax") {
    list(q_low = min(x, na.rm = TRUE), q_high = max(x, na.rm = TRUE))
  } else {
    list(q_low = as.numeric(quantile(x, 0.00135, na.rm = TRUE, type = 6)),
         q_high = as.numeric(quantile(x, 0.99865, na.rm = TRUE, type = 6)))
  }
}

# Full index computation, given data + all chosen settings
compute_indices <- function(x, n, loc_method, disp_method, pctl_method,
                            L, U, side = "two", statistical_control = FALSE) {
  
  subgroups <- make_subgroups(x, n)
  Xmid <- calc_location(x, subgroups, loc_method)
  disp <- calc_dispersion(x, subgroups, disp_method, n, pctl_method)
  
  if (disp$kind == "quantile") {
    deltaL <- Xmid - disp$q_low
    deltaU <- disp$q_high - Xmid
    delta  <- disp$q_high - disp$q_low
    q_low <- disp$q_low; q_high <- disp$q_high
  } else {
    sigma <- disp$sigma
    delta <- 6 * sigma
    deltaL <- 3 * sigma
    deltaU <- 3 * sigma
    q_low <- Xmid - deltaL
    q_high <- Xmid + deltaU
  }
  
  res <- list(
    Xmid = Xmid, delta = delta, deltaL = deltaL, deltaU = deltaU,
    X0135 = q_low, X99865 = q_high, subgroups = subgroups,
    n_subgroups = length(subgroups)
  )
  
  idx_prefix <- if (statistical_control) "C" else "P"
  
  if (side == "two") {
    res$index_p  <- (U - L) / delta
    res$index_pkL <- (Xmid - L) / deltaL
    res$index_pkU <- (U - Xmid) / deltaU
    res$index_pk  <- min(res$index_pkL, res$index_pkU)
  } else if (side == "upper") {
    res$index_pkU <- (U - Xmid) / deltaU
    res$index_pk  <- res$index_pkU
  } else if (side == "lower") {
    res$index_pkL <- (Xmid - L) / deltaL
    res$index_pk  <- res$index_pkL
  }
  res$idx_prefix <- idx_prefix
  res
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("ISO 22514-2 Process Capability / Performance Calculator"),
  tags$head(tags$style(HTML("
    .well { background-color: #f7f9fb; }
    .index-box { padding: 14px; border-radius: 8px; background: #eef4fb; margin-bottom: 10px; }
    .index-value { font-size: 26px; font-weight: 700; }
    .index-label { font-size: 13px; color: #555; }
  "))),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Data"),
      radioButtons("data_source", NULL,
                   choices = c("Upload CSV (one numeric column, time order)" = "upload",
                               "Paste values" = "paste",
                               "Load example data" = "example"),
                   selected = "example"),
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("file", "CSV file", accept = ".csv"),
        checkboxInput("header", "File has header row", TRUE)
      ),
      conditionalPanel(
        condition = "input.data_source == 'paste'",
        textAreaInput("pasted", "One value per line or comma-separated",
                      rows = 6, placeholder = "20.01\n19.98\n20.05\n...")
      ),
      conditionalPanel(
        condition = "input.data_source == 'example'",
        selectInput("example_model", "Example model to simulate", choices = model_choices),
        actionButton("gen_example", "Generate example data", class = "btn-sm")
      ),
      
      numericInput("subgroup_size", "Subgroup size n (use 1 for individual values)",
                   value = 5, min = 1, step = 1),
      
      hr(),
      h4("2. Time-dependent distribution model"),
      selectInput("model", "Select model (ISO 22514-2, Table 1)", choices = model_choices),
      checkboxInput("in_control", "Process is in statistical control (report Cp/Cpk instead of Pp/Ppk)", FALSE),
      helpText(textOutput("model_note")),
      
      hr(),
      h4("3. Calculation method"),
      selectInput("loc_method", "Location method (l)",
                  choices = c(
                    "1 - mean of all individual values" = 1,
                    "2 - median of combined data (X50%)" = 2,
                    "3 - mean of subgroup means" = 3,
                    "4 - mean of subgroup medians" = 4
                  ), selected = 2),
      selectInput("disp_method", "Dispersion method (d)",
                  choices = c(
                    "1 - quantile-based (X0.135% / X99.865%)" = 1,
                    "2 - pooled subgroup variance" = 2,
                    "3 - mean subgroup s / c4" = 3,
                    "4 - mean subgroup range / d2" = 4,
                    "5 - total sample standard deviation" = 5
                  ), selected = 1),
      conditionalPanel(
        condition = "input.disp_method == 1",
        selectInput("pctl_method", "How to obtain X0.135% / X99.865%",
                    choices = c("Empirical percentile from data" = "empirical",
                                "Fit normal distribution" = "normal_fit",
                                "Min / max (large samples, N>=1000)" = "minmax"))
      ),
      uiOutput("method_warning"),
      
      hr(),
      h4("4. Specification limits"),
      radioButtons("side", "Specification type",
                   choices = c("Two-sided (L and U)" = "two",
                               "One-sided - upper only" = "upper",
                               "One-sided - lower only" = "lower"),
                   selected = "two"),
      conditionalPanel(condition = "input.side != 'upper'", numericInput("L", "Lower specification limit (L)", value = 19.8, step = 0.01)),
      conditionalPanel(condition = "input.side != 'lower'", numericInput("U", "Upper specification limit (U)", value = 20.2, step = 0.01)),
      
      hr(),
      actionButton("calc", "Calculate", class = "btn-primary btn-lg", width = "100%")
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Results",
                 br(),
                 uiOutput("index_boxes"),
                 h4("Full report (ISO 22514-2, Clause 7 style)"),
                 DT::DTOutput("report_table"),
                 br(),
                 downloadButton("download_report", "Download report (CSV)")
        ),
        tabPanel("Run chart", br(), plotOutput("run_chart", height = "450px")),
        tabPanel("Histogram", br(), plotOutput("histogram", height = "450px")),
        tabPanel("Data", br(), DT::DTOutput("data_table")),
        tabPanel("About the models",
                 br(),
                 h4("Time-dependent distribution models (own summary)"),
                 p("ISO 22514-2 classifies processes into eight models based on whether location and dispersion are constant, and whether any change is random or systematic. Select a model in the sidebar to see its description; the app does not enforce a specific method for a model but will flag calculation methods that the standard's own notes caution against for that situation."),
                 tags$ul(
                   lapply(names(model_notes), function(m) {
                     tags$li(strong(m), ": ", model_notes[[m]])
                   })
                 )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Example data generator -------------------------------------------
  example_data <- eventReactive(input$gen_example, {
    set.seed(42)
    m <- input$example_model
    n_total <- 300
    if (m == "A1") {
      x <- rnorm(n_total, 20, 0.05)
    } else if (m == "A2") {
      x <- rgamma(n_total, shape = 2, scale = 0.05) + 19.9
    } else if (m == "B") {
      sigma_t <- seq(0.03, 0.12, length.out = n_total)
      x <- rnorm(n_total, 20, sigma_t)
    } else if (m == "C1") {
      mu_t <- rnorm(n_total / 5, 20, 0.06)
      x <- unlist(lapply(mu_t, function(mu) rnorm(5, mu, 0.03)))
    } else if (m == "C2") {
      mu_t <- rgamma(n_total / 5, shape = 2, scale = 0.05) + 19.85
      x <- unlist(lapply(mu_t, function(mu) rnorm(5, mu, 0.03)))
    } else if (m == "C3") {
      trend <- seq(20.15, 19.9, length.out = n_total)
      x <- rnorm(n_total, trend, 0.03)
    } else if (m == "C4") {
      shifts <- rep(c(20.0, 20.06, 19.95, 20.03), each = n_total / 4)
      x <- rnorm(n_total, shifts, 0.02)
    } else { # D
      shifts <- rep(c(20.0, 20.08, 19.92), each = n_total / 3)
      sigma_t <- rep(c(0.03, 0.06, 0.1), each = n_total / 3)
      x <- rnorm(n_total, shifts, sigma_t)
    }
    round(x, 4)
  }, ignoreNULL = FALSE)
  
  raw_data <- reactive({
    src <- input$data_source
    if (src == "upload") {
      shiny::req(input$file)
      df <- read.csv(input$file$datapath, header = input$header)
      as.numeric(df[[1]])
    } else if (src == "paste") {
      txt <- input$pasted
      shiny::req(nzchar(txt))
      vals <- unlist(strsplit(txt, "[,\\n\\r\\t ]+"))
      vals <- vals[nzchar(vals)]
      as.numeric(vals)
    } else {
      example_data()
    }
  })
  
  data_vec <- reactive({
    x <- raw_data()
    x <- x[!is.na(x)]
    shiny::validate(shiny::need(length(x) >= max(4, input$subgroup_size * 2),
                                "Need at least a few subgroups' worth of numeric data. Please provide more values."))
    x
  })
  
  output$data_table <- DT::renderDT({
    DT::datatable(data.frame(index = seq_along(data_vec()), value = data_vec()),
                  options = list(pageLength = 15))
  })
  
  output$model_note <- renderText({
    model_notes[[input$model]]
  })
  
  output$method_warning <- renderUI({
    msgs <- dispersion_guidance(input$model, as.integer(input$disp_method))
    if (length(msgs) == 0) return(NULL)
    tags$div(style = "color:#a15c00; background:#fff6e5; border:1px solid #f0d38a; padding:8px; border-radius:6px; margin-top:6px; font-size:12.5px;",
             lapply(msgs, function(m) tags$p(m, style="margin:2px 0;")))
  })
  
  results <- eventReactive(input$calc, {
    x <- data_vec()
    n <- max(1, as.integer(input$subgroup_size))
    disp_m <- as.integer(input$disp_method)
    
    if (input$side == "two") {
      shiny::req(input$L, input$U)
      shiny::validate(shiny::need(input$U > input$L, "Upper specification limit must be greater than lower."))
    }
    if (disp_m %in% c(2, 3, 4)) {
      shiny::validate(shiny::need(n >= 2, "Dispersion methods 2, 3 and 4 require a subgroup size n >= 2 (they use within-subgroup variance/range). Use n = 1 with dispersion method 1 or 5 instead."))
    }
    if (disp_m == 4) {
      shiny::validate(shiny::need(n <= 25, "Dispersion method 4 (range/d2) is only tabulated here for subgroup sizes up to 25. Use method 2, 3, 1 or 5 for larger subgroups."))
    }
    
    pctl_m <- if (is.null(input$pctl_method)) "empirical" else input$pctl_method
    
    compute_indices(
      x = x, n = n,
      loc_method = as.integer(input$loc_method),
      disp_method = disp_m,
      pctl_method = pctl_m,
      L = if (!is.null(input$L)) input$L else NA,
      U = if (!is.null(input$U)) input$U else NA,
      side = input$side,
      statistical_control = input$in_control
    )
  })
  
  output$index_boxes <- renderUI({
    r <- tryCatch(results(), error = function(e) NULL)
    if (is.null(r)) return(p("Set your options and click Calculate."))
    
    prefix <- r$idx_prefix
    boxes <- list()
    
    if (input$side == "two") {
      boxes <- list(
        tags$div(class = "index-box", tags$div(class = "index-value", sprintf("%.3f", r$index_p)),
                 tags$div(class = "index-label", paste0(prefix, "p"))),
        tags$div(class = "index-box", tags$div(class = "index-value", sprintf("%.3f", r$index_pk)),
                 tags$div(class = "index-label", paste0(prefix, "pk (min)"))),
        tags$div(class = "index-box", tags$div(class = "index-value", sprintf("%.3f", r$index_pkL)),
                 tags$div(class = "index-label", paste0(prefix, "pk-L"))),
        tags$div(class = "index-box", tags$div(class = "index-value", sprintf("%.3f", r$index_pkU)),
                 tags$div(class = "index-label", paste0(prefix, "pk-U")))
      )
    } else {
      boxes <- list(
        tags$div(class = "index-box", tags$div(class = "index-value", sprintf("%.3f", r$index_pk)),
                 tags$div(class = "index-label", paste0(prefix, "pk (", input$side, ")")))
      )
    }
    fluidRow(lapply(boxes, function(b) column(3, b)))
  })
  
  output$report_table <- DT::renderDT({
    r <- tryCatch(results(), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    prefix <- r$idx_prefix
    rows <- list(
      c("Process location (Xmid)", sprintf("%.4f", r$Xmid)),
      c("X0.135%", sprintf("%.4f", r$X0135)),
      c("X99.865%", sprintf("%.4f", r$X99865)),
      c("Dispersion Delta", sprintf("%.4f", r$delta)),
      c("Delta L", sprintf("%.4f", r$deltaL)),
      c("Delta U", sprintf("%.4f", r$deltaU))
    )
    if (input$side == "two") {
      rows <- c(rows, list(
        c(paste0(prefix, "p"), sprintf("%.3f", r$index_p)),
        c(paste0(prefix, "pk-L"), sprintf("%.3f", r$index_pkL)),
        c(paste0(prefix, "pk-U"), sprintf("%.3f", r$index_pkU)),
        c(paste0(prefix, "pk (minimum)"), sprintf("%.3f", r$index_pk))
      ))
    } else {
      rows <- c(rows, list(c(paste0(prefix, "pk"), sprintf("%.3f", r$index_pk))))
    }
    rows <- c(rows, list(
      c("Calculation method", paste0("M", input$loc_method, ",", input$disp_method)),
      c("Time distribution model", input$model),
      c("Number of values used", length(data_vec())),
      c("Number of subgroups", r$n_subgroups),
      c("Subgroup size (n)", input$subgroup_size),
      c("Statistical control assumed", ifelse(input$in_control, "Yes (Capability, C)", "No (Performance, P)"))
    ))
    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(df) <- c("Item", "Value")
    DT::datatable(df, options = list(dom = "t", paging = FALSE), rownames = FALSE)
  })
  
  output$download_report <- downloadHandler(
    filename = function() paste0("iso22514-2_report_", Sys.Date(), ".csv"),
    content = function(file) {
      r <- results()
      prefix <- r$idx_prefix
      rows <- list(
        c("Process location (Xmid)", r$Xmid),
        c("X0.135%", r$X0135),
        c("X99.865%", r$X99865),
        c("Dispersion Delta", r$delta),
        c("Delta L", r$deltaL),
        c("Delta U", r$deltaU),
        c("Calculation method", paste0("M", input$loc_method, ",", input$disp_method)),
        c("Time distribution model", input$model),
        c("Number of values used", length(data_vec())),
        c("Number of subgroups", r$n_subgroups),
        c("Subgroup size (n)", input$subgroup_size)
      )
      if (input$side == "two") {
        rows <- c(rows, list(
          c(paste0(prefix, "p"), r$index_p),
          c(paste0(prefix, "pk-L"), r$index_pkL),
          c(paste0(prefix, "pk-U"), r$index_pkU),
          c(paste0(prefix, "pk"), r$index_pk)
        ))
      } else {
        rows <- c(rows, list(c(paste0(prefix, "pk"), r$index_pk)))
      }
      df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
      names(df) <- c("Item", "Value")
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$run_chart <- renderPlot({
    r <- tryCatch(results(), error = function(e) NULL)
    x <- data_vec()
    df <- data.frame(idx = seq_along(x), value = x)
    p <- ggplot(df, aes(idx, value)) +
      geom_line(color = "#2f7ed8", alpha = 0.5) +
      geom_point(color = "#1a4f8c", size = 1.3) +
      labs(x = "Observation number", y = "Characteristic value", title = "Run chart") +
      theme_minimal(base_size = 13)
    
    if (!is.null(r)) {
      p <- p + geom_hline(yintercept = r$Xmid, color = "black", linetype = "dashed") +
        geom_hline(yintercept = r$X0135, color = "#3a9d3a", linetype = "dotted") +
        geom_hline(yintercept = r$X99865, color = "#3a9d3a", linetype = "dotted")
      if (!is.null(input$L) && input$side != "upper") p <- p + geom_hline(yintercept = input$L, color = "red")
      if (!is.null(input$U) && input$side != "lower") p <- p + geom_hline(yintercept = input$U, color = "red")
    }
    p
  })
  
  output$histogram <- renderPlot({
    r <- tryCatch(results(), error = function(e) NULL)
    x <- data_vec()
    df <- data.frame(value = x)
    p <- ggplot(df, aes(value)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#6fb85e", color = "white", alpha = 0.85) +
      geom_density(color = "#2f7ed8", linewidth = 1) +
      labs(x = "Characteristic value", y = "Density", title = "Histogram of resulting distribution") +
      theme_minimal(base_size = 13)
    
    if (!is.null(r)) {
      p <- p + geom_vline(xintercept = r$Xmid, color = "black", linetype = "dashed") +
        geom_vline(xintercept = r$X0135, color = "#3a9d3a", linetype = "dotted") +
        geom_vline(xintercept = r$X99865, color = "#3a9d3a", linetype = "dotted")
      if (!is.null(input$L) && input$side != "upper") p <- p + geom_vline(xintercept = input$L, color = "red")
      if (!is.null(input$U) && input$side != "lower") p <- p + geom_vline(xintercept = input$U, color = "red")
    }
    p
  })
}

shinyApp(ui, server)
