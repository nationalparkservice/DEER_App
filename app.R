# Deer Density Lab – USCR vs REM vs TTE
# Models run via run_USCR(), run_REM(), run_TTE() in R/sim_and_models.R

# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------

needs <- c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr", "stringr", "secr", "data.table",
  "leaflet", "ggrepel", "shinyjs",
  "nimble", "parallel", "MCMCvis", "lubridate", "sf", "tibble",
  "future", "promises"
)

# Redwood-inspired palette (hex only — no extra color package)
redwood_colors <- c(
  "#303018", "#604830", "#609048", "#90A860", "#786048", "#B8C4A8"
)

missing <- needs[!needs %in% rownames(installed.packages())]
if (length(missing)) {
  message("DEER App: installing missing packages: ", paste(missing, collapse = ", "))
  utils::install.packages(missing, dependencies = TRUE)
}
failed <- character()
for (pkg in needs) {
  if (!requireNamespace(pkg, quietly = TRUE)) failed <- c(failed, pkg)
}
if (length(failed)) {
  stop(
    "DEER App could not load: ", paste(failed, collapse = ", "),
    "\nInstall manually, e.g. install.packages(c(",
    paste0("\"", failed, "\"", collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(secr)
  library(data.table)
  library(leaflet)
  library(ggrepel)
  library(tibble)
  library(sf)
  library(nimble)
  library(shinyjs)
  library(future)
  library(promises)
})

# Conservative async worker pool for server-side background REM jobs.
# Keep this small because the model code can also use parallel chains internally.
available_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.na(available_cores) || available_cores < 2) {
  async_workers <- 1L
} else {
  async_workers <- min(2L, available_cores - 1L)
}
future::plan(future::multisession, workers = async_workers)

# -------------------------------------------------------------------
# Helper files
# -------------------------------------------------------------------

source("R/sim_and_models.R")   # run_USCR(), run_REM(), run_TTE()
source("R/data_checks.R")      # QC + summary helpers
source("R/camera_array_helpers.R")

# Pin the active USCR implementation to the current file contents so stale
# objects from prior sessions cannot route the app through an older wrapper.
uscr_model_env <- new.env(parent = globalenv())
sys.source("R/sim_and_models.R", envir = uscr_model_env)
run_USCR_app <- uscr_model_env$run_USCR

# -------------------------------------------------------------------
# Simulation helpers
# -------------------------------------------------------------------
# Simulation and model-input helpers now come directly from R/sim_and_models.R:
#   simulate_camera_counts(), get_counts_matrix(), capthist_to_events(),
#   sim_model_inputs(), build_nps_model_inputs()

# -------------------------------------------------------------------
# Generic helpers for summaries & WAIC-combo (like Rmd)
# -------------------------------------------------------------------

get_waic_value <- function(fit) {
  if (is.null(fit$waic)) return(NA_real_)
  if (is.list(fit$waic) && !is.null(fit$waic$WAIC)) {
    as.numeric(fit$waic$WAIC)
  } else {
    as.numeric(fit$waic)
  }
}

summarize_method <- function(fit) {
  D_mi2 <- fit$samples_all[, "D_mi2"]
  D_km2 <- D_mi2 / 2.59
  
  list(
    mean_km2  = mean(D_km2),
    sd_km2    = sd(D_km2),
    q2.5_km2  = quantile(D_km2, 0.025),
    q97.5_km2 = quantile(D_km2, 0.975),
    waic      = get_waic_value(fit)
  )
}

extract_density_draws_km2 <- function(fit) {
  if (is.null(fit) || is.null(fit$samples_all) || !"D_mi2" %in% colnames(fit$samples_all)) {
    return(numeric(0))
  }
  x <- as.numeric(fit$samples_all[, "D_mi2"]) / 2.59
  x[is.finite(x)]
}

prepare_combo_density_data <- function(fits) {
  model_order <- c("REM", "TTE", "USCR")
  fits <- fits[model_order[model_order %in% names(fits)]]
  fits <- fits[!vapply(fits, is.null, logical(1))]
  if (!length(fits)) return(NULL)

  model_names <- names(fits)
  waic_tbl <- tibble::tibble(
    model = model_names,
    waic = vapply(fits, get_waic_value, numeric(1))
  )

  waic_weights_available <- FALSE
  if (length(model_names) == 1L) {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(deltaWAIC = 0, rel_lik = 1, w = 1)
    waic_weights_available <- TRUE
  } else if (any(!is.finite(waic_tbl$waic))) {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(
        deltaWAIC = NA_real_,
        rel_lik = NA_real_,
        w = NA_real_
      )
  } else {
    waic_tbl <- waic_tbl %>%
      dplyr::mutate(
        deltaWAIC = waic - min(waic),
        rel_lik = exp(-0.5 * deltaWAIC),
        w = rel_lik / sum(rel_lik)
      )
    waic_weights_available <- TRUE
  }

  density_draws <- lapply(fits, function(fit) sort(extract_density_draws_km2(fit)))
  n_draws <- min(lengths(density_draws))
  if (!is.finite(n_draws) || n_draws < 1L) return(NULL)

  mat <- do.call(
    cbind,
    lapply(density_draws, function(x) x[seq_len(n_draws)])
  )
  colnames(mat) <- model_names

  density_est <- mat
  summary_rows <- character()
  if (ncol(mat) > 1L) {
    density_est <- cbind(mat, `unweighted mean` = rowMeans(mat))
    summary_rows <- "unweighted mean"
    if (isTRUE(waic_weights_available)) {
      density_est <- cbind(
        density_est,
        `weighted mean` = as.numeric(mat %*% waic_tbl$w)
      )
      summary_rows <- c(summary_rows, "weighted mean")
    }
  }

  list(
    fits = fits,
    model_names = model_names,
    waic = waic_tbl,
    waic_weights_available = waic_weights_available,
    model_draws = density_draws,
    density_est = density_est,
    summary_rows = summary_rows
  )
}

build_combo_table_from_fits <- function(fits, density_threshold = 20) {
  combo_data <- prepare_combo_density_data(fits)
  if (is.null(combo_data)) return(NULL)

  waic_tbl <- combo_data$waic
  model_names <- combo_data$model_names
  density_est <- combo_data$density_est
  summary_rows <- combo_data$summary_rows
  waic_weights_available <- combo_data$waic_weights_available

  means <- apply(density_est, 2, mean)
  lower <- apply(density_est, 2, stats::quantile, probs = 0.025)
  upper <- apply(density_est, 2, stats::quantile, probs = 0.975)
  prob_threshold <- apply(density_est, 2, function(x) mean(x > density_threshold))
  threshold_label <- paste0(
    "P(D > ",
    format(density_threshold, trim = TRUE, scientific = FALSE),
    " animals/km²)"
  )

  table_out <- tibble::tibble(
    Method = c(model_names, summary_rows),
    `Mean density (animals/km²)` = as.numeric(means),
    `Lower 2.5%` = as.numeric(lower),
    `Upper 97.5%` = as.numeric(upper),
    `WAIC weight` = c(waic_tbl$w, rep(NA_real_, length(summary_rows)))
  )
  table_out[[threshold_label]] <- as.numeric(prob_threshold)
  table_out <- table_out[, c(
    "Method",
    "Mean density (animals/km²)",
    "Lower 2.5%",
    "Upper 97.5%",
    threshold_label,
    "WAIC weight"
  )]

  list(
    table = table_out,
    waic = waic_tbl,
    models_available = model_names,
    waic_weights_available = waic_weights_available
  )
}

combo_plot_colors <- c(
  "REM" = "#8C4F2D",
  "TTE" = "#3F7F5F",
  "USCR" = "#4F6FA3",
  "unweighted mean" = "#8A8677",
  "weighted mean" = "#2F3E46"
)

build_combo_interval_plot <- function(fits, title_text) {
  combo_data <- prepare_combo_density_data(fits)
  if (is.null(combo_data)) return(NULL)

  interval_df <- tibble::tibble(
    Method = colnames(combo_data$density_est),
    mean = as.numeric(apply(combo_data$density_est, 2, mean)),
    lower = as.numeric(apply(combo_data$density_est, 2, stats::quantile, probs = 0.025)),
    upper = as.numeric(apply(combo_data$density_est, 2, stats::quantile, probs = 0.975))
  ) %>%
    dplyr::mutate(
      Summary = dplyr::case_when(
        Method %in% combo_data$model_names ~ "Model",
        Method == "weighted mean" ~ "Weighted mean",
        TRUE ~ "Unweighted mean"
      ),
      Method = factor(Method, levels = rev(Method)),
      point_shape = dplyr::case_when(
        Summary == "Model" ~ 16,
        Summary == "Weighted mean" ~ 15,
        TRUE ~ 17
      )
    )

  interval_levels <- rev(colnames(combo_data$density_est))
  interval_colors <- combo_plot_colors[interval_levels]
  interval_colors <- interval_colors[!is.na(interval_colors)]

  ggplot(interval_df, aes(x = mean, y = Method, color = Method)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 1, alpha = 0.95) +
    geom_point(aes(shape = point_shape), size = 3.2, show.legend = FALSE) +
    scale_color_manual(
      values = interval_colors,
      drop = FALSE
    ) +
    scale_shape_identity() +
    labs(
      title = title_text,
      x = "Density (animals/km²)",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

build_combo_density_overlay_plot <- function(fits, title_text) {
  combo_data <- prepare_combo_density_data(fits)
  if (is.null(combo_data)) return(NULL)

  density_df <- purrr::imap_dfr(
    combo_data$model_draws,
    ~ tibble::tibble(
      Method = .y,
      density_km2 = .x
    )
  )
  if (!nrow(density_df)) return(NULL)

  density_colors <- combo_plot_colors[c("REM", "TTE", "USCR")]

  ggplot(density_df, aes(x = density_km2, color = Method, fill = Method)) +
    geom_density(alpha = 0.16, linewidth = 1) +
    scale_color_manual(
      values = density_colors,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = density_colors,
      drop = FALSE
    ) +
    labs(
      title = title_text,
      x = "Density (animals/km²)",
      y = "Posterior density",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

#' Simulated data: only USCR is fit — single-model summary table (no WAIC averaging)
build_sim_combo_table_uscr_only <- function(uscr_fit, density_threshold = 20) {
  if (is.null(uscr_fit) || is.null(uscr_fit$samples_all)) return(NULL)
  D <- uscr_fit$samples_all[, "D_mi2"] / 2.59
  waic_val <- get_waic_value(uscr_fit)
  waic_tbl <- tibble::tibble(
    model     = "USCR",
    waic      = waic_val,
    deltaWAIC = 0,
    rel_lik   = 1,
    w         = 1
  )
  threshold_label <- paste0(
    "P(D > ",
    format(density_threshold, trim = TRUE, scientific = FALSE),
    " animals/km²)"
  )
  table_out <- tibble::tibble(
    Method                    = "USCR",
    `Mean density (animals/km²)` = mean(D),
    `Lower 2.5%`              = stats::quantile(D, 0.025),
    `Upper 97.5%`             = stats::quantile(D, 0.975),
    `WAIC weight`             = 1
  )
  table_out[[threshold_label]] <- mean(D > density_threshold)
  table_out <- table_out[, c(
    "Method",
    "Mean density (animals/km²)",
    "Lower 2.5%",
    "Upper 97.5%",
    threshold_label,
    "WAIC weight"
  )]
  list(table = table_out, waic = waic_tbl)
}

#' Posterior summary for CSV export (parameter names, mean, 95% CI)
posterior_summary_df <- function(fit) {
  if (is.null(fit) || is.null(fit$samples_all)) return(NULL)
  m <- as.data.frame(fit$samples_all)
  if ("D_mi2" %in% names(m)) {
    m$D_km2 <- m$D_mi2 / 2.59
    m$D_mi2 <- NULL
  }
  tibble::tibble(
    parameter = names(m),
    mean      = sapply(m, mean),
    q025      = sapply(m, function(x) stats::quantile(x, 0.025)),
    q975      = sapply(m, function(x) stats::quantile(x, 0.975))
  )
}

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------

ui <- page_fillable(
  shinyjs::useShinyjs(),
  theme = bs_theme(
    version = 5, 
    bootswatch = "flatly",
    primary = redwood_colors[3],  # #609048 - green
    secondary = redwood_colors[2], # #604830 - medium brown
    success = redwood_colors[4],   # #90A860 - light green
    bg = "#fafafa",
    fg = redwood_colors[1]         # #303018 - dark brown/green
  ),
  tags$head(
    tags$script(src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"),
    tags$script(HTML("
      window.MathJax = { 
        tex: { 
          inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
          displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
          processEscapes: true,
          processEnvironments: true
        } 
      };
    ")),
    tags$style(HTML("
      :root {
        --rw1: #A45A52;
        --rw2: #52372E;
        --rw3: #8C9B7A;
        --rw4: #EDE8E2;
        --rw5: #3C4B3E;
        --ink: #213026;
        --muted: #5c675d;
        --radius: 14px;
        --shadow: 0 10px 26px rgba(33,48,38,.10);
      }
      h1, h2, h3, h4 {
        color: var(--ink) !important;
      }
      .hero h1, .hero h2, .hero h3, .hero h4 {
        color: #609048 !important;
      }
      h1 {
        font-size: clamp(1.9rem, 2.4vw + 1rem, 2.6rem);
        line-height: 1.15;
        margin: 0.25rem 0;
      }
      h2 {
        font-size: clamp(1.2rem, 1.2vw + .9rem, 1.6rem);
        margin: 1.2rem 0 0.5rem;
      }
      h3 {
        font-size: 1.05rem;
        margin: 0.8rem 0 0.35rem;
      }
      .nav-link {
        color: var(--rw2) !important;
        font-size: 1.08rem;
      }
      .nav-link.active {
        color: var(--rw1) !important;
        font-weight: 600;
      }
      .btn-primary {
        background-color: #609048 !important;
        border-color: #609048 !important;
        color: white !important;
      }
      .btn-primary:hover {
        background-color: #90A860 !important;
        border-color: #90A860 !important;
        color: white !important;
      }
      .btn-danger {
        background-color: #786048 !important;
        border-color: #786048 !important;
        color: white !important;
      }
      .btn-danger:hover {
        background-color: #604830 !important;
        border-color: #604830 !important;
        color: white !important;
      }
      code {
        background-color: var(--rw4) !important;
        color: var(--rw2) !important;
        padding: 2px 6px !important;
        border-radius: 3px !important;
      }
      .card {
        border-left: 3px solid var(--rw1) !important;
      }
      details summary {
        background: var(--rw4) !important;
        border-left: 3px solid var(--rw1) !important;
      }
      details summary:hover {
        background: #e8e4df !important;
      }
      hr {
        border-color: var(--rw3) !important;
        opacity: 0.3;
      }
      .hero {
        background: transparent;
        color: var(--ink);
        padding: 44px 0 24px;
        margin: -20px -15px 24px -15px;
      }
      .badge {
        display: inline-block;
        padding: 0.18rem 0.55rem;
        border-radius: 999px;
        background: rgba(255,255,255,.16);
        border: 1px solid rgba(255,255,255,.28);
        font-size: 0.8rem;
        letter-spacing: 0.02em;
        margin-bottom: 0.5rem;
      }
      .divider {
        height: 2px;
        background: linear-gradient(90deg, transparent, var(--rw3), transparent);
        margin: 16px 0;
        border: none;
      }
      .about-card {
        background: #ffffff;
        border: 1px solid #e7e6e2;
        border-radius: var(--radius);
        padding: 14px 16px;
        box-shadow: var(--shadow);
        margin-bottom: 14px;
      }
      .tag {
        display: inline-block;
        background: var(--rw1);
        color: #fff;
        padding: 0.22rem 0.55rem;
        border-radius: 10px;
        font-size: 0.8rem;
        margin: 0.05rem 0.4rem 0.35rem 0;
      }
      .tag-uscr {
        background: #609048 !important;
      }
      .tag-rem {
        background: #786048 !important;
      }
      .tag-tte {
        background: #604830 !important;
      }
      .pcbox {
        background: #fff;
        border: 1px solid #e5e3df;
        border-radius: 12px;
        padding: 10px 12px;
        box-shadow: 0 8px 20px rgba(33,48,38,.06);
        margin-top: 0.5rem;
      }
      .pcbox .pcsec + .pcsec {
        margin-top: 0.5rem;
        padding-top: 0.5rem;
        border-top: 1px dashed #ddd6cd;
      }
      .pcbox h4 {
        margin: 0.1rem 0 0.25rem;
        font-size: 0.95rem;
        color: #3C4B3E;
      }
      .pcbox ul {
        margin: 0.1rem 0 0.1rem 1rem;
      }
      .kicker {
        font-weight: 600;
        color: var(--rw5);
        letter-spacing: 0.02em;
      }
      .ref {
        display: block;
        font-size: 0.86rem;
        color: var(--muted);
        margin-top: 0.25rem;
      }
      .lead {
        font-size: 1.05rem;
      }
      .small {
        font-size: 0.92rem;
        color: var(--muted);
      }
      ol.simple {
        padding-left: 1rem;
        margin: 0.6rem 0;
      }
      ol.simple > li {
        margin: 0.35rem 0;
        padding-left: 0.3rem;
      }
      ol.simple strong {
        color: var(--rw5);
      }
      a {
        color: var(--rw1);
        text-decoration: none;
        border-bottom: 1px dotted rgba(0,0,0,.18);
      }
      .nav-tabs {
        justify-content: center !important;
      }
      .banner-wrap {
        text-align: center;
        padding: 6px 0 12px;
      }
      .banner-image {
        width: min(100%, 1480px);
        height: auto;
        display: block;
        margin: 0 auto;
      }
    "))
  ),
  tags$div(
    class = "banner-wrap",
    tags$img(
      src = "deer_app_banner.png",
      alt = "DEER App banner",
      class = "banner-image"
    )
  ),
  
  navset_tab(
    id = "deer_tabs",
        
        # ------------------------- ABOUT ------------------------------
        nav_panel(
          "About",
          tags$div(
            class = "hero",
            tags$div(
              style = "max-width: 980px; margin: 0 auto; padding: 0 16px;",
              tags$h1(style = "color: var(--ink); margin-top: 0.5rem; text-align: center;",
                "DEER App"
              ),
              tags$h2(style = "color: var(--ink); margin-top: 0.5rem; font-size: 1.2rem; font-weight: 400; text-align: center;",
                "Density Estimation from Encounter Rates"
              ),
              tags$p(class = "small", style = "opacity: 0.9; color: var(--muted); text-align: center;",
                "USCR · REM · TTE — three unmarked camera methods for estimating animal density, reported as animals/km²."
              )
            )
          ),
          tags$div(
            style = "max-width: 980px; margin: 0 auto; padding: 0 16px 42px;",
            
            # Welcome section
            tags$section(
              tags$h1("Welcome to the DEER app!"),
              tags$h2("What it does"),
              tags$p(class = "lead",
                "The DEER app helps estimate ", tags$strong("animal density at sites"), " from ",
                tags$strong("unmarked camera detections"), ". It fits ", tags$strong("three complementary models"),
                ", reports ", tags$strong("credible intervals"), ", and lets users ",
                tags$strong("compare, diagnose, and average results"), ". The current upload workflow expects data on ",
                tags$strong("camera deployment"), " (e.g., where and when cameras were set and began recording) and ",
                tags$strong("images"), " (e.g., date and time animals were detected at each camera). The broader goal is a ",
                tags$strong("consistent, reproducible workflow"), " for estimating ", tags$strong("animal density"),
                " from camera-trap data without requiring individually identifiable animals."
              )
            ),
            
            tags$div(class = "divider"),
            
            # How to use
            tags$section(
              tags$h2("How to use the DEER App"),
              tags$div(
                class = "about-card",
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 1: Upload field data"),
                tags$p(
                  "The main analysis workflow uses uploaded field data. The app fits ",
                  tags$strong("Unmarked Spatial Capture-Recapture (USCR)"),
                  ", ",
                  tags$strong("Random Encounter Model (REM)"),
                  ", and ",
                  tags$strong("Time-to-Event (TTE)"),
                  " models from camera-trap data."
                ),
                tags$div(
                  style = "margin: 1rem 0;",
                  tags$p("You will need two data files:"),
                  tags$ul(
                    tags$li(tags$strong("Deployment CSV"), " — camera deployment information such as locations, dates, and detection distances."),
                    tags$li(tags$strong("Images CSV"), " — detection records such as timestamps, species, and ", tags$strong("Cluster ID"), " values. In this workflow, ", tags$strong("Cluster ID"), " means the unique identifier for an independent encounter event.")
                  ),
                  tags$p(
                    "The current upload pipeline expects the ", tags$strong("exact column names"), " described in the Add your data tab. Default camera spacing and camera counts are based on the ",
                    tags$a(
                      href = "https://irma.nps.gov/DataStore/Reference/Profile/2307775",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "northeastern National Park Service deer monitoring program"
                    ),
                    ", but the app analyzes the data you provide rather than requiring one exact array design or minimum camera count to function."
                  ),
                  tags$p(
                    "Click the ", tags$strong("'Add your data'"), " tab for column requirements and examples. The app will automatically:"
                  ),
                  tags$ul(
                    tags$li("Check required columns and data types."),
                    tags$li("Optionally trim each camera deployment length to meet study design criteria.")
                  )
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 2: Adjust settings (optional)"),
                tags$p(
                  "The default settings work for many use cases. To change priors or detection geometry, open the ",
                  tags$strong("'Model settings'"), " tab, choose ", tags$strong("Advanced"), ", and edit the fields."
                ),
                tags$ul(
                  tags$li("Modify priors for movement speed, viewshed or detection parameters, and camera heterogeneity."),
                  tags$li("Change the fallback camera detection angle, or provide camera-specific `Camera Detection Angle` values in the deployment file.")
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 3: Run the models"),
                tags$p(
                  tags$strong("Uploaded field data:"), " use the ", tags$strong("USCR"), ", ", tags$strong("REM"), ", and ", tags$strong("TTE"),
                  " tabs. Each tab has its own run button and troubleshooting panel."
                ),
                tags$p(
                  tags$strong("Simulated data:"), " run ", tags$strong("USCR"), ", ", tags$strong("REM"), ", and ", tags$strong("TTE"), " from their own tabs after generating the shared simulated dataset."
                ),
                tags$ul(
                  tags$li("Click the ", tags$strong("'Run'"), " button for your data type."),
                  tags$li("Watch stage labels, run-status notes, and troubleshooting text where available."),
                  tags$li("Use the red ", tags$strong("'Stop'"), " button if necessary and when a model supports stopping."),
                  tags$li("Results appear below the buttons once each model completes.")
                ),
                tags$p(
                  tags$em("Note:"), " USCR can take 30 minutes or longer for high-density sites, larger arrays, or larger ", tags$code("M"),
                  " settings, so plan to be patient and, if needed, leave the computer running. REM and TTE are usually much faster."
                ),
                tags$p(
                  "For longer runs, users may prefer to download the app and run it locally. Download and setup instructions are available in the ",
                  tags$a(
                    href = "https://github.com/kcring/DEER_app",
                    target = "_blank",
                    style = "color: var(--rw1); text-decoration: underline;",
                    "GitHub repository"
                  ),
                  "."
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 0.5rem;", "Step 4: Compare results"),
                tags$p("Use the ", tags$strong("'Compare & combine'"), " tab to:"),
                tags$ul(
                  tags$li("For ", tags$strong("uploaded field data"), ": compare whichever model fits have finished so far, and compute WAIC-weighted summaries when WAIC is available for all completed fits."),
                  tags$li("For ", tags$strong("simulated data"), ": compare whichever shared-simulation model fits have finished so far, as in the uploaded-data workflow."),
                  tags$li("Compare model performance using WAIC values for uploaded field data."),
                  tags$li("Download posterior summaries with parameter names, means, and 95% credible intervals.")
                ),
                tags$h2(style = "font-size: 1.5rem; font-weight: 600; margin-top: 1rem;", "Step 5: Optional design tools"),
                tags$p(
                  "The ", tags$strong("'Simulate data'"), " and ", tags$strong("'Camera array design'"), " tabs can be used as design tools."
                ),
                tags$ul(
                  tags$li("The shared simulator can generate practice datasets for USCR, REM, and TTE."),
                  tags$li("The camera-array tool can help lay out cameras for a new site from uploaded spatial layers.")
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # The three models
            tags$section(
              tags$h2("Quick model descriptions"),
              tags$div(
                style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 14px;",
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-uscr", "USCR"),
                  tags$p(
                    tags$strong("Unmarked Spatial Capture–Recapture (USCR)"), " — ", tags$em("Chandler & Royle, 2013"), ".",
                    tags$br(),
                    tags$strong("Spatial pattern → density."), " Where detections fall across the array informs density estimates inside the defined state-space."
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1214/12-AOAS610",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Chandler, R.B. & Royle, J.A. (2013)"
                    )
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-rem", "REM"),
                  tags$p(
                    tags$strong("Random Encounter Model (REM)"), " — ", tags$em("Rowcliffe et al., 2008"), ".",
                    tags$br(),
                    tags$strong("Events per unit time → density."), " Encounter rates are converted to density using animal speed and the camera viewshed."
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1111/j.1365-2664.2008.01473.x",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Rowcliffe, J.M. et al. (2008)"
                    )
                  )
                ),
                tags$div(
                  class = "about-card",
                  tags$span(class = "tag tag-tte", "TTE"),
                  tags$p(
                    tags$strong("Time‑to‑Event (TTE)"), " — ", tags$em("Moeller et al., 2018"), ".",
                    tags$br(),
                    tags$strong("Waiting time → density."), " Shorter effective time between encounters implies higher density when speed and viewshed are held fixed."
                  ),
                  tags$span(
                    class = "ref",
                    "Reference: ",
                    tags$a(
                      href = "https://doi.org/10.1002/ecs2.2331",
                      target = "_blank",
                      style = "color: var(--rw1); text-decoration: underline;",
                      "Moeller, A.K. et al. (2018, Ecosphere)"
                    )
                  )
                )
              )
            ),
            
            tags$div(class = "divider"),
            
            # What you get
            tags$section(
              tags$h2("What you get"),
              tags$div(
                class = "about-card",
              tags$p(
                tags$span(class = "kicker", "Current default output:"),
                " ",
                tags$strong("Animals per square kilometer (km²)"),
                " with credible intervals, plus per-model estimates, model averaging, and downloadable posterior summaries for transparency."
              )
              )
            ),
            
            tags$div(class = "divider"),
            
            # Acknowledgments
            tags$section(
              tags$h2("Acknowledgments"),
              tags$div(
                class = "about-card",
                tags$h3(style = "font-size: 1.2rem; font-weight: 600; margin-top: 0.5rem;", "Model Development"),
                tags$p(
                  "The underlying models and code were created by ",
                  tags$strong("Dr. Amanda Van Buskirk"),
                  " with guidance from ",
                  tags$a(
                    href = "https://www.davis.wvu.edu/faculty-staff/directory/christopher-rota",
                    target = "_blank",
                    style = "color: var(--rw1); text-decoration: underline;",
                    tags$strong("Dr. Christopher Rota")
                  ),
                  " at ",
                  tags$strong("West Virginia University, Davis College of Agriculture and Natural Resources"),
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  "Scientific collaboration and feedback on model integration and quality control also came from ",
                  tags$strong("Dr. Laura C. Gigliotti"),
                  ", ",
                  tags$strong("U.S. Geological Survey, West Virginia Cooperative Fish and Wildlife Research Unit, West Virginia University"),
                  "."
                )
              ),
              tags$div(
                class = "about-card",
                tags$h3(style = "font-size: 1.2rem; font-weight: 600; margin-top: 0.5rem;", "Shiny App Development"),
                tags$p(
                  "This Shiny application was developed by ",
                  tags$strong("PhD Candidate Kacie Ring"),
                  " as part of the ",
                  tags$a(href = "https://esa.org/programs/scip/", target = "_blank",
                         style = "color: var(--rw1); text-decoration: underline;",
                         tags$strong("SciComm in the Parks fellowship")),
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  "The ", 
                  tags$a(href = "https://esa.org/programs/scip/", target = "_blank", 
                         style = "color: var(--rw1); text-decoration: underline;",
                         tags$strong("SciComm in the Parks fellowship")), 
                  ", a collaborative effort between the ", 
                  tags$strong("Ecological Society of America (ESA)"), 
                  " and the ", 
                  tags$strong("National Park Service (NPS)"), 
                  "."
                ),
                tags$p(
                  style = "margin-top: 0.75rem;",
                  tags$strong("Fellowship Support:"), tags$br(),
                  "• ", tags$strong("Dr. Brian Mitchell"), " (NPS) — Fellowship Liaison", tags$br(),
                  "• ", tags$strong("Jasjeet Dhanota"), " (ESA) — Mentor", tags$br(),
                  "• ", tags$strong("Mary Joy Mulumba"), " (ESA) — Program Assistant"
                )
              )
            ),
            
            tags$div(class = "divider")
          )
        ),
        
        # ---------------------- SIMULATE DATA -------------------------
        nav_panel(
          "Simulate data",
          markdown(paste(
            "Generate a shared simulated dataset under **SECR** (*spatially explicit capture–recapture*) assumptions and use it to run **USCR, REM, and TTE** from their model tabs for comparison-ready simulated analyses.",
            "",
            "The underlying generator is spatial and most closely matches USCR, but the app reformats that same simulated dataset so REM and TTE can also be compared on it.",
            "",
            "Adjust the simulation parameters below, then analyze the simulated data using the model tabs.",
            sep = "\n"
          )),
          
          h3("Simulation parameters"),
          fluidRow(
            column(4,
              sliderInput("n_side", "Grid dimension (n x n cameras)",
                          min = 3, max = 8, value = 5, step = 1),
              sliderInput("spacing", "Camera spacing (m)",
                          min = 150, max = 500, value = 300, step = 10),
              sliderInput("days", "Survey length (days)",
                          min = 7, max = 60, value = 21, step = 1)
            ),
            column(4,
              sliderInput("Dtrue", "True density for simulation (animals/km²)",
                          min = 5, max = 80, value = 25, step = 1),
              sliderInput("home_range_km2", "Mean 95% home-range size (km²)",
                          min = 0.10, max = 3.00, value = 0.89, step = 0.05),
              sliderInput("lambda0", "Detection rate λ₀ (expected detections/day at camera center)",
                          min = 0.05, max = 0.6, value = 0.20, step = 0.01)
            ),
            column(4,
              numericInput("seed", "Random seed", value = 1, min = 1),
              numericInput("sim_buffer_m", "State-space buffer around camera array (m)",
                           value = 850, min = 50, step = 50),
              p(class = "small", style = "color: var(--muted);",
                "This buffer controls the area available to simulated animal activity centers. Smaller animals may need a smaller buffer; larger-ranging animals may need a larger one."),
              tags$details(
                class = "about-card",
                style = "margin-top: 0.75rem;",
                tags$summary(
                  style = "cursor: pointer; font-weight: 600;",
                  "REM/TTE model-input defaults for the shared simulator"
                ),
                tags$div(
                  style = "margin-top: 0.75rem;",
                  sliderInput("r_m_sim", "Detection radius for REM/TTE model inputs (m)",
                              min = 8, max = 25, value = 12, step = 1),
                  p(class = "small", style = "color: var(--muted); margin-bottom: 0;",
                    "Detection radius is not used to generate the spatial detections themselves. It is only carried into the shared simulated model-input table used by REM and TTE. The fallback detection angle for simulated REM/TTE runs is set under ",
                    tags$strong("Model settings"),
                    ".")
                )
              ),
              br(),
              actionButton("run_sim", "Simulate data", class = "btn-primary")
            )
          ),
          tags$details(
            class = "about-card",
            style = "margin-top: 1rem;",
            tags$summary(
              style = "cursor: pointer; font-weight: 600;",
              "Simulation parameter guide"
            ),
            tags$div(
              style = "margin-top: 0.85rem;",
              tags$p(
                tags$strong("Grid dimension"),
                " sets the number of cameras along one side of the toy grid. A 5 x 5 grid gives 25 cameras."
              ),
              tags$p(
                tags$strong("Camera spacing"),
                " is the distance between neighboring cameras. Cameras should be close enough together for an individual to be detected at multiple adjacent cameras."
              ),
              tags$p(
                tags$strong("Survey length"),
                " is the number of deployed days simulated for each camera."
              ),
              tags$p(
                tags$strong("True density"),
                " is the underlying animal density used to generate the toy data. This input and the app summaries are both reported in animals/km²."
              ),
              tags$p(
                tags$strong("Mean 95% home-range size"),
                " is converted in the background to the USCR spatial scale parameter σ. Larger home ranges spread detections across more cameras and correspond to broader space use."
              ),
              tags$p(
                tags$strong("Baseline detection rate λ₀"),
                " is the expected number of detections per deployed day if an animal's activity center were at the camera. Larger values produce more detections everywhere in the array."
              ),
              tags$p(
                tags$strong("State-space buffer"),
                " is the distance added around the camera array when placing simulated activity centers. It should be large enough to contain realistic space use for the species being simulated."
              ),
              tags$p(
                tags$strong("Random seed"),
                " reproduces the same simulated dataset when you rerun the grid with identical settings."
              ),
              tags$p(
                tags$strong("Detection radius r"),
                " is not used to generate the spatial detections themselves. It is carried into the shared model-input table used by REM and TTE."
              )
            )
          ),
          
          hr(),
          h3("Simulated visualizations"),
          fluidRow(
            column(6,
              h4("Simulated field grid"),
              plotOutput("grid_plot", height = "420px")
            ),
            column(6,
              h4("Simulated animal distribution"),
              plotOutput("sim_deer_distribution_plot", height = "420px")
            )
          ),
          h3("Camera data table"),
          DTOutput("camera_table"),
        ),

        # ------------------- CAMERA ARRAY DESIGN --------------------
        nav_panel(
          "Camera array design",
          markdown(paste(
            "Design a camera layout for a site from uploaded spatial layers.",
            "Upload a boundary and any optional exclusion layers, choose spacing and buffer settings, and the app will suggest a camera array you can preview and export.",
            "",
            "Use this as a planning tool: the app can lay out the full suggested camera network, or a smaller final set if you only have a limited number of cameras available to deploy.",
            sep = "\n"
          )),
          tags$div(
            id = "camera_array_form",
            fluidRow(
            column(
              4,
              h3("Step 1: Upload site layers"),
              p(
                "Boundary is required. Optional layers can keep cameras away from roads, trails, buildings, parking lots, or other excluded areas."
              ),
              tags$p(
                class = "small",
                style = "color: var(--muted);",
                "If your site has multiple units, you can upload multiple boundary or optional-layer files in the same input and the app will merge them."
              ),
              tags$p(
                class = "small",
                style = "color: var(--muted);",
                "Recommended option: upload one zipped layer per input. You can also upload shapefile parts together (.shp, .dbf, .shx, .prj, and optional .cpg), or use KML, KMZ, or GeoPackage files."
              ),
              fileInput(
                "camera_boundary_files",
                "Boundary file(s)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_roads_files",
                "Roads layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_trails_files",
                "Trails layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_buildings_files",
                "Buildings layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_parking_files",
                "Parking lots layer (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              ),
              fileInput(
                "camera_exclusion_files",
                "Additional exclusion layer(s) (optional)",
                multiple = TRUE,
                accept = c(".zip", ".shp", ".dbf", ".shx", ".prj", ".cpg", ".kml", ".kmz", ".gpkg")
              )
            ),
            column(
              4,
              h3("Step 2: Choose layout settings"),
              textInput("camera_prefix", "Camera name prefix", value = "CAM"),
              sliderInput(
                "camera_spacing_m",
                "Maximum target spacing between nearby cameras (m)",
                min = 100, max = 600, value = 400, step = 10
              ),
              sliderInput(
                "camera_buffer_m",
                "Buffer around roads / trails / buildings / parking (m)",
                min = 0, max = 100, value = 10, step = 1
              ),
              numericInput(
                "camera_budget",
                "Number of cameras to deploy (0 = keep all suggested cameras)",
                value = 0,
                min = 0,
                step = 1
              ),
              numericInput(
                "camera_alternates",
                "Alternate cameras to keep when budget is smaller than the suggested layout",
                value = 2,
                min = 0,
                step = 1
              ),
              actionButton("run_camera_array", "Generate camera array", class = "btn-primary"),
              actionButton("reset_camera_array", "Clear uploads and reset", style = "margin-left: 10px;")
            ),
            column(
              4,
              h3("What this produces"),
              tags$ul(
                tags$li("A projected design area with optional exclusion buffers removed."),
                tags$li("Suggested camera points spaced across the allowable area."),
                tags$li("A final camera set, plus alternates if you enter a smaller budget."),
                tags$li("Downloadable camera coordinates in latitude/longitude and UTM.")
              ),
              tags$p(
                class = "small",
                style = "color: var(--muted);",
                "If the suggested layout uses more cameras than you have available, enter your available camera count and the app will keep a well-spread final subset plus optional alternates. If that smaller set cannot preserve the original spacing goal, the app will flag it in the log and summary."
              )
            )
          )),
          hr(),
          h3("Camera-array log"),
          verbatimTextOutput("camera_array_log"),
          fluidRow(
            column(
              5,
              h3("Layout summary"),
              verbatimTextOutput("camera_array_summary"),
              downloadButton("download_camera_array_csv", "Download camera CSV"),
              tags$span(" "),
              downloadButton("download_camera_array_gpkg", "Download camera GeoPackage"),
              tags$span(" "),
              downloadButton("download_camera_area_gpkg", "Download design-area GeoPackage")
            ),
            column(
              7,
              h3("Camera-array map"),
              leafletOutput("camera_array_map", height = "520px")
            )
          ),
          h3("Camera coordinates"),
          DTOutput("camera_array_table")
        ),
        
        # ---------------------- ADD YOUR DATA -------------------------
        nav_panel(
          "Add your data",
          markdown(paste(
            "Upload a **deployment CSV** (where and when cameras were set and recording) and an **images CSV** (timestamps, species, counts, and Cluster IDs).",
            "The current upload checker expects the exact column names listed below.",
            "",
            sep = "\n"
          )),
          h3("The app will:"),
          markdown(paste(
            "1. Check required columns and data types;",
            "2. Flag image timestamps that fall outside the deployment window;",
            "3. Optionally trim each camera deployment length to meet study design criteria.",
            "",
            "Download processed CSVs if needed, then configure MCMC and priors in the **Model settings** tab and run models in the USCR/REM/TTE tabs.",
            "",
            "---",
            "",
            "## Required data columns",
            "",
            "### Deployment file columns (**bold = required**, others optional for quality control):",
            "",
            "- **`Site Name`** — Camera site identifier",
            "- `Site` — Optional broader site identifier",
            "- `Camera ID` — Optional camera identifier",
            "- `SD Card ID` — Optional SD card identifier",
            "- **`Start Date`** — Deployment start date (`MM/DD/YYYY`)",
            "- **`Start Time`** — Deployment start time (`HH:MM:SS`)",
            "- **`End Date`** — Deployment end date (`MM/DD/YYYY`)",
            "- **`End Time`** — Deployment end time (`HH:MM:SS`)",
            "- **`Latitude`** — Camera latitude (decimal degrees)",
            "- **`Longitude`** — Camera longitude (decimal degrees)",
            "- `Camera Model` — Optional camera make/model field for recordkeeping (for example `Browning Strike Force Pro`)",
            "- `Camera Detection Angle` — Optional full detection angle in degrees for that camera. If left blank, the app uses the fallback angle from Model settings.",
            "- `Camera Height` — Optional camera height (meters)",
            "- `Camera Orientation` — Optional cardinal direction (for example `N`, `NE`, `E`, `SE`, `S`, `SW`, `W`, `NW`) or 0-359 degrees",
            "- **`Camera Functioning`** — Camera status; `Yes/No`, `TRUE/FALSE`, `T/F`, and `1/0` values are accepted",
            "- **`Camera Malfunction Date`** — Keep this value blank for cameras that did not malfunction. Enter the failure date (`MM/DD/YYYY`) if the camera failed.",
            "- **`Detection Distance`** — Effective detection radius in meters",
            "",
            "### **Images file** required columns:",
            "",
            "- **`Site Name`** — Camera site identifier. These identifiers must match identifiers in the deployment file.",
            "- **`Timestamp`** — Image timestamp. Preferred format: `MM/DD/YYYY HH:MM:SS`.",
            "- **`Species`** — Species identifier. Pipe-delimited values such as `deer|squirrel` are allowed for multi-species rows, and the order must match `Sighting Count`.",
            "- **`Cluster ID`** — Unique identifier for independent detection events",
            "- **`Sighting Count`** — Number of individuals in the image; pipe-delimited values such as `1|3` are allowed for multi-species rows and should follow the same order as `Species`.",
            "- `Image URL` — Optional image reference/link column. If included, store it as plain text in a single column, one value per row.",
            "",
            "**Notes:** Deployment QC expects date-only fields in `MM/DD/YYYY`. Image timestamps are parsed more flexibly if needed, including 2-digit years (for example `2/3/25`).",
            "The app expects the exact column names listed above. Files from other tagging workflows can usually be renamed to match.",
            "Cross-year winter surveys (for example `12/2025` to `01/2026`) are supported; the app uses the actual deployment dates/times and image timestamps, so no separate `Survey Year` field is required.",
            sep = "\n"
          )),
          
          h3("Step 1: Deployment file"),
          fileInput("deployment_csv", "Upload deployment CSV", accept = ".csv"),
          h4("Deployment check log"),
          verbatimTextOutput("deployment_check_log"),
          h4("Preview of cleaned deployment data"),
          DTOutput("deployment_preview"),
          downloadButton("download_deployment_checked", "Download cleaned deployment CSV"),
          
          hr(),
          
          h3("Step 2: Images file"),
          fileInput("images_csv", "Upload images CSV", accept = ".csv"),
          checkboxInput(
            "apply_56day_trim",
            "After validation, trim each camera deployment length to a user-set number of days.",
            value = TRUE
          ),
          conditionalPanel(
            "input.apply_56day_trim",
            numericInput(
              "trim_days",
              "Deployment length to keep (days)",
              value = 56,
              min = 1,
              step = 1
            )
          ),
          h4("Images check log"),
          verbatimTextOutput("images_check_log"),
          h4("Preview of processed images data"),
          DTOutput("images_preview"),
          downloadButton("download_images_checked", "Download processed images CSV"),
          
          hr(),
          
          h3("Model settings"),
          p(
            "Configure priors and fallback detection-angle settings (",
            tags$strong("θ"), ") in the ",
            tags$strong("Model settings"), " tab."
          )
        ),
        
        nav_panel(
          "Model settings",
          markdown(paste(
            "### Run defaults",
            "",
            "**Current app defaults**:",
            "",
            "- **Number of chains**: 2 by default; the app enforces a minimum of 2 chains for model runs",
            "- **REM/TTE**: adaptive runs start from 6000 iterations, 1000 burn-in, and thin = 5",
            "- **USCR**: adaptive runs start from 6000 iterations, 1000 burn-in, thin = 5, and M = 300",
            "- **Convergence criteria**: R̂ < 1.1 for all monitored parameters",
            "",
            "The app handles chain counts and reruns internally, so users do not need to manually edit iterations, thinning, or data-augmentation size in the interface.",
            "If convergence is still poor after a run, review the run-status panel and consider whether the priors need to better match your site or species.",
            "",
            "### Prior distributions",
            "",
            "The priors below are grouped by model type.",
            "The current defaults are calibrated for white-tailed deer using home-range size and daily movement information from the literature.",
            "These defaults may not be well calibrated for other species or field conditions.",
            "Model output can be sensitive to prior choice, so users working on other species should adjust the priors to better match their system.",
            "",
            "### Fallback camera geometry",
            "",
            "Fallback detection angle is 55° (based on Browning-style camera specifications). This value is used for simulated runs and for uploaded cameras that do not include a `Camera Detection Angle` value.",
            "",
            sep = "\n"
          )),
          
          h4("Fallback camera geometry"),
          fluidRow(
            column(12,
              sliderInput("theta", "Fallback detection angle θ (degrees)",
                          min = 20, max = 80, value = 55, step = 1)
            )
          ),
          
          hr(),
          
          radioButtons(
            "mode", "Settings mode",
            choices  = c("Simple (recommended)" = "Default", "Advanced" = "Advanced"),
            selected = "Default",
            inline   = TRUE
          ),
          
          conditionalPanel(
            "input.mode == 'Default'",
            tags$div(
              class = "about-card",
              tags$p(
                tags$strong("Using the app defaults."),
                " These settings are meant to be a practical starting point for most users."
              ),
              tags$ul(
                tags$li("2 chains are used for all model runs."),
                tags$li("REM and TTE start from 6000 iterations, 1000 burn-in, and thin = 5."),
                tags$li("USCR starts from 6000 iterations, 1000 burn-in, thin = 5, and M = 300."),
                tags$li("The app checks convergence internally and can rerun models when needed."),
                tags$li("Most users only need to change the fallback camera angle above or switch to Advanced for species- or site-specific priors.")
              )
            )
          ),
          
          conditionalPanel(
            "input.mode == 'Advanced'",
            tagList(
              tags$p(
                class = "small",
                style = "color: var(--muted);",
                "Most users can leave these advanced controls alone. The app manages chain counts and reruns internally, so the advanced controls focus on priors rather than MCMC tuning."
              ),
              tags$details(
                tags$summary(tags$strong("Priors: REM & TTE")),
                markdown(paste(
                  "- **D** ~ Uniform(0, D_max): Uniform prior for density. `D_max` should be well above the maximum value of the posterior draws from `D`.",
                  "- **log(v)** ~ Normal(mean, SD): Informative prior on daily movement speed (log scale). The current default corresponds to an average daily movement rate of `3.09 km/day` with an approximate 95% interval of `1.60` to `6.00 km/day`.",
                  "- **sd_eps** ~ Gamma(shape, rate): Prior on the standard deviation of the camera-level lognormal random effect in the Poisson-lognormal count model. The default `Gamma(1, 1)` prior has mean 1 and variance 1. These shared `sd_eps` controls are used by REM, TTE, and USCR.",
                  "",
                  "These defaults are meant for white-tailed deer. If you are working with another species, adjust the priors to match that species and study design.",
                  "",
                  sep = "\n"
                )),
                fluidRow(
                  column(6,
                    numericInput("D_max", "Max density prior (D max, animals/km²)",
                                 value = 200, min = 10, step = 10),
                    numericInput("log_v_mean", "log(v) mean (km/day)",
                                 value = 1.130, step = 0.1),
                    numericInput("log_v_sd", "log(v) SD",
                                 value = 0.3372, min = 0.01, step = 0.05)
                  ),
                  column(6,
                    numericInput("sd_eps_shape", "sd_eps gamma shape",
                                 value = 1, min = 0.1, step = 0.1),
                    numericInput("sd_eps_rate", "sd_eps gamma rate",
                                 value = 1, min = 0.1, step = 0.1)
                  )
                )
              ),
              
              tags$details(
                tags$summary(tags$strong("Priors: USCR")),
                markdown(paste(
                  "- **M = 300**: Starting upper bound on total abundance, not density. Larger values increase run time, and the app checks whether this starting value is large enough during fitting.",
                  "- **log(σ)** ~ Normal(mean, SD): Informative prior on animal space use (log scale). The current default corresponds to an average 95% circular home-range size of `0.89 km²` with an approximate 95% interval of `0.48` to `1.62 km²`.",
                  "- **log(λ₀)** ~ Normal(mean, SD): Prior on expected detections per day when an animal's activity center is at the camera. The default mean of `0` corresponds to about one expected detection per deployed day at the detector because `exp(0) = 1`.",
                  "- **State-space buffer**: Distance added around the camera array when defining the USCR state-space. Smaller species may need a smaller buffer; larger-ranging species may need a larger one.",
                  "- **sd_eps** ~ Gamma(shape, rate): Prior on the standard deviation of the camera-level lognormal random effect. USCR uses the same shared `sd_eps` gamma controls listed in the REM & TTE section above.",
                  "",
                  "These defaults are meant for white-tailed deer. If you are working with another species, adjust the priors to better match that species and study design.",
                  "",
                  sep = "\n"
                )),
                fluidRow(
                  column(6,
                    numericInput("log_sigma_mean", "log(σ) mean",
                                 value = -1.5269, step = 0.1),
                    numericInput("log_sigma_sd", "log(σ) SD",
                                 value = 0.1535, min = 0.01, step = 0.05),
                    numericInput("log_lam0_mean", "log(λ₀) mean",
                                 value = 0, step = 0.1)
                  ),
                  column(6,
                    numericInput("uscr_buffer_m", "USCR state-space buffer (m)",
                                 value = 1200, min = 50, step = 50),
                    numericInput("log_lam0_sd", "log(λ₀) SD",
                                 value = 1, min = 0.1, step = 0.1)
                  )
                )
              )
            )
          ),
          
          # Hidden run defaults used in both modes
          tags$div(
            style = "display:none;",
            tags$div(
              numericInput("n_chains", NULL, value = 2),
              numericInput("iter_rem_tte", NULL, value = 6000),
              numericInput("burnin_rem_tte", NULL, value = 1000),
              numericInput("thin_rem_tte", NULL, value = 5),
              numericInput("iter_uscr", NULL, value = 6000),
              numericInput("burnin_uscr", NULL, value = 1000),
              numericInput("thin_uscr", NULL, value = 5),
              numericInput("M_uscr", NULL, value = 300),
              numericInput("D_max", NULL, value = 200),
              numericInput("log_v_mean", NULL, value = 1.130),
              numericInput("log_v_sd", NULL, value = 0.3372),
              numericInput("log_sigma_mean", NULL, value = -1.5269),
              numericInput("log_sigma_sd", NULL, value = 0.1535),
              numericInput("uscr_buffer_m", NULL, value = 1200),
              numericInput("log_lam0_mean", NULL, value = 0),
              numericInput("log_lam0_sd", NULL, value = 1),
              numericInput("sd_eps_shape", NULL, value = 1),
              numericInput("sd_eps_rate", NULL, value = 1)
            )
          )
        ),
        
        # ---------------------- DATA SUMMARY --------------------------
        nav_panel(
          "Data summary",
          p(
            class = "small",
            style = "color: var(--muted);",
            "This tab summarizes uploaded deployment and image data. Simulated datasets are summarized in their own model tabs and the Compare & combine tab."
          ),
          h4("Site deployment summary"),
          DTOutput("deploy_summary_table"),
          leafletOutput("camera_map", height = "300px"),
          hr(),
          h4("Uploaded detections by species"),
          DTOutput("image_summary_table"),
          plotOutput("species_bar_plot", height = "300px"),
          hr(),
          h4("Selected species detections by site"),
          selectInput(
            "summary_species",
            "Species to summarize",
            choices = character(0)
          ),
          DTOutput("deer_summary_table"),
          plotOutput("deer_bubble_plot", height = "300px"),
          plotOutput("deer_daily_plot", height = "300px")
        ),
        
        # ---------------------- USCR MODEL TAB ------------------------
        nav_panel(
          "USCR model",
          tags$div(
            id = "uscr-content",
            HTML('
            <h2>Model 1 — USCR (Unmarked Spatial Capture–Recapture)</h2>
            <p><strong>The gist:</strong> Estimates <strong>animal density in the site/state‑space</strong> by learning from <strong>where and when</strong> animals are detected <strong>across a camera array</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera observation model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\mu_j), \\qquad \\log \\mu_j = \\log(\\mathrm{days}_j) + \\log(\\Lambda_j) + \\epsilon_j$$</p>
                <p>Expected encounter rate from augmented individuals:</p>
                <p>$$\\Lambda_j = \\sum_{i=1}^{M} z_i\\,\\lambda_0\\,\\exp\\!\\Big(-\\frac{d_{ij}^2}{2\\sigma^2}\\Big), \\qquad z_i \\sim \\mathrm{Bernoulli}(\\psi)$$</p>
                <p>Activity centers and camera effects:</p>
                <p>$$s_i \\sim \\mathrm{Uniform}(S), \\qquad \\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)$$</p>
                <p>Here \\(S\\) is the buffered state-space. Density is obtained from \\(N = \\sum_i z_i\\) divided by the study-area size. With real coordinates, area comes from the buffered camera state-space; on the toy simulator, area comes from the rectangular simulated state-space. The state-space buffer can be adjusted in Model settings. The app currently reports \\(D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — total independent animal detections at camera \\(j\\) across the deployment.</li>
                  <li>\\(\\mathrm{days}_j\\) — number of deployed days for camera \\(j\\).</li>
                  <li>\\(\\mu_j\\) — expected number of animal detections at camera \\(j\\) during the deployment.</li>
                  <li>\\(\\lambda_0\\) — expected detections per deployed day if an animal\'s activity center were at the camera.</li>
                  <li>\\(\\sigma\\) — spatial scale parameter (km), related to home-range size and how quickly detection falls with distance.</li>
                  <li>\\(s_i\\) — activity center of augmented individual \\(i\\), located within the buffered state-space \\(S\\).</li>
                  <li>\\(d_{ij}\\) — projected distance (km) from augmented individual \\(i\\)\'s activity center to camera \\(j\\).</li>
                  <li>\\(M\\) — upper limit of the prior distribution for total abundance in the augmented population.</li>
                  <li>\\(z_i\\) — indicator that augmented individual \\(i\\) is part of the real population.</li>
                </ul>
                <h3>Default priors used in the app (can be changed in Model settings)</h3>
                <ul>
                  <li>\\(M = 300\\), the starting upper bound on total abundance used for data augmentation. If the posterior presses against \\(M\\), the app increases it and reruns the model.</li>
                  <li>\\(\\log \\sigma \\sim \\mathcal{N}(-1.5269,\\,0.1535)\\), which corresponds to an average 95% circular home-range size of about 0.89 km².</li>
                  <li>\\(\\log \\lambda_0 \\sim \\mathcal{N}(0,\\,1)\\), centered near one expected detection per deployed day when the activity center is at the camera.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathrm{Gamma}(1,1)\\), giving the camera-level random-effect standard deviation a prior mean of 1 and variance of 1.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('uscr-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Uploaded field data"),
          div(
            actionButton("run_uscr_nps", "Run USCR on uploaded data", class = "btn-primary"),
            actionButton("pause_resume_uscr_nps", "Pause/Resume", class = "btn-warning", style = "margin-left: 10px;"),
            actionButton("stop_uscr_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          checkboxInput(
            "uscr_fast_mode_nps",
            "Fast test mode for uploaded USCR",
            value = FALSE
          ),
          p(
            class = "small",
            style = "max-width: 52rem; color: var(--muted); margin-top: -0.2rem;",
            "Fast test mode is for quick checks only. It uses lighter settings and caps adaptive reruns so you can test whether the uploaded-data pipeline works without waiting for a full production run."
          ),
          p(
            class = "small",
            style = "max-width: 52rem; color: var(--muted); margin-top: -0.2rem;",
            "Pause waits for the current USCR round to finish, then saves a provisional estimate from the most recent completed round. The same button switches to Resume when a paused USCR fit is available. Provisional paused fits are useful for troubleshooting, but are excluded from Compare & combine by default."
          ),
          br(),
          verbatimTextOutput("uscr_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("uscr_nps_debug"),
          
          hr(),
          h4("Simulated data"),
          p(
            style = "max-width: 52rem;",
            "USCR runs on the",
            tags$strong("shared spatial simulator"),
            " generated in the Simulate data tab, so it can be compared and combined with simulated REM and TTE."
          ),
          div(
            actionButton("run_uscr_sim", "Run USCR on simulated data", class = "btn-primary"),
            actionButton("stop_uscr_sim", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("uscr_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("uscr_sim_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- REM MODEL TAB -------------------------
        nav_panel(
          "REM model",
          tags$div(
            id = "rem-content",
            HTML('
            <h2>Model 2 — REM (Random Encounter Model)</h2>
            <p><strong>The gist:</strong> Converts <strong>how often animals pass a camera</strong> into density, correcting for <strong>movement speed</strong> and <strong>camera view geometry</strong>.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera Poisson model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j)$$</p>
                <p>$$\\log \\lambda_j = \\log D + \\log(\\mathrm{days}_j) + \\log v + \\log r_j + \\log\\!\\Big(\\frac{2 + \\theta_{\\mathrm{rad},j}}{\\pi}\\Big) + \\epsilon_j$$</p>
                <p>$$\\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)$$</p>
                <p>The app currently reports density as \\(D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — number of independent animal detection events at camera \\(j\\).</li>
                  <li>\\(\\mathrm{days}_j\\) — deployed days for camera \\(j\\).</li>
                  <li>\\(D\\) — animal density in animals/km².</li>
                  <li>\\(v\\) — animal movement speed (km/day).</li>
                  <li>\\(r_j\\) — effective detection radius for camera \\(j\\) (km).</li>
                  <li>\\(\\theta_{\\mathrm{rad},j}\\) — full detection angle used for camera \\(j\\) in the REM formula (radians).</li>
                  <li>\\(\\epsilon_j\\) — camera-level random effect for overdispersed counts.</li>
                </ul>
                <h3>Default priors/inputs used in the app</h3>
                <ul>
                  <li>\\(D \\sim \\mathcal{U}(0, D_{\\max})\\), a uniform prior for density. \\(D_{\\max}\\) should be well above the maximum value of the posterior draws from \\(D\\).</li>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), which corresponds to an average daily movement rate of about 3.09 km/day.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathrm{Gamma}(1, 1)\\), the default prior on the standard deviation of the camera-level lognormal random effect.</li>
                  <li>Uploaded field data can provide a camera-specific `Camera Detection Angle`; otherwise the app uses the fallback angle from Model settings.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('rem-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Uploaded field data"),
          div(
            actionButton("run_rem_nps", "Run REM on uploaded data", class = "btn-primary"),
            actionButton("stop_rem_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("rem_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("rem_nps_debug"),
          
          hr(),
          h4("Simulated data"),
          p(
            style = "max-width: 52rem;",
            "REM runs on the",
            tags$strong("shared spatial simulator"),
            " generated in the Simulate data tab, so it can be compared and combined with simulated USCR and TTE."
          ),
          div(
            actionButton("run_rem_sim", "Run REM on simulated data", class = "btn-primary"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("rem_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("rem_sim_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- TTE MODEL TAB -------------------------
        nav_panel(
          "TTE model",
          tags$div(
            id = "tte-content",
            HTML('
            <h2>Model 3 — TTE (Time‑to‑Event)</h2>
            <p><strong>The gist:</strong> Uses <strong>animal detection events per camera</strong>, scaled by <strong>movement-based time units</strong> and the <strong>viewshed area</strong>. Shorter effective time between encounters implies higher density.</p>
            <details>
              <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
              <div style="margin: 1rem 0; padding-left: 1rem;">
                <p>Per-camera Poisson model:</p>
                <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j \\;=\\; \\log D \\;+\\; \\log(\\mathrm{days}_j)\\;+\\;\\log(U_j)\\;+\\;\\log(A_j)\\;+\\;\\epsilon_j$$</p>
                <p>Movement-based encounter multiplier, viewshed area, and camera effect:</p>
                <p>$$U_j \\;=\\; \\frac{v}{f(\\theta_j)\\,r_j}, \\qquad A_j \\;=\\; \\pi r_j^2 \\frac{\\theta_{\\mathrm{deg},j}}{360}, \\qquad \\epsilon_j \\sim \\mathcal{N}(0, sd_\\epsilon)$$</p>
                <p>with \\(f(\\theta_j) = 0.3324 + 0.005580\\,\\theta_j - 1.454\\times10^{-5}\\,\\theta_j^2\\), where \\(\\theta_j\\) is the full camera detection angle in degrees.</p>
                <p>The app currently reports density as \\(D_{\\mathrm{km}^2}\\).</p>
                <h3>Variables</h3>
                <ul>
                  <li>\\(y_j\\) — number of animal detection events for camera \\(j\\).</li>
                  <li>\\(\\mathrm{days}_j\\) — total deployed days for camera \\(j\\).</li>
                  <li>\\(U_j\\) — movement-based encounter-rate multiplier for camera \\(j\\).</li>
                  <li>\\(v\\) — movement speed (km/day); prior on \\(\\log v\\) as below.</li>
                  <li>\\(r_j\\) — effective detection radius for camera \\(j\\) (km).</li>
                  <li>\\(\\theta_{\\mathrm{deg},j}\\) — full detection angle used for camera \\(j\\) (degrees).</li>
                  <li>\\(A_j\\) — viewshed area for camera \\(j\\) (km²).</li>
                  <li>\\(\\epsilon_j\\) — camera-level random effect for overdispersed counts.</li>
                </ul>
                <h3>Default priors/inputs used in the app</h3>
                <ul>
                  <li>\\(D \\sim \\mathcal{U}(0, D_{\\max})\\), a uniform prior for density. \\(D_{\\max}\\) should be well above the maximum value of the posterior draws from \\(D\\).</li>
                  <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), which corresponds to an average daily movement rate of about 3.09 km/day.</li>
                  <li>\\(sd_\\epsilon \\sim \\mathrm{Gamma}(1, 1)\\), the default prior on the standard deviation of the camera-level lognormal random effect.</li>
                  <li>Uploaded field data can provide a camera-specific `Camera Detection Angle`; otherwise the app uses the fallback angle from Model settings.</li>
                </ul>
              </div>
            </details>
            '),
            tags$script(HTML("
              function processMathJax(element) {
                if (window.MathJax) {
                  if (window.MathJax.typesetPromise) {
                    MathJax.typesetPromise([element]).catch(function(err) {
                      console.log('MathJax error:', err);
                    });
                  } else if (window.MathJax.typeset) {
                    MathJax.typeset([element]);
                  }
                }
              }
              setTimeout(function() {
                var content = document.getElementById('tte-content');
                if (content) {
                  processMathJax(content);
                  var details = content.querySelector('details');
                  if (details) {
                    details.addEventListener('toggle', function(e) {
                      if (e.target.open) {
                        setTimeout(function() {
                          processMathJax(content);
                        }, 300);
                      }
                    });
                  }
                }
              }, 500);
            "))
          ),
          hr(),
          h4("Uploaded field data"),
          div(
            actionButton("run_tte_nps", "Run TTE on uploaded data", class = "btn-primary"),
            actionButton("stop_tte_nps", "Stop", class = "btn-danger", style = "margin-left: 10px;"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_nps_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("tte_nps_debug"),
          
          hr(),
          h4("Simulated data"),
          p(
            style = "max-width: 52rem;",
            "TTE runs on the",
            tags$strong("shared spatial simulator"),
            " generated in the Simulate data tab, so it can be compared and combined with simulated USCR and REM."
          ),
          div(
            actionButton("run_tte_sim", "Run TTE on simulated data", class = "btn-primary"),
            style = "margin-bottom: 10px;"
          ),
          br(),
          verbatimTextOutput("tte_sim_text"),
          h5("Run status & troubleshooting"),
          verbatimTextOutput("tte_sim_debug"),
          tags$div(
            style = "text-align: center; margin-top: 40px; padding-top: 20px; border-top: 1px solid #e0e0e0;",
            tags$h3(style = "margin: 0; color: var(--rw3);", "DEER App"),
            tags$p(style = "margin: 5px 0 0 0; color: var(--muted);", "Density Estimation from Encounter Rates")
          )
        ),
        
        # ---------------------- COMPARE & COMBINE --------------------
        nav_panel(
          "Compare & combine",
          markdown(paste(
            "**Uploaded field data:** the table updates as REM, TTE, and USCR finish. If only some models have completed, the table will still summarize the completed fits.",
            "",
            "1. Compute ΔWAIC and WAIC weights when multiple models are available and WAIC is available for all completed fits;",
            "2. Combine posterior draws of density (animals/km²) across the completed fits;",
            "3. Report model-specific densities, and when possible also report unweighted and WAIC-weighted summaries plus the probability that density exceeds the user-set threshold. If WAIC is missing for one or more completed fits, only the unweighted summary is shown.",
            "",
            "**Simulated data:** if you run USCR, REM, and TTE from the **shared spatial simulator**, the table below will compare and combine those completed simulated fits too.",
            "",
            "Run the models from their tabs first, then check the tables.",
            sep = "\n"
          )),
          numericInput(
            "combo_density_threshold",
            "Density threshold for probability summary (animals/km²)",
            value = 20,
            min = 0,
            step = 1
          ),
          
          h4("Uploaded field data – model comparison and combined results (animals/km²)"),
          DTOutput("nps_combo_table"),
          h5("Visual summaries"),
          plotOutput("nps_combo_interval_plot", height = "320px"),
          div(
            downloadButton("dl_nps_combo_interval_png", "Download uploaded-data interval plot (PNG)"),
            style = "margin-bottom: 12px;"
          ),
          plotOutput("nps_combo_density_plot", height = "320px"),
          div(
            downloadButton("dl_nps_combo_density_png", "Download uploaded-data posterior overlay (PNG)"),
            style = "margin-bottom: 12px;"
          ),
          p("Posterior summaries (all monitored parameters, mean, 95% CI) for each completed model:"),
          downloadButton("dl_nps_all_csv", "Download available uploaded-data posterior summaries (CSV)"),

          hr(),
          h4("Simulated data – shared spatial simulator comparison (animals/km²)"),
          DTOutput("sim_combo_table"),
          h5("Visual summaries"),
          plotOutput("sim_combo_interval_plot", height = "320px"),
          div(
            downloadButton("dl_sim_combo_interval_png", "Download simulated interval plot (PNG)"),
            style = "margin-bottom: 12px;"
          ),
          plotOutput("sim_combo_density_plot", height = "320px"),
          div(
            downloadButton("dl_sim_combo_density_png", "Download simulated posterior overlay (PNG)"),
            style = "margin-bottom: 12px;"
          ),
          p(class = "small", style = "color: var(--muted);",
            "Only fits from the current shared spatial simulator are combined here."),
          downloadButton("dl_sim_uscr_csv", "Download available shared-simulation posterior summaries (CSV)")
        )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------

server <- function(input, output, session) {
  
  # ================================================================
  # SIMULATION: grid + inputs for models
  # ================================================================
  
  sim <- eventReactive(input$run_sim, {
    simulate_camera_counts(
      n_side    = input$n_side,
      spacing_m = input$spacing,
      days      = input$days,
      D_per_km2 = input$Dtrue,
      lambda0   = input$lambda0,
      home_range_km2 = input$home_range_km2,
      buffer_m  = input$sim_buffer_m,
      seed      = input$seed
    )
  }, ignoreInit = TRUE)
  
  sim_data <- reactive({
    req(sim())
    sim_model_inputs(
      sim = sim(),
      detection_radius_m = input$r_m_sim
    )
  })
  
  observeEvent(input$run_sim, {
    sim_id <- paste0("shared_sim_", format(Sys.time(), "%Y%m%d%H%M%S"))
    current_shared_sim_id(sim_id)
    uscr_sim_fit(NULL)
    rem_sim_fit(NULL)
    tte_sim_fit(NULL)
    uscr_sim_dataset_id(NULL)
    rem_sim_dataset_id(NULL)
    tte_sim_dataset_id(NULL)
  }, ignoreInit = TRUE)
  
  output$grid_plot <- renderPlot({
    req(sim())
    traps_df <- as.data.frame(secr::traps(sim()$ch)) %>%
      mutate(camera = row_number())
    mask_df <- as.data.frame(sim()$mask)
    
    ggplot() +
      geom_point(data = mask_df, aes(x = x, y = y), alpha = 0.1, color = redwood_colors[5]) +
      geom_point(data = traps_df, aes(x = x, y = y), size = 3, color = redwood_colors[3]) +
      geom_label(
        data = traps_df,
        aes(x = x, y = y, label = camera),
        nudge_y = 18,
        size = 3.8,
        color = redwood_colors[1],
        fill = "white",
        label.size = 0.2,
        label.padding = grid::unit(0.12, "lines")
      ) +
      coord_equal() +
      theme_minimal() +
      theme(axis.text = element_text(size = 11)) +
      labs(
        x = "m", y = "m",
        title = sprintf(
          "%d cameras, spacing %dm (simulated density: %.1f animals/km²)",
          sim()$truth$n_cams,
          sim()$truth$spacing_m,
          sim()$truth$D_per_km2
        )
      )
  })
  
  output$sim_deer_distribution_plot <- renderPlot({
    req(sim_data())
    d <- sim_data()
    
    # Combine camera locations with counts
    deer_dist <- d$out %>%
      mutate(
        x = utm_e * 1000,  # Convert km to m for consistency with grid plot
        y = utm_n * 1000,
        total_animals = d$camera_counts
      )
    
    td_max <- max(deer_dist$total_animals, na.rm = TRUE)
    size_breaks <- pretty(c(0, max(1, td_max)), n = 4)
    size_breaks <- unique(size_breaks[size_breaks >= 0 & size_breaks <= max(1, td_max)])
    ggplot(deer_dist, aes(x = x, y = y)) +
      geom_point(
        aes(size = total_animals),
        shape = 21,
        fill = redwood_colors[3],
        color = redwood_colors[1],
        stroke = 0.35,
        alpha = 0.75
      ) +
      ggrepel::geom_label_repel(
        aes(label = Site),
        size = 3,
        max.overlaps = Inf,
        fill = "white",
        label.size = 0.2,
        box.padding = 0.25,
        point.padding = 0.2
      ) +
      scale_size_area(
        name   = "Total detections",
        max_size = 14,
        breaks = size_breaks
      ) +
      coord_equal() +
      labs(
        title = "Total Animal Detections by Camera",
        x = "m",
        y = "m"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
  
  output$camera_table <- renderDT({
    req(sim())
    ch <- sim()$ch
    cams <- as.data.frame(secr::traps(ch)) %>%
      mutate(camera_id = paste0("C", row_number()))
    
    counts_mat <- get_counts_matrix(ch)
    totcounts  <- rowSums(counts_mat, na.rm = TRUE)
    
    cams %>%
      mutate(total_counts = totcounts) %>%
      select(camera_id, x, y, total_counts) %>%
      datatable(options = list(pageLength = 8))
  })
  
  # ================================================================
  # CAMERA ARRAY DESIGN
  # ================================================================

  camera_array_result <- reactiveVal(NULL)
  camera_array_log <- reactiveVal("No camera array has been generated yet.")

  observeEvent(input$reset_camera_array, {
    shinyjs::reset("camera_array_form")
    camera_array_result(NULL)
    camera_array_log("Camera-array inputs cleared. Upload a boundary to start a new design.")
  })

  observeEvent(input$run_camera_array, {
    req(input$camera_boundary_files)

    log_lines <- character()
    add_log <- function(...) {
      log_lines <<- c(log_lines, paste0(...))
      camera_array_log(paste(log_lines, collapse = "\n"))
    }

    camera_array_result(NULL)
    camera_array_log("Starting camera-array generation...")

    showNotification(
      "Reading uploaded spatial layers for the camera array...",
      type = "message",
      duration = NULL,
      id = "camera_array_status"
    )

    res <- tryCatch(
      {
        add_log("Reading boundary file(s)...")
        boundary_sf <- read_uploaded_spatial_group(
          input$camera_boundary_files,
          label = "boundary",
          required = TRUE
        )

        showNotification(
          "Reading optional roads, trails, buildings, parking, and exclusion layers...",
          type = "message",
          duration = NULL,
          id = "camera_array_status"
        )
        roads_sf <- read_uploaded_spatial_group(input$camera_roads_files, "roads layer")
        trails_sf <- read_uploaded_spatial_group(input$camera_trails_files, "trails layer")
        buildings_sf <- read_uploaded_spatial_group(input$camera_buildings_files, "buildings layer")
        parking_sf <- read_uploaded_spatial_group(input$camera_parking_files, "parking layer")
        exclusion_sf <- read_uploaded_spatial_group(input$camera_exclusion_files, "exclusion layer")

        add_log("Boundary loaded successfully.")
        optional_counts <- c(
          roads = if (is.null(roads_sf)) 0 else nrow(roads_sf),
          trails = if (is.null(trails_sf)) 0 else nrow(trails_sf),
          buildings = if (is.null(buildings_sf)) 0 else nrow(buildings_sf),
          parking = if (is.null(parking_sf)) 0 else nrow(parking_sf),
          exclusions = if (is.null(exclusion_sf)) 0 else nrow(exclusion_sf)
        )
        add_log(
          "Optional layer features loaded: ",
          paste(names(optional_counts), optional_counts, sep = "=", collapse = ", ")
        )

        showNotification(
          "Building the allowable camera-design area...",
          type = "message",
          duration = NULL,
          id = "camera_array_status"
        )
        design_area <- build_camera_design_area(
          boundary_sf = boundary_sf,
          roads = roads_sf,
          trails = trails_sf,
          buildings = buildings_sf,
          parking_lots = parking_sf,
          exclusion_areas = exclusion_sf,
          buffer_dist_m = input$camera_buffer_m
        )

        boundary_area_ha <- as.numeric(sum(sf::st_area(design_area$boundary))) / 10000
        allowed_area_ha <- as.numeric(sum(sf::st_area(design_area$allowed_area))) / 10000
        allowed_area_m2 <- as.numeric(sum(sf::st_area(design_area$allowed_area)))
        excluded_area_ha <- if (!is.null(design_area$excluded_area) && nrow(design_area$excluded_area)) {
          as.numeric(sum(sf::st_area(design_area$excluded_area))) / 10000
        } else {
          0
        }

        add_log("Projected design area created.")
        add_log("Boundary area: ", format_num(boundary_area_ha, 1), " ha")
        add_log("Allowable camera area: ", format_num(allowed_area_ha, 1), " ha")
        if (excluded_area_ha > 0) {
          add_log("Excluded / buffered area removed: ", format_num(excluded_area_ha, 1), " ha")
        }

        showNotification(
          "Generating candidate camera locations...",
          type = "message",
          duration = NULL,
          id = "camera_array_status"
        )
        prefix <- trimws(input$camera_prefix)
        if (!nzchar(prefix)) prefix <- "CAM"

        candidates <- generate_camera_candidates(
          allowed_area = design_area$allowed_area,
          spacing_m = input$camera_spacing_m,
          camera_prefix = prefix
        )
        selected <- select_camera_subset(
          candidates_sf = candidates,
          budget = input$camera_budget,
          n_alternates = input$camera_alternates,
          spacing_m = input$camera_spacing_m
        )

        add_log("Suggested camera candidates: ", nrow(candidates))
        add_log(
          "Final cameras retained: ",
          sum(selected$status == "final"),
          if (any(selected$status == "alternate")) {
            paste0(" (plus ", sum(selected$status == "alternate"), " alternates)")
          } else {
            ""
          }
        )
        approx_spacing_needed_m <- if (!is.na(input$camera_budget) &&
                                       input$camera_budget > 0 &&
                                       allowed_area_m2 > 0) {
          sqrt(allowed_area_m2 / input$camera_budget)
        } else {
          NA_real_
        }
        if (is.finite(approx_spacing_needed_m) && input$camera_budget < nrow(candidates)) {
          add_log(
            "Approximate even spacing supported by ",
            input$camera_budget,
            " final cameras across the allowable area: ",
            format_num(approx_spacing_needed_m, 0),
            " m"
          )
          if (approx_spacing_needed_m > input$camera_spacing_m) {
            add_log(
              "NOTE: keeping only ",
              input$camera_budget,
              " final cameras cannot preserve a maximum adjacent spacing of ",
              input$camera_spacing_m,
              " m across the full allowable area. A spacing closer to ",
              format_num(approx_spacing_needed_m, 0),
              " m would be needed."
            )
          }
        }
        budget_note <- NA_character_
        if (!is.na(input$camera_budget) && input$camera_budget > nrow(candidates) && nrow(candidates) > 0) {
          budget_note <- paste0(
            "The current spacing generated ",
            nrow(candidates),
            " candidate cameras. To place more cameras than that, rerun with a smaller spacing value."
          )
          add_log(
            "NOTE: ", budget_note
          )
        }
        spacing_note <- attr(selected, "spacing_note", exact = TRUE)
        if (!is.null(spacing_note) && nzchar(spacing_note)) {
          add_log("NOTE: ", spacing_note)
          if (is.finite(approx_spacing_needed_m)) {
            add_log(
              "Tip: if you want the array to be generated directly at that coarser spacing, try rerunning with spacing near ",
              format_num(approx_spacing_needed_m, 0),
              " m and leave the camera-count field at 0."
            )
          }
        }

        showNotification(
          "Preparing camera-array exports and map outputs...",
          type = "message",
          duration = NULL,
          id = "camera_array_status"
        )
        candidates_export <- camera_points_for_export(candidates)
        selected_export <- camera_points_for_export(selected)

        list(
          boundary = design_area$boundary,
          allowed_area = design_area$allowed_area,
          excluded_area = design_area$excluded_area,
          candidates = candidates_export,
          selected = selected_export,
          summary = list(
            spacing_m = input$camera_spacing_m,
            buffer_m = input$camera_buffer_m,
            budget = input$camera_budget,
            alternates = input$camera_alternates,
            boundary_area_ha = boundary_area_ha,
            allowed_area_ha = allowed_area_ha,
            excluded_area_ha = excluded_area_ha,
            approx_spacing_needed_m = approx_spacing_needed_m,
            budget_note = budget_note,
            spacing_note = spacing_note %||% NA_character_,
            max_nn_m = attr(selected, "max_nearest_neighbor_m", exact = TRUE) %||% NA_real_,
            mean_nn_m = attr(selected, "mean_nearest_neighbor_m", exact = TRUE) %||% NA_real_
          )
        )
      },
      error = function(e) {
        removeNotification("camera_array_status")
        add_log("ERROR: ", e$message)
        showNotification(
          paste("Camera array design failed:", e$message),
          type = "error",
          duration = NULL
        )
        NULL
      }
    )

    if (!is.null(res)) {
      removeNotification("camera_array_status")
      camera_array_result(res)
      add_log("Camera-array design completed successfully.")
      showNotification("Camera array generated successfully.", type = "message", duration = 6)
    }
  })

  output$camera_array_log <- renderText({
    camera_array_log()
  })

  output$camera_array_summary <- renderText({
    res <- camera_array_result()
    if (is.null(res)) {
      return("No camera array has been generated yet.")
    }

    summary <- res$summary
    selected <- res$selected

    paste(
      c(
        paste("Suggested spacing:", summary$spacing_m, "m"),
        paste("Buffer around optional layers:", summary$buffer_m, "m"),
        paste("Boundary area:", format_num(summary$boundary_area_ha, 1), "ha"),
        paste("Allowable area:", format_num(summary$allowed_area_ha, 1), "ha"),
        paste("Candidates generated:", nrow(res$candidates)),
        paste("Final cameras:", sum(selected$status == "final")),
        paste("Alternate cameras:", sum(selected$status == "alternate")),
        if (is.finite(summary$approx_spacing_needed_m)) paste("Approximate even spacing supported by final camera count:", format_num(summary$approx_spacing_needed_m, 0), "m"),
        if (!is.na(summary$budget_note) && nzchar(summary$budget_note)) paste("Budget note:", summary$budget_note),
        if (is.finite(summary$mean_nn_m)) paste("Mean nearest-neighbor spacing (final cameras):", format_num(summary$mean_nn_m, 0), "m"),
        if (is.finite(summary$max_nn_m)) paste("Largest nearest-neighbor spacing (final cameras):", format_num(summary$max_nn_m, 0), "m"),
        if (!is.na(summary$spacing_note) && nzchar(summary$spacing_note)) paste("Spacing note:", summary$spacing_note)
      ),
      collapse = "\n"
    )
  })

  output$camera_array_map <- renderLeaflet({
    res <- camera_array_result()
    req(res)

    boundary_ll <- sf::st_transform(res$boundary, 4326)
    allowed_ll <- sf::st_transform(res$allowed_area, 4326)
    excluded_ll <- if (!is.null(res$excluded_area) && nrow(res$excluded_area)) {
      sf::st_transform(res$excluded_area, 4326)
    } else {
      NULL
    }
    candidates_ll <- sf::st_transform(res$candidates, 4326)
    selected_ll <- sf::st_transform(res$selected, 4326)

    m <- leaflet() |>
      addTiles() |>
      addPolygons(
        data = boundary_ll,
        color = redwood_colors[1],
        weight = 2,
        fill = FALSE,
        group = "Boundary"
      ) |>
      addPolygons(
        data = allowed_ll,
        color = redwood_colors[3],
        weight = 1,
        fillColor = redwood_colors[4],
        fillOpacity = 0.2,
        group = "Allowed area"
      ) |>
      addCircleMarkers(
        data = candidates_ll,
        radius = 3,
        color = "#666666",
        stroke = FALSE,
        fillOpacity = 0.7,
        group = "Candidates",
        popup = ~paste0("<strong>", camera_id, "</strong>")
      )

    if (!is.null(excluded_ll)) {
      m <- m |>
        addPolygons(
          data = excluded_ll,
          color = "#8b5a2b",
          weight = 1,
          fillColor = "#c59d74",
          fillOpacity = 0.3,
          group = "Excluded area"
        )
    }

    final_pts <- selected_ll[selected_ll$status == "final", , drop = FALSE]
    alt_pts <- selected_ll[selected_ll$status == "alternate", , drop = FALSE]

    if (nrow(final_pts)) {
      m <- m |>
        addCircleMarkers(
          data = final_pts,
          radius = 6,
          color = redwood_colors[3],
          fillColor = redwood_colors[3],
          fillOpacity = 1,
          weight = 1,
          group = "Final cameras",
          popup = ~paste0("<strong>", camera_id, "</strong><br>Status: final")
        )
    }

    if (nrow(alt_pts)) {
      m <- m |>
        addCircleMarkers(
          data = alt_pts,
          radius = 6,
          color = "#d97904",
          fillColor = "#f4b266",
          fillOpacity = 1,
          weight = 1,
          group = "Alternates",
          popup = ~paste0("<strong>", camera_id, "</strong><br>Status: alternate")
        )
    }

    m |>
      addLayersControl(
        overlayGroups = c("Boundary", "Allowed area", "Excluded area", "Candidates", "Final cameras", "Alternates"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })

  output$camera_array_table <- DT::renderDT({
    res <- camera_array_result()
    req(res)

    export_df <- res$selected |>
      sf::st_drop_geometry() |>
      dplyr::select(camera_id, candidate_camera_id, status, longitude, latitude, utm_easting_m, utm_northing_m) |>
      dplyr::rename(
        `Camera ID` = camera_id,
        `Original candidate ID` = candidate_camera_id,
        Status = status,
        Longitude = longitude,
        Latitude = latitude,
        `UTM Easting (m)` = utm_easting_m,
        `UTM Northing (m)` = utm_northing_m
      )

    DT::datatable(export_df, options = list(pageLength = 10, autoWidth = TRUE), rownames = FALSE)
  })

  output$download_camera_array_csv <- downloadHandler(
    filename = function() {
      paste0("camera_array_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      readr::write_csv(sf::st_drop_geometry(res$selected), file)
    }
  )

  output$download_camera_array_gpkg <- downloadHandler(
    filename = function() {
      paste0("camera_array_", Sys.Date(), ".gpkg")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      sf::st_write(res$selected, file, layer = "camera_points", delete_dsn = TRUE, quiet = TRUE)
    }
  )

  output$download_camera_area_gpkg <- downloadHandler(
    filename = function() {
      paste0("camera_design_area_", Sys.Date(), ".gpkg")
    },
    content = function(file) {
      res <- camera_array_result()
      req(res)
      sf::st_write(res$boundary, file, layer = "boundary", delete_dsn = TRUE, quiet = TRUE)
      sf::st_write(res$allowed_area, file, layer = "allowed_area", append = TRUE, quiet = TRUE)
      if (!is.null(res$excluded_area) && nrow(res$excluded_area)) {
        sf::st_write(res$excluded_area, file, layer = "excluded_area", append = TRUE, quiet = TRUE)
      }
    }
  )
  
  # ================================================================
  # NPS QC: deployment + images
  # ================================================================
  
  deployment_checked <- reactiveVal(NULL)
  deployment_issues  <- reactiveVal(NULL)
  
  observeEvent(input$deployment_csv, {
    req(input$deployment_csv)
    
    deployment_raw <- readr::read_csv(
      input$deployment_csv$datapath,
      show_col_types = FALSE
    )
    
    deployment_raw <- clean_deployment_import(deployment_raw)
    
    msgs  <- character()
    warns <- character()
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            check_deployments(deployment_raw)
          },
          message = function(m) {
            msgs  <<- c(msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            warns <<- c(warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        deployment_checked(NULL)
        deployment_issues(list(
          messages = msgs,
          warnings = c(warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (!is.null(res)) {
      deployment_checked(res)
      deployment_issues(list(
        messages = msgs,
        warnings = warns
      ))
    }
  })
  
  output$deployment_check_log <- renderText({
    issues <- deployment_issues()
    req(issues)
    
    lines <- c()
    if (length(issues$messages)) {
      lines <- c(lines, "Messages:",
                 paste0("  - ", trimws(issues$messages)), "")
    }
    if (length(issues$warnings)) {
      lines <- c(lines, "Warnings / issues:",
                 paste0("  - ", trimws(issues$warnings)))
    }
    if (!length(lines)) lines <- "No issues reported."
    paste(lines, collapse = "\n")
  })
  
  output$deployment_preview <- DT::renderDT({
    req(deployment_checked())
    DT::datatable(deployment_checked(), options = list(pageLength = 10))
  })
  
  output$download_deployment_checked <- downloadHandler(
    filename = function() {
      base <- tools::file_path_sans_ext(input$deployment_csv$name)
      paste0(base, "_CHECKED.csv")
    },
    content = function(file) {
      req(deployment_checked())
      readr::write_csv(deployment_checked(), file)
    }
  )
  
  # ---- Images QC ----
  
  images_checked <- reactiveVal(NULL)
  images_issues  <- reactiveVal(NULL)
  
  observeEvent(input$images_csv, {
    req(input$images_csv)
    
    if (is.null(deployment_checked())) {
      showNotification(
        "Please upload and check a deployment CSV before checking images.",
        type = "error"
      )
      return(NULL)
    }
    
    images_raw <- readr::read_csv(
      input$images_csv$datapath,
      show_col_types = FALSE,
      col_types = readr::cols(
        `Sighting Count` = readr::col_character(),
        Species = readr::col_character()
      )
    )
    
    images_raw <- clean_images_import(images_raw)
    
    # Re-run deployment QC with images (malfunction date rules use site overlap)
    dep_msgs  <- character()
    dep_warns <- character()
    dep_res <- tryCatch(
      {
        withCallingHandlers(
          {
            check_deployments(deployment_checked(), images_raw)
          },
          message = function(m) {
            dep_msgs <<- c(dep_msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            dep_warns <<- c(dep_warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
      },
      error = function(e) {
        showNotification(
          paste("Deployment re-check failed:", e$message),
          type = "error"
        )
        images_checked(NULL)
        images_issues(list(
          messages = dep_msgs,
          warnings = c(dep_warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (is.null(dep_res)) return(NULL)
    
    deployment_checked(dep_res)
    deployment_issues(list(messages = dep_msgs, warnings = dep_warns))
    
    msgs  <- dep_msgs
    warns <- dep_warns
    
    checked <- NULL
    
    res <- tryCatch(
      {
        withCallingHandlers(
          {
            checked <- check_images(
              images      = images_raw,
              deployments = dep_res
            )
          },
          message = function(m) {
            msgs  <<- c(msgs, m$message)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            warns <<- c(warns, w$message)
            invokeRestart("muffleWarning")
          }
        )
        checked
      },
      error = function(e) {
        images_checked(NULL)
        images_issues(list(
          messages = msgs,
          warnings = c(warns, paste0("ERROR: ", e$message))
        ))
        return(NULL)
      }
    )
    
    if (is.null(res)) return(NULL)
    
    if (isTRUE(input$apply_56day_trim)) {
      trim_msgs <- character()
      trim_days <- input$trim_days
      if (is.null(trim_days) || is.na(trim_days)) trim_days <- 56L
      trim_days <- max(1L, as.integer(trim_days))
      
      trimmed <- withCallingHandlers(
        {
          trim_images_to_days(res, max_days = trim_days)
        },
        message = function(m) {
          trim_msgs <<- c(trim_msgs, m$message)
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          warns <<- c(warns, w$message)
          invokeRestart("muffleWarning")
        }
      )
      
      images_checked(trimmed)
      images_issues(list(
        messages = c(msgs, trim_msgs),
        warnings = warns
      ))
    } else {
      images_checked(res)
      images_issues(list(
        messages = c(
          msgs,
          "Skipped the optional deployment-length trim; using all validated images that fall within the deployment windows."
        ),
        warnings = warns
      ))
    }
  })
  
  output$images_check_log <- renderText({
    issues <- images_issues()
    req(issues)
    
    lines <- c()
    if (length(issues$messages)) {
      lines <- c(lines, "Messages:",
                 paste0("  - ", trimws(issues$messages)), "")
    }
    if (length(issues$warnings)) {
      lines <- c(lines, "Warnings / issues:",
                 paste0("  - ", trimws(issues$warnings)))
    }
    if (!length(lines)) lines <- "No issues reported."
    paste(lines, collapse = "\n")
  })
  
  output$images_preview <- DT::renderDT({
    req(images_checked())
    DT::datatable(images_checked(), options = list(pageLength = 10))
  })
  
  output$download_images_checked <- downloadHandler(
    filename = function() {
      base <- tools::file_path_sans_ext(input$images_csv$name)
      if (isTRUE(input$apply_56day_trim)) {
        trim_days <- input$trim_days
        if (is.null(trim_days) || is.na(trim_days)) trim_days <- 56L
        paste0(base, "_CHECKED_TRIM", max(1L, as.integer(trim_days)), ".csv")
      } else {
        paste0(base, "_CHECKED.csv")
      }
    },
    content = function(file) {
      req(images_checked())
      readr::write_csv(images_checked(), file)
    }
  )
  
  # ================================================================
  # DATA SUMMARY (uploaded data)
  # ================================================================
  
  deploy_summary <- reactive({
    req(deployment_checked())
    summarize_deployments(deployment_checked())
  })
  
  images_standardized <- reactive({
    req(images_checked())
    standardize_deer_species(images_checked())
  })
  
  image_summary <- reactive({
    summarize_images_by_species(images_standardized())
  })
  
  available_species <- reactive({
    req(images_standardized())
    species <- sort(unique(images_standardized()$Species))
    species <- species[!is.na(species) & nzchar(species)]
    species
  })
  
  observe({
    species <- available_species()
    current <- isolate(input$summary_species)
    selected <- if (!is.null(current) && current %in% species) {
      current
    } else if ("Deer" %in% species) {
      "Deer"
    } else if (length(species) > 0) {
      species[[1]]
    } else {
      character(0)
    }
    
    updateSelectInput(
      session,
      "summary_species",
      choices = species,
      selected = selected
    )
  })
  
  species_objects <- reactive({
    req(images_standardized(), input$summary_species)
    
    list(
      species_summary = species_summary_per_site(images_standardized(), input$summary_species),
      species_counts  = species_counts_per_camera(
        images_standardized(),
        input$summary_species,
        deployments = deployment_checked()
      ),
      daily_species   = species_daily_detections(images_standardized(), input$summary_species)
    )
  })
  
  output$deploy_summary_table <- DT::renderDT({
    req(deploy_summary())
    summary_df <- deploy_summary() |>
      dplyr::rename(
        `Number of Cameras` = n_cameras,
        `Mean Number of Days Deployed` = mean_days,
        `Total Number of Days Deployed` = total_camera_days
      )
    DT::datatable(
      summary_df,
      options = list(dom = "t", paging = FALSE, autoWidth = TRUE),
      rownames = FALSE
    )
  })
  
  output$camera_map <- leaflet::renderLeaflet({
    req(deployment_checked())
    dep <- deployment_checked()
    
    validate(
      need(all(c("Longitude", "Latitude", "Site Name") %in% names(dep)),
           "Deployment file must contain Longitude, Latitude, and Site Name.")
    )
    
    icon.fa <- leaflet::makeAwesomeIcon(
      icon = "camera", markerColor = "green",
      library = "fa", iconColor = "black"
    )
    
    leaflet::leaflet(dep) %>%
      leaflet::addTiles() %>%
      leaflet::addAwesomeMarkers(
        lng   = ~Longitude,
        lat   = ~Latitude,
        icon  = icon.fa,
        popup = ~as.character(`Site Name`),
        label = ~as.character(`Site Name`)
      )
  })
  
  output$image_summary_table <- DT::renderDT({
    req(image_summary())
    summary_df <- image_summary() |>
      dplyr::rename(
        `Total Images` = total_images,
        `Total Detections` = total_detections,
        `Cameras with Detections` = cameras_detected
      )
    DT::datatable(
      summary_df,
      options = list(pageLength = 10, autoWidth = TRUE, dom = "tp"),
      rownames = FALSE
    )
  })
  
  output$species_bar_plot <- renderPlot({
    df <- image_summary()
    validate(
      need(nrow(df) > 0, "No images available to summarize.")
    )
    
    ggplot(df, aes(x = reorder(Species, -total_detections),
                   y = total_detections)) +
      geom_bar(stat = "identity", fill = redwood_colors[3]) +
      labs(x = "Species", y = "Total detections") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$deer_summary_table <- DT::renderDT({
    species_name <- req(input$summary_species)
    species_df <- species_objects()$species_summary
    validate(
      need(nrow(species_df) > 0, paste0("No '", species_name, "' images found in dataset."))
    )
    species_df <- species_df |>
      dplyr::rename(
        `Site Name` = `Site Name`,
        `Number of Images` = images,
        `Number of Detections` = detections
      )
    DT::datatable(
      species_df,
      options = list(pageLength = 10, autoWidth = TRUE, dom = "tp"),
      rownames = FALSE
    )
  })
  
  output$deer_bubble_plot <- renderPlot({
    species_name <- req(input$summary_species)
    dc <- species_objects()$species_counts
    validate(
      need(nrow(dc) > 0, paste0("No '", species_name, "' images found in dataset."))
    )

    dc <- dc |>
      dplyr::mutate(
        Longitude = suppressWarnings(as.numeric(Longitude)),
        Latitude = suppressWarnings(as.numeric(Latitude))
      ) |>
      dplyr::filter(is.finite(Longitude), is.finite(Latitude))

    validate(
      need(
        nrow(dc) > 0,
        paste0(
          "No valid site coordinates are available for the selected species. ",
          "Check the deployment file Latitude/Longitude columns."
        )
      )
    )

    td_max <- max(dc$total_detections, na.rm = TRUE)
    ggplot(dc, aes(x = Longitude, y = Latitude)) +
      geom_point(aes(size = total_detections),
                 color = redwood_colors[3], alpha = 0.7) +
      ggrepel::geom_text_repel(
        aes(label = `Site Name`),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_size_continuous(
        name   = paste("Total", species_name),
        range  = c(2, 14),
        limits = c(0, NA),
        breaks = pretty(c(0, max(1, td_max)))
      ) +
      coord_fixed() +
      labs(
        title = paste("Total", species_name, "detections by site"),
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })
  
  output$deer_daily_plot <- renderPlot({
    species_name <- req(input$summary_species)
    dd <- species_objects()$daily_species
    validate(
      need(nrow(dd) > 0, paste0("No '", species_name, "' images found in dataset."))
    )
    
    ggplot(dd, aes(x = Date, y = `Site Name`, size = detections)) +
      geom_point(
        color = redwood_colors[3], alpha = 0.7,
        position = position_jitter(width = 0, height = 0.1)
      ) +
      scale_size_continuous(name = paste(species_name, "count"), range = c(2, 8)) +
      labs(
        title = paste("Daily", species_name, "detections by site"),
        x = "Date",
        y = "Site name"
      ) +
      theme_minimal() +
      theme(
        axis.text.y     = element_text(size = 8),
        legend.position = "right"
      )
  })
  
  # ================================================================
  # MODEL FIT OBJECTS (reactiveVal) – SIM + NPS PER MODEL
  # ================================================================
  
  uscr_sim_fit <- reactiveVal(NULL)
  rem_sim_fit  <- reactiveVal(NULL)
  tte_sim_fit  <- reactiveVal(NULL)
  current_shared_sim_id <- reactiveVal(NULL)
  uscr_sim_dataset_id <- reactiveVal(NULL)
  rem_sim_dataset_id  <- reactiveVal(NULL)
  tte_sim_dataset_id  <- reactiveVal(NULL)

  uscr_nps_fit <- reactiveVal(NULL)
  rem_nps_fit  <- reactiveVal(NULL)
  tte_nps_fit  <- reactiveVal(NULL)
  
  make_model_debug <- function(model, data_source, supported_sources, preprocess_note) {
    list(
      model = model,
      data_source = data_source,
      supported_sources = supported_sources,
      preprocess_note = preprocess_note,
      status = "idle",
      stage = "Not started",
      started_at = NULL,
      finished_at = NULL,
      guidance = "No run has been started yet.",
      raw_error = NULL,
      context = character(),
      pause_queued = FALSE,
      stop_queued = FALSE,
      history = character()
    )
  }
  
  update_model_debug <- function(rv,
                                 status = NULL,
                                 stage = NULL,
                                 started_at = NULL,
                                 finished_at = NULL,
                                 guidance = NULL,
                                 raw_error = NULL,
                                 pause_queued = NULL,
                                 stop_queued = NULL,
                                 context = NULL,
                                 log_entry = NULL) {
    state <- rv()
    if (is.null(state)) state <- list()
    
    if (!is.null(status)) state$status <- status
    if (!is.null(stage)) state$stage <- stage
    if (!is.null(started_at)) state$started_at <- started_at
    if (!is.null(finished_at)) state$finished_at <- finished_at
    if (!is.null(guidance)) state$guidance <- guidance
    if (!is.null(raw_error)) state$raw_error <- raw_error
    if (!is.null(pause_queued)) state$pause_queued <- pause_queued
    if (!is.null(stop_queued)) state$stop_queued <- stop_queued
    if (!is.null(context)) state$context <- context
    if (!is.null(log_entry)) {
      state$history <- c(
        state$history,
        paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", log_entry)
      )
    }
    
    rv(state)
    invisible(state)
  }
  
  format_num <- function(x, digits = 2) {
    if (length(x) != 1 || is.na(x) || !is.finite(x)) return("NA")
    format(round(x, digits), nsmall = digits, trim = TRUE)
  }

  get_camera_angle_vector <- function(out, fallback_deg = 55) {
    fallback_deg <- as.numeric(fallback_deg[[1]])
    if (!is.finite(fallback_deg)) fallback_deg <- 55
    fallback_deg <- max(1, min(360, fallback_deg))

    if (!"Camera Detection Angle" %in% names(out)) {
      return(rep(fallback_deg, nrow(out)))
    }

    theta_raw <- suppressWarnings(as.numeric(out$`Camera Detection Angle`))
    theta_use <- theta_raw
    theta_use[is.na(theta_use) | theta_use <= 0 | theta_use > 360] <- fallback_deg
    theta_use
  }

  summarize_camera_angle_source <- function(out, fallback_deg = 55) {
    theta_use <- get_camera_angle_vector(out, fallback_deg = fallback_deg)
    has_column <- "Camera Detection Angle" %in% names(out)
    theta_raw <- if (has_column) suppressWarnings(as.numeric(out$`Camera Detection Angle`)) else numeric(0)
    n_fallback <- if (has_column) {
      sum(is.na(theta_raw) | theta_raw <= 0 | theta_raw > 360, na.rm = FALSE)
    } else {
      nrow(out)
    }

    range_text <- paste0(
      format_num(min(theta_use, na.rm = TRUE), 1),
      " to ",
      format_num(max(theta_use, na.rm = TRUE), 1),
      " degrees"
    )

    if (!has_column) {
      return(paste("Camera detection angle:", range_text, "(all cameras using fallback from Model settings)"))
    }

    if (n_fallback > 0) {
      return(paste(
        "Camera detection angle:",
        range_text,
        paste0("(", n_fallback, " camera(s) using the fallback angle from Model settings)")
      ))
    }

    paste("Camera detection angle:", range_text, "(from uploaded deployment data)")
  }
  
  summarize_uscr_context <- function(d, source_label = "Uploaded", sim_truth = NULL, fast_mode = FALSE) {
    effective_iter <- if (isTRUE(fast_mode)) min(as.integer(input$iter_uscr), 4000L) else as.integer(input$iter_uscr)
    effective_burnin <- if (isTRUE(fast_mode)) min(as.integer(input$burnin_uscr), 500L) else as.integer(input$burnin_uscr)
    effective_thin <- if (isTRUE(fast_mode)) max(1L, min(as.integer(input$thin_uscr), 4L)) else as.integer(input$thin_uscr)
    effective_M <- if (isTRUE(fast_mode)) min(as.integer(input$M_uscr), 200L) else as.integer(input$M_uscr)

    det_dist <- if ("Detection Distance" %in% names(d$out)) {
      paste0(
        format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
        " to ",
        format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
        " m"
      )
    } else {
      "Not available"
    }
    
    lines <- c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste("Detection distance range:", det_dist),
      paste(
        "Requested USCR full run:",
        paste0(
          "iter=", effective_iter,
          ", burn-in=", effective_burnin,
          ", thin=", effective_thin,
          ", chains=", app_n_chains(),
          ", M=", effective_M
        )
      ),
      paste("USCR state-space buffer:", format_num(input$uscr_buffer_m, 0), "m"),
      if (isTRUE(fast_mode)) {
        "USCR fast test mode is ON: adaptive reruns are capped and lighter run settings are used for faster exploratory runs."
      } else {
        "Adaptive tuning may start with shorter runs but now uses the same chain count shown above."
      },
      if (source_label == "Uploaded") {
        "Supported sources here: uploaded field data and the shared spatial simulator."
      } else {
        "Supported sources here: the shared spatial simulator and uploaded field data."
      },
      "Preprocessing: total animal events per camera, camera-days, and buffered spatial state space."
    )
    
    if (!is.null(sim_truth)) {
      lines <- c(
        lines,
        paste("True simulated density:", format_num(sim_truth$D_per_km2, 1), "animals/km^2"),
        if (!is.null(sim_truth$home_range_km2)) {
          paste("Mean simulated 95% home-range size:", format_num(sim_truth$home_range_km2, 2), "km^2")
        },
        if (!is.null(sim_truth$buffer_m)) {
          paste("Simulated state-space buffer:", format_num(sim_truth$buffer_m, 0), "m")
        }
      )
    }
    
    lines
  }
  
  summarize_rem_context <- function(d) {
    c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste(
        "Detection distance range:",
        paste0(
          format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " to ",
          format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " m"
        )
      ),
      summarize_camera_angle_source(d$out, fallback_deg = input$theta),
      paste(
        "Requested REM run:",
        paste0(
          "iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", app_n_chains()
        )
      ),
      "Supported sources here: uploaded field data only.",
      "Preprocessing: REM uses total animal events per camera and camera-days."
    )
  }
  
  summarize_tte_context <- function(d) {
    c(
      paste("Cameras:", nrow(d$out)),
      paste("Total animal events passed to TTE:", sum(d$camera_counts, na.rm = TRUE)),
      paste("Total camera-days:", format_num(sum(d$camera_days, na.rm = TRUE), 1)),
      paste(
        "Detection distance range:",
        paste0(
          format_num(min(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " to ",
          format_num(max(d$out$`Detection Distance`, na.rm = TRUE), 1),
          " m"
        )
      ),
      summarize_camera_angle_source(d$out, fallback_deg = input$theta),
      paste(
        "Requested TTE run:",
        paste0(
          "iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", app_n_chains()
        )
      ),
      "Supported sources here: uploaded field data only.",
      "Current app preprocessing: TTE is receiving total animal events per camera and camera-days."
    )
  }
  
  summarize_shared_sim_context <- function(d, model_label, sim_truth = NULL) {
    model_settings <- switch(
      model_label,
      "USCR" = NULL,
      "REM" = paste(
        "Requested REM run:",
        paste0(
          "iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", app_n_chains(),
          ", fallback angle=", input$theta,
          " degrees"
        )
      ),
      "TTE" = paste(
        "Requested TTE run:",
        paste0(
          "iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", app_n_chains(),
          ", fallback angle=", input$theta,
          " degrees"
        )
      ),
      NULL
    )

    c(
      paste("Shared simulator model fit:", model_label),
      summarize_uscr_context(d, source_label = "Simulated", sim_truth = sim_truth),
      if (!is.null(model_settings)) model_settings,
      "This fit comes from the shared spatial simulator and is eligible for simulated Compare & combine."
    )
  }
  
  friendly_model_error <- function(model, data_source, error_text) {
    tips <- c(
      paste0(model, " on ", data_source, " did not finish."),
      "Read the raw error below, then try the model-specific guidance."
    )
    
    if (grepl("unused arguments", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "This usually means an older function definition is still loaded. Restart R and launch this app from the current project folder.")
    }
    if (grepl("zero camera-days", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "One or more cameras ended up with zero effort days. Check Start Date/End Date, malfunction handling, and the optional deployment-length trimming choice.")
    }
    if (grepl("camera_counts and camera_days must have length nrow\\(out\\)", error_text)) {
      tips <- c(tips, "The camera table, counts, and camera-days vectors got out of sync. Re-upload the data and check the uploaded-data model-input step.")
    }
    if (grepl("cannot allocate|vector memory|std::bad_alloc", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "This looks like a memory issue. Try fewer chains, fewer iterations, or a smaller USCR M for a test run.")
    }
    if (grepl("Package '.*' is required|there is no package called", error_text, ignore.case = TRUE)) {
      tips <- c(tips, "A required R package is missing in this session. Install the package named in the raw error, restart R, and rerun.")
    }
    if (grepl("NA|NaN|Inf", error_text)) {
      tips <- c(tips, "Missing or invalid numeric values reached the model. Check detection distances, camera-days, timestamps, and coordinates in the upload summary tabs.")
    }
    
    if (identical(model, "USCR")) {
      tips <- c(tips, "USCR does short adaptive tuning runs before the full run, and the final run is checked against the same convergence target. If it is too slow for testing, reduce chains, iterations, or M.")
    }
    if (identical(model, "REM")) {
      tips <- c(tips, "REM expects per-camera animal events, camera-days, and detection distances.")
      tips <- c(tips, "REM now doubles iterations automatically until the maximum Rhat is acceptable, so longer runs can be expected when convergence is slow.")
    }
    if (identical(model, "TTE")) {
      tips <- c(tips, "TTE expects per-camera animal events, camera-days, and detection distances. The current build passes total animal events per camera to TTE.")
      tips <- c(tips, "TTE now doubles iterations automatically until the maximum Rhat is acceptable, so longer runs can be expected when convergence is slow.")
    }
    
    paste(tips, collapse = "\n")
  }

  make_encounter_status_callback <- function(rv, model_label, set_progress = NULL, notification_id = NULL) {
    function(stage, detail = NULL, value = NULL) {
      stage_label <- switch(
        stage,
        setup = "Preparing model code",
        tuning = "Adaptive tuning / convergence checks",
        final_run = "Final MCMC run",
        stage
      )
      progress_detail <- detail %||% paste(model_label, "stage:", stage_label)

      if (is.function(set_progress)) {
        if (!is.null(value)) {
          set_progress(value = value, detail = progress_detail)
        } else {
          set_progress(detail = progress_detail)
        }
      }

      if (!is.null(notification_id) && nzchar(notification_id)) {
        showNotification(
          progress_detail,
          type = "message",
          duration = NULL,
          id = notification_id
        )
      }

      update_model_debug(
        rv,
        stage = stage_label,
        log_entry = progress_detail
      )
    }
  }
  
  format_model_debug <- function(state, fit = NULL) {
    lines <- c(
      paste("Model:", state$model),
      paste("Data source:", state$data_source),
      paste("Supported in app:", state$supported_sources),
      paste("Status:", state$status),
      paste("Stage:", state$stage)
    )
    
    if (!is.null(state$started_at)) {
      lines <- c(lines, paste("Started:", format(state$started_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(state$finished_at)) {
      lines <- c(lines, paste("Finished:", format(state$finished_at, "%Y-%m-%d %H:%M:%S")))
    }
    if (!is.null(state$started_at)) {
      end_time <- if (identical(state$status, "running")) Sys.time() else state$finished_at %||% Sys.time()
      elapsed_min <- as.numeric(difftime(end_time, state$started_at, units = "mins"))
      lines <- c(lines, paste("Elapsed:", format_num(elapsed_min, 1), "minutes"))
    }

    if (isTRUE(state$pause_queued) && identical(state$status, "running")) {
      lines <- c(
        lines,
        "Pause request: queued",
        "The app will save a provisional USCR estimate after the current round finishes."
      )
    }
    if (isTRUE(state$stop_queued) && identical(state$status, "running")) {
      lines <- c(
        lines,
        "Stop request: queued",
        "The app will stop the USCR run after the current round finishes."
      )
    }
    
    lines <- c(lines, "", "Guidance:", state$guidance)
    
    if (length(state$context)) {
      lines <- c(lines, "", "Input summary:", paste0("- ", state$context))
    }
    
    if (!is.null(fit) && !is.null(fit$settings)) {
      setting_lines <- c(
        paste0("- Final iterations: ", fit$settings$iter),
        paste0("- Burn-in: ", fit$settings$burnin),
        paste0("- Thin: ", fit$settings$thin),
        paste0("- Chains: ", fit$settings$n_chains)
      )
      if (!is.null(fit$settings$M)) {
        setting_lines <- c(paste0("- Final M: ", fit$settings$M), setting_lines)
      }
      lines <- c(lines, "", "Run settings used:", setting_lines)
      if (!is.null(fit$final_rhat_max) && is.finite(fit$final_rhat_max)) {
        lines <- c(lines, paste0("- Final max Rhat: ", format_num(fit$final_rhat_max, 3)))
      }
    }

    if (!is.null(fit) && isTRUE(fit$provisional)) {
      lines <- c(
        lines,
        "",
        "Provisional result:",
        paste0("- Stage saved: ", fit$provisional_stage %||% "latest completed round"),
        paste0("- Reason: ", fit$provisional_reason %||% "Paused before a final converged fit was produced."),
        "- This paused fit is excluded from Compare & combine until you resume or rerun USCR."
      )
    }
    
    if (!is.null(fit) && !is.null(fit$tuning_history) && nrow(fit$tuning_history) > 0) {
      tuning_df <- fit$tuning_history
      tuning_lines <- apply(
        tuning_df,
        1,
        function(row) {
          parts <- c(
            paste0("round ", row[["round"]]),
            paste0("iter=", row[["niter"]]),
            paste0("thin=", row[["thin"]]),
            paste0("chains=", row[["n_chains"]]),
            paste0("rhat=", format_num(as.numeric(row[["rhat_max"]]), 3))
          )
          if ("M" %in% names(tuning_df)) {
            parts <- append(parts, paste0("M=", row[["M"]]), after = 1)
          }
          if ("M_too_small" %in% names(tuning_df)) {
            parts <- c(parts, paste0("M_too_small=", row[["M_too_small"]]))
          }
          paste(parts, collapse = ", ")
        }
      )
      lines <- c(lines, "", paste0(state$model, " tuning history:"), paste0("- ", tuning_lines))
    }

    if (!is.null(fit) && !is.null(fit$final_run_history) && nrow(fit$final_run_history) > 0) {
      final_df <- fit$final_run_history
      final_lines <- apply(
        final_df,
        1,
        function(row) {
          parts <- c(
            paste0("round ", row[["round"]]),
            paste0("iter=", row[["niter"]]),
            paste0("thin=", row[["thin"]]),
            paste0("chains=", row[["n_chains"]]),
            paste0("rhat=", format_num(as.numeric(row[["rhat_max"]]), 3))
          )
          if ("M" %in% names(final_df)) {
            parts <- append(parts, paste0("M=", row[["M"]]), after = 1)
          }
          if ("M_too_small" %in% names(final_df)) {
            parts <- c(parts, paste0("M_too_small=", row[["M_too_small"]]))
          }
          paste(parts, collapse = ", ")
        }
      )
      lines <- c(lines, "", paste0(state$model, " final-run checks:"), paste0("- ", final_lines))
    }
    
    if (!is.null(state$raw_error)) {
      lines <- c(lines, "", "Raw error:", state$raw_error)
    }
    
    if (length(state$history)) {
      lines <- c(lines, "", "Recent log:", paste0("- ", tail(state$history, 8)))
    }
    
    paste(lines, collapse = "\n")
  }
  
  `%||%` <- function(x, y) {
    if (is.null(x)) y else x
  }
  
  uscr_sim_debug <- reactiveVal(
    make_model_debug("USCR", "Simulated grid", "Simulated grid; uploaded field data", "Total animal events per camera, camera-days, and buffered state space.")
  )
  rem_sim_debug <- reactiveVal(
    make_model_debug("REM", "Shared spatial simulator", "Shared spatial simulator", "REM uses total animal events per camera and camera-days derived from the shared spatial simulator.")
  )
  tte_sim_debug <- reactiveVal(
    make_model_debug("TTE", "Shared spatial simulator", "Shared spatial simulator", "TTE uses total animal events per camera and camera-days derived from the shared spatial simulator.")
  )
  uscr_nps_debug <- reactiveVal(
    make_model_debug("USCR", "Uploaded field data", "Uploaded field data; simulated grid", "Total animal events per camera, camera-days, and buffered state space.")
  )
  rem_nps_debug <- reactiveVal(
    make_model_debug("REM", "Uploaded field data", "Uploaded field data only", "REM uses total animal events per camera and camera-days.")
  )
  tte_nps_debug <- reactiveVal(
    make_model_debug("TTE", "Uploaded field data", "Uploaded field data only", "Current build passes total animal events per camera and camera-days.")
  )
  
  # Status tracking for model runs
  uscr_sim_running <- reactiveVal(FALSE)
  rem_sim_running <- reactiveVal(FALSE)
  tte_sim_running <- reactiveVal(FALSE)
  uscr_nps_running <- reactiveVal(FALSE)
  rem_nps_running <- reactiveVal(FALSE)
  tte_nps_running <- reactiveVal(FALSE)
  rem_nps_background <- reactiveVal(FALSE)
  
  # Stop flags for interrupting model runs
  stop_uscr_sim <- reactiveVal(FALSE)
  stop_uscr_nps <- reactiveVal(FALSE)
  pause_uscr_nps <- reactiveVal(FALSE)
  stop_rem_nps <- reactiveVal(FALSE)
  stop_tte_nps <- reactiveVal(FALSE)
  uscr_nps_resume_settings <- reactiveVal(NULL)
  
  nps_model_inputs <- reactive({
    req(deployment_checked(), images_checked())
    trim_days <- NULL
    if (isTRUE(input$apply_56day_trim)) {
      trim_days <- input$trim_days
      if (is.null(trim_days) || is.na(trim_days)) trim_days <- 56L
      trim_days <- max(1L, as.integer(trim_days))
    }
    build_nps_model_inputs(
      deployment_checked(),
      images_checked(),
      max_days = trim_days
    )
  })
  
  # --- Stop handlers ---
  
  # Use a global environment to store stop flags that can be checked during execution
  stop_flags_env <- new.env()
  
  observeEvent(input$stop_uscr_sim, {
    stop_uscr_sim(TRUE)
    stop_flags_env$stop_uscr_sim <- TRUE
    showNotification("Stopping USCR (sim) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_uscr_nps, {
    stop_uscr_nps(TRUE)
    stop_flags_env$stop_uscr_nps <- TRUE
    update_model_debug(
      uscr_nps_debug,
      stop_queued = TRUE,
      guidance = paste(
        "Stop has been queued for uploaded-data USCR.",
        "The current USCR round still has to finish before the run can terminate cleanly."
      ),
      log_entry = "USCR stop requested; waiting for the current round to finish."
    )
    showNotification("Stopping USCR (uploaded data) run...", 
                     type = "warning", duration = 3)
  })

  refresh_uscr_pause_resume_button <- function() {
    if (isTRUE(uscr_nps_running())) {
      updateActionButton(session, "pause_resume_uscr_nps", label = "Pause after current round")
      shinyjs::enable("pause_resume_uscr_nps")
    } else if (!is.null(uscr_nps_resume_settings())) {
      updateActionButton(session, "pause_resume_uscr_nps", label = "Resume from paused settings")
      shinyjs::enable("pause_resume_uscr_nps")
    } else {
      updateActionButton(session, "pause_resume_uscr_nps", label = "Pause/Resume")
      shinyjs::disable("pause_resume_uscr_nps")
    }
  }

  observe({
    uscr_nps_running()
    uscr_nps_resume_settings()
    refresh_uscr_pause_resume_button()
  })

  observeEvent(input$pause_resume_uscr_nps, {
    if (isTRUE(uscr_nps_running())) {
      if (isTRUE(pause_uscr_nps())) {
        showNotification(
          "USCR pause is already queued for the end of the current round.",
          type = "warning",
          duration = 4
        )
        return(NULL)
      }
      pause_uscr_nps(TRUE)
      stop_flags_env$pause_uscr_nps <- TRUE
      update_model_debug(
        uscr_nps_debug,
        pause_queued = TRUE,
        guidance = paste(
          "Pause has been queued for uploaded-data USCR.",
          "The current USCR round still has to finish before a provisional estimate can be saved."
        ),
        stop_queued = FALSE,
        log_entry = "USCR pause requested; waiting for the current round to finish."
      )
      showNotification(
        "Pausing USCR (uploaded data) after the current round finishes...",
        type = "warning",
        duration = 5
      )
      return(NULL)
    }

    resume_override <- uscr_nps_resume_settings()
    if (is.null(resume_override)) {
      showNotification(
        "No paused USCR round is available to resume from yet.",
        type = "warning",
        duration = 5
      )
      return(NULL)
    }

    start_uscr_nps_run(resume_override = resume_override, resumed_from_pause = TRUE)
  })
  
  observeEvent(input$stop_rem_nps, {
    if (isTRUE(rem_nps_background())) {
      update_model_debug(
        rem_nps_debug,
        status = "running",
        stage = "Background run in progress",
        guidance = "This REM run is already running in a background worker. It will keep running until it finishes or fails; server-side cancellation is not available yet for this path.",
        log_entry = "Stop requested, but background REM cancellation is not available."
      )
      showNotification(
        "REM background jobs cannot be cancelled after launch yet. This run will keep going in the background.",
        type = "warning",
        duration = 8
      )
      return(NULL)
    }
    stop_rem_nps(TRUE)
    stop_flags_env$stop_rem_nps <- TRUE
    showNotification("Stopping REM (uploaded data) run...", 
                     type = "warning", duration = 3)
  })
  
  observeEvent(input$stop_tte_nps, {
    stop_tte_nps(TRUE)
    stop_flags_env$stop_tte_nps <- TRUE
    showNotification("Stopping TTE (uploaded data) run...", 
                     type = "warning", duration = 3)
  })
  
  # Helper function to check stop flag and throw error if set
  check_stop_flag <- function(flag_name) {
    if (exists(flag_name, envir = stop_flags_env) && 
        isTRUE(get(flag_name, envir = stop_flags_env))) {
      err <- simpleError("Model execution stopped by user")
      class(err) <- c("model_stop_requested", class(err))
      stop(err)
    }
  }

  check_interrupt_action <- function(stop_flag_name = NULL, pause_flag_name = NULL) {
    if (!is.null(stop_flag_name) &&
        exists(stop_flag_name, envir = stop_flags_env) &&
        isTRUE(get(stop_flag_name, envir = stop_flags_env))) {
      return("stop")
    }
    if (!is.null(pause_flag_name) &&
        exists(pause_flag_name, envir = stop_flags_env) &&
        isTRUE(get(pause_flag_name, envir = stop_flags_env))) {
      return("pause")
    }
    NULL
  }

  app_n_chains <- function() {
    max(2L, as.integer(input$n_chains))
  }
  
  uscr_run_args <- function(waic = TRUE, status_callback = NULL, fast_mode = FALSE,
                            override = NULL, interrupt_callback = NULL) {
    args <- list(
      iter = input$iter_uscr,
      burnin = input$burnin_uscr,
      thin = input$thin_uscr,
      n_chains = app_n_chains(),
      M = input$M_uscr,
      log_sigma_mean = input$log_sigma_mean,
      log_sigma_sd = input$log_sigma_sd,
      buffer_m = input$uscr_buffer_m,
      log_lam0_mean = input$log_lam0_mean,
      log_lam0_sd = input$log_lam0_sd,
      sd_eps_shape = input$sd_eps_shape,
      sd_eps_rate = input$sd_eps_rate,
      adaptive = TRUE,
      compute_WAIC = waic,
      diagnostic_mode = FALSE,
      tuning_n_chains = app_n_chains(),
      parallel_chains = isTRUE(app_n_chains() > 1),
      interrupt_callback = interrupt_callback,
      status_callback = status_callback,
      verbose = FALSE
    )

    if (isTRUE(fast_mode)) {
      args$iter <- min(as.integer(input$iter_uscr), 4000L)
      args$burnin <- min(as.integer(input$burnin_uscr), 500L)
      args$thin <- max(1L, min(as.integer(input$thin_uscr), 4L))
      args$M <- min(as.integer(input$M_uscr), 200L)
      args$iter_tune <- max(args$burnin + 500L, 1500L)
      args$thin_tune <- 1L
      args$max_adapt_rounds <- 3L
      args$iter_cap <- 4000L
      args$M_cap <- 200L
    }

    if (!is.null(override)) {
      if (!is.null(override$iter)) args$iter <- max(as.integer(args$iter), as.integer(override$iter))
      if (!is.null(override$thin)) args$thin <- max(as.integer(args$thin), as.integer(override$thin))
      if (!is.null(override$M)) args$M <- max(as.integer(args$M), as.integer(override$M))
      if (!is.null(override$burnin)) args$burnin <- max(as.integer(args$burnin), as.integer(override$burnin))
    }
    
    # Keep the app compatible with whichever run_USCR() definition is
    # currently loaded, including older variants without these extras.
    valid_names <- names(formals(run_USCR_app))
    args[names(args) %in% valid_names]
  }

  # --- USCR: simulated ---
  
  observeEvent(input$run_uscr_sim, {
    req(sim_data())
    d <- sim_data()
    uscr_sim_debug(make_model_debug(
      "USCR",
      "Simulated grid",
      "Simulated grid; uploaded field data",
      "Total animal events per camera, camera-days, and buffered state space."
    ))
    update_model_debug(
      uscr_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "USCR supports simulated and uploaded field data. Watch this panel for tuning rounds and the final run stage.",
      raw_error = NULL,
      context = summarize_uscr_context(d, source_label = "Simulated", sim_truth = sim()$truth),
      log_entry = "Simulated USCR run requested."
    )
    
    # Reset stop flag
    stop_uscr_sim(FALSE)
    stop_flags_env$stop_uscr_sim <- FALSE
    uscr_sim_running(TRUE)
    uscr_sim_fit(NULL)  # Clear previous results
    uscr_sim_dataset_id(current_shared_sim_id())
    
    showNotification(
      "Running USCR on simulated data.",
      type = "message",
      duration = NULL,
      id = "uscr_sim_status"
    )
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        check_stop_flag("stop_uscr_sim")
        
        update_model_debug(
          uscr_sim_debug,
          stage = "Preparing run",
          log_entry = "Input checks passed. Starting USCR setup."
        )
        
        # Check stop flag again before running
        check_stop_flag("stop_uscr_sim")
        
        status_callback <- function(stage, detail = NULL, value = NULL) {
          check_stop_flag("stop_uscr_sim")
          stage_label <- switch(
            stage,
            setup = "Preparing state space and model code",
            tuning = "Adaptive tuning",
            final_run = "Final MCMC run / convergence check",
            stage
          )
          progress_detail <- detail %||% paste("USCR stage:", stage_label)
          update_model_debug(
            uscr_sim_debug,
            status = "running",
            stage = stage_label,
            log_entry = progress_detail
          )
          showNotification(
            progress_detail,
            type = "message",
            duration = NULL,
            id = "uscr_sim_status"
          )
        }
        
        do.call(
          run_USCR_app,
          c(
            list(
              out = d$out,
              camera_counts = d$camera_counts,
              camera_days = d$camera_days
            ),
            uscr_run_args(waic = FALSE, status_callback = status_callback)
          )
        )
      },
      error = function(e) {
        removeNotification("uscr_sim_status")
        if (stop_uscr_sim() || grepl("stopped by user", e$message, ignore.case = TRUE)) {
          update_model_debug(
            uscr_sim_debug,
            status = "stopped",
            stage = "Stopped by user",
            finished_at = Sys.time(),
            guidance = "The simulated USCR run was stopped manually before completion.",
            raw_error = e$message,
            log_entry = "Simulated USCR run stopped by user."
          )
          showNotification("USCR (sim) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            uscr_sim_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("USCR", "simulated data", e$message),
            raw_error = e$message,
            log_entry = paste("Simulated USCR failed:", e$message)
          )
          showNotification(
            paste("USCR (sim) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        # If nimble compilation failed, errors can be inspected via nimble::printErrors() in the console
        return(NULL)
      },
      finally = {
        removeNotification("uscr_sim_status")
        uscr_sim_running(FALSE)
        stop_uscr_sim(FALSE)  # Reset stop flag
        stop_flags_env$stop_uscr_sim <- FALSE
      }
    )
    uscr_sim_fit(fit)
    if (!is.null(fit) && !stop_uscr_sim()) {
      removeNotification("uscr_sim_status")
      update_model_debug(
        uscr_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "Simulated USCR completed successfully. Review the summary above and the tuning history below.",
        log_entry = "Simulated USCR run completed."
      )
      showNotification("USCR (sim) complete!", type = "message")
    }
  })
  
  # --- REM: simulated ---
  
  observeEvent(input$run_rem_sim, {
    shared_truth <- sim()
    shared_d <- if (!is.null(shared_truth)) sim_data() else NULL

    if (is.null(shared_d)) {
      showNotification(
        "Run the shared spatial simulator in the Simulate data tab first.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    
    rem_sim_dataset_id(current_shared_sim_id())
    rem_sim_debug(make_model_debug(
      "REM",
      "Shared spatial simulator",
      "Shared spatial simulator",
      "REM uses total animal events per camera and camera-days derived from the shared spatial simulator."
    ))
    update_model_debug(
      rem_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "REM is running on the shared spatial simulator so it can be compared and combined with simulated USCR and TTE.",
      raw_error = NULL,
      context = summarize_shared_sim_context(shared_d, "REM", sim_truth = shared_truth$truth),
      log_entry = "Simulated REM shared-spatial run requested."
    )
    
    rem_sim_running(TRUE)
    rem_sim_fit(NULL)
    
    showNotification(
      "Running REM on shared simulated spatial data...",
      type = "message",
      duration = NULL,
      id = "rem_sim_status"
    )
    
    fit <- tryCatch(
      {
        status_callback <- make_encounter_status_callback(
          rem_sim_debug,
          "REM",
          notification_id = "rem_sim_status"
        )
        status_callback(
          stage = "setup",
          detail = "Preparing REM inputs from the shared spatial simulator and initial model build...",
          value = 0.1
        )
        
        run_REM(
          y            = shared_d$camera_counts,
          r_km         = shared_d$out$`Detection Distance` / 1000,
          camera_days  = shared_d$camera_days,
          theta_deg    = input$theta,
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = app_n_chains(),
          D_max        = input$D_max,
          log_v_mean   = input$log_v_mean,
          log_v_sd     = input$log_v_sd,
          sd_eps_shape = input$sd_eps_shape,
          sd_eps_rate  = input$sd_eps_rate,
          status_callback = status_callback
        )
      },
      error = function(e) {
        removeNotification("rem_sim_status")
        update_model_debug(
          rem_sim_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("REM", "simulated shared spatial data", e$message),
          raw_error = e$message,
          log_entry = paste("Simulated REM run failed:", e$message)
        )
        showNotification(
          paste("REM (simulated) failed:", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
      },
      finally = {
        removeNotification("rem_sim_status")
        rem_sim_running(FALSE)
      }
    )
    rem_sim_fit(fit)
    if (!is.null(fit)) {
      removeNotification("rem_sim_status")
      update_model_debug(
        rem_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "REM shared-simulator run completed successfully. Review the summary above and use Compare & combine after the other simulated models finish.",
        log_entry = "Simulated REM run completed."
      )
      showNotification("REM (simulated) complete!", type = "message")
    }
  })
  
  # --- TTE: simulated ---
  
  observeEvent(input$run_tte_sim, {
    shared_truth <- sim()
    shared_d <- if (!is.null(shared_truth)) sim_data() else NULL

    if (is.null(shared_d)) {
      showNotification(
        "Run the shared spatial simulator in the Simulate data tab first.",
        type = "error",
        duration = 6
      )
      return(NULL)
    }
    
    tte_sim_dataset_id(current_shared_sim_id())
    tte_sim_debug(make_model_debug(
      "TTE",
      "Shared spatial simulator",
      "Shared spatial simulator",
      "TTE uses total animal events per camera and camera-days derived from the shared spatial simulator."
    ))
    update_model_debug(
      tte_sim_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "TTE is running on the shared spatial simulator so it can be compared and combined with simulated USCR and REM.",
      raw_error = NULL,
      context = summarize_shared_sim_context(shared_d, "TTE", sim_truth = shared_truth$truth),
      log_entry = "Simulated TTE shared-spatial run requested."
    )
    
    tte_sim_running(TRUE)
    tte_sim_fit(NULL)
    
    showNotification(
      "Running TTE on shared simulated spatial data...",
      type = "message",
      duration = NULL,
      id = "tte_sim_status"
    )
    
    fit <- tryCatch(
      {
        status_callback <- make_encounter_status_callback(
          tte_sim_debug,
          "TTE",
          notification_id = "tte_sim_status"
        )
        status_callback(
          stage = "setup",
          detail = "Preparing TTE inputs from the shared spatial simulator and initial model build...",
          value = 0.1
        )
        
        run_TTE(
          y            = shared_d$camera_counts,
          r_km         = shared_d$out$`Detection Distance` / 1000,
          camera_days  = shared_d$camera_days,
          theta_deg    = input$theta,
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = app_n_chains(),
          D_max        = input$D_max,
          log_v_mean   = input$log_v_mean,
          log_v_sd     = input$log_v_sd,
          sd_eps_shape = input$sd_eps_shape,
          sd_eps_rate  = input$sd_eps_rate,
          status_callback = status_callback
        )
      },
      error = function(e) {
        removeNotification("tte_sim_status")
        update_model_debug(
          tte_sim_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("TTE", "simulated shared spatial data", e$message),
          raw_error = e$message,
          log_entry = paste("Simulated TTE run failed:", e$message)
        )
        showNotification(
          paste("TTE (simulated) failed:", e$message),
          type = "error",
          duration = NULL
        )
        return(NULL)
      },
      finally = {
        removeNotification("tte_sim_status")
        tte_sim_running(FALSE)
      }
    )
    tte_sim_fit(fit)
    if (!is.null(fit)) {
      removeNotification("tte_sim_status")
      update_model_debug(
        tte_sim_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "TTE shared-simulator run completed successfully. Review the summary above and use Compare & combine after the other simulated models finish.",
        log_entry = "Simulated TTE run completed."
      )
      showNotification("TTE (simulated) complete!", type = "message")
    }
  })
  
  # --- USCR: NPS ---
  
  start_uscr_nps_run <- function(resume_override = NULL, resumed_from_pause = FALSE) {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    uscr_fast_mode <- isTRUE(input$uscr_fast_mode_nps)
    uscr_nps_debug(make_model_debug(
      "USCR",
      "Uploaded field data",
      "Uploaded field data; simulated grid",
      "Total animal events per camera, camera-days, and buffered state space."
    ))
    update_model_debug(
      uscr_nps_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = if (resumed_from_pause) {
        "USCR resumed from the latest paused round settings. This panel will show whether the run is still in setup, adaptive tuning, or the final run."
      } else {
        "USCR supports uploaded field data. This panel will show whether the run is still in setup, adaptive tuning, or the final run."
      },
      raw_error = NULL,
      pause_queued = FALSE,
      stop_queued = FALSE,
      context = summarize_uscr_context(d, source_label = "Uploaded", fast_mode = uscr_fast_mode),
      log_entry = if (resumed_from_pause) {
        "Uploaded-data USCR resume requested from paused settings."
      } else if (uscr_fast_mode) {
        "Uploaded-data USCR run requested with fast test mode enabled."
      } else {
        "Uploaded-data USCR run requested."
      }
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_uscr_nps(FALSE)
    pause_uscr_nps(FALSE)
    stop_flags_env$stop_uscr_nps <- FALSE
    stop_flags_env$pause_uscr_nps <- FALSE
    uscr_nps_running(TRUE)
    uscr_nps_fit(NULL)  # Clear previous results
    
    showNotification(
      if (resumed_from_pause) {
        "Resuming USCR on uploaded data from paused settings."
      } else if (uscr_fast_mode) {
        "Running USCR on uploaded data (fast test mode)."
      } else {
        "Running USCR on uploaded data."
      },
      type = "message",
      duration = NULL,
      id = "uscr_nps_status"
    )
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_uscr_nps()) {
          showNotification("USCR (uploaded data) run cancelled.", type = "warning")
          return(NULL)
        }
        
        update_model_debug(
          uscr_nps_debug,
          stage = "Preparing run",
          log_entry = "Input checks passed. Starting USCR setup."
        )
        
        # Check stop flag again before running
        if (stop_uscr_nps()) {
          showNotification("USCR (uploaded data) run cancelled.", type = "warning")
          return(NULL)
        }

        interrupt_callback <- function() {
          check_interrupt_action(
            stop_flag_name = "stop_uscr_nps",
            pause_flag_name = "pause_uscr_nps"
          )
        }
        
        status_callback <- function(stage, detail = NULL, value = NULL) {
          check_stop_flag("stop_uscr_nps")
          stage_label <- switch(
            stage,
            setup = "Preparing state space and model code",
            tuning = "Adaptive tuning",
            final_run = "Final MCMC run / convergence check",
            stage
          )
          progress_detail <- detail %||% paste("USCR stage:", stage_label)
          update_model_debug(
            uscr_nps_debug,
            status = "running",
            stage = stage_label,
            log_entry = progress_detail
          )
          showNotification(
            progress_detail,
            type = "message",
            duration = NULL,
            id = "uscr_nps_status"
          )
        }
        
        do.call(
          run_USCR_app,
          c(
            list(
              out = d$out,
              camera_counts = d$camera_counts,
              camera_days = d$camera_days
            ),
            uscr_run_args(
              waic = TRUE,
              status_callback = status_callback,
              fast_mode = uscr_fast_mode,
              override = resume_override,
              interrupt_callback = interrupt_callback
            )
          )
        )
      },
      error = function(e) {
        removeNotification("uscr_nps_status")
        if (stop_uscr_nps()) {
          update_model_debug(
            uscr_nps_debug,
        status = "stopped",
        stage = "Stopped by user",
        finished_at = Sys.time(),
            guidance = "The uploaded-data USCR run was stopped manually before completion.",
            raw_error = e$message,
            pause_queued = FALSE,
            stop_queued = FALSE,
            log_entry = "Uploaded-data USCR run stopped by user."
          )
          showNotification("USCR (uploaded data) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            uscr_nps_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("USCR", "uploaded field data", e$message),
            raw_error = e$message,
            pause_queued = FALSE,
            stop_queued = FALSE,
            log_entry = paste("Uploaded-data USCR failed:", e$message)
          )
          showNotification(
            paste("USCR (uploaded data) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        removeNotification("uscr_nps_status")
        uscr_nps_running(FALSE)
        stop_uscr_nps(FALSE)  # Reset stop flag
        pause_uscr_nps(FALSE)
        stop_flags_env$stop_uscr_nps <- FALSE
        stop_flags_env$pause_uscr_nps <- FALSE
      }
    )
    uscr_nps_fit(fit)
    if (!is.null(fit) && isTRUE(fit$provisional)) {
      uscr_nps_resume_settings(fit$resumable_settings)
      removeNotification("uscr_nps_status")
      update_model_debug(
        uscr_nps_debug,
        status = "paused",
        stage = "Paused after completed round",
        finished_at = Sys.time(),
        guidance = paste(
          "A provisional USCR estimate was saved from the most recent completed round.",
          "Use it for troubleshooting only, then click 'Pause/Resume' to continue from these settings."
        ),
        pause_queued = FALSE,
        stop_queued = FALSE,
        log_entry = "Uploaded-data USCR paused and provisional fit saved."
      )
      showNotification(
        "USCR paused after the latest completed round. A provisional estimate was saved.",
        type = "warning",
        duration = 8
      )
    } else if (!is.null(fit) && !stop_uscr_nps()) {
      uscr_nps_resume_settings(NULL)
      removeNotification("uscr_nps_status")
      update_model_debug(
        uscr_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "Uploaded-data USCR completed successfully. Review the summary above and the tuning history below.",
        pause_queued = FALSE,
        stop_queued = FALSE,
        log_entry = "Uploaded-data USCR run completed."
      )
      showNotification("USCR (uploaded data) complete!", type = "message")
    }
  }

  observeEvent(input$run_uscr_nps, {
    if (isTRUE(uscr_nps_running())) {
      showNotification("USCR (uploaded data) is already running.", type = "warning", duration = 4)
      return(NULL)
    }
    start_uscr_nps_run()
  })
  
  # --- REM: NPS ---
  
  observeEvent(input$run_rem_nps, {
    if (rem_nps_running()) {
      showNotification(
        "REM (uploaded data) is already running in this session. Wait for it to finish before starting another REM job.",
        type = "warning",
        duration = 6
      )
      return(NULL)
    }
    req(nps_model_inputs())
    d <- nps_model_inputs()
    requested_chains <- app_n_chains()
    rem_context <- c(
      summarize_rem_context(d),
      paste("Chains requested in UI:", requested_chains),
      "Background worker run: yes; the app enforces a minimum of 2 chains."
    )
    rem_args <- list(
      y            = d$camera_counts,
      r_km         = d$out$`Detection Distance` / 1000,
      camera_days  = d$camera_days,
      theta_deg    = get_camera_angle_vector(d$out, fallback_deg = input$theta),
      iter         = input$iter_rem_tte,
      burnin       = input$burnin_rem_tte,
      thin         = input$thin_rem_tte,
      n_chains     = requested_chains,
      D_max        = input$D_max,
      log_v_mean   = input$log_v_mean,
      log_v_sd     = input$log_v_sd,
      sd_eps_shape = input$sd_eps_shape,
      sd_eps_rate  = input$sd_eps_rate
    )
    rem_nps_debug(make_model_debug(
      "REM",
      "Uploaded field data",
      "Uploaded field data only",
      "REM uses total animal events per camera and camera-days."
    ))
    update_model_debug(
      rem_nps_debug,
      status = "running",
      stage = "Queued background run",
      started_at = Sys.time(),
      guidance = paste(
        "REM on uploaded field data is running in a background worker so this session should stay responsive.",
        "This path does not provide live round-by-round tuning updates while the worker is running,",
        "but the requested run settings are shown below and the full tuning history will appear after completion."
      ),
      raw_error = NULL,
      context = rem_context,
      log_entry = paste(
        "Uploaded-data REM run requested.",
        paste0(
          "Requested settings: iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", requested_chains,
          "."
        )
      )
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_rem_nps(FALSE)
    stop_flags_env$stop_rem_nps <- FALSE
    rem_nps_running(TRUE)
    rem_nps_background(TRUE)
    rem_nps_fit(NULL)  # Clear previous results
    
    update_model_debug(
      rem_nps_debug,
      stage = "Background run in progress",
      log_entry = paste(
        "Submitted REM background job.",
        "Chains used =", requested_chains,
        "(minimum enforced by the app).",
        "Live tuning updates are not available for this background path.",
        "REM may internally increase iterations until max Rhat is acceptable."
      )
    )
    
    showNotification(
      paste(
        "REM background run started.",
        paste0(
          "Requested settings: iter=", input$iter_rem_tte,
          ", burn-in=", input$burnin_rem_tte,
          ", thin=", input$thin_rem_tte,
          ", chains=", requested_chains,
          "."
        ),
        "This session should stay responsive, but live tuning updates are not available for this path."
      ),
      type = "message",
      duration = NULL,
      id = "rem_nps_progress"
    )
    
    rem_future <- future::future({
      do.call(run_REM, rem_args)
    })
    
    promises::as.promise(rem_future) %...>% (function(fit) {
      removeNotification("rem_nps_progress")
      rem_nps_fit(fit)
      rem_nps_running(FALSE)
      rem_nps_background(FALSE)
      stop_rem_nps(FALSE)
      stop_flags_env$stop_rem_nps <- FALSE
      update_model_debug(
        rem_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "REM completed successfully in a background worker. Review the summary above and the tuning history below.",
        log_entry = "Uploaded-data REM background run completed."
      )
      showNotification("REM (uploaded data) complete!", type = "message", duration = 5)
      invisible(NULL)
    }) %...!% (function(e) {
      removeNotification("rem_nps_progress")
      rem_nps_fit(NULL)
      rem_nps_running(FALSE)
      rem_nps_background(FALSE)
      stop_rem_nps(FALSE)
      stop_flags_env$stop_rem_nps <- FALSE
      update_model_debug(
        rem_nps_debug,
        status = "error",
        stage = "Failed",
        finished_at = Sys.time(),
        guidance = friendly_model_error("REM", "uploaded field data", conditionMessage(e)),
        raw_error = conditionMessage(e),
        log_entry = paste("Uploaded-data REM background run failed:", conditionMessage(e))
      )
      showNotification(
        paste("REM (uploaded data) failed:", conditionMessage(e)),
        type = "error",
        duration = NULL
      )
      invisible(NULL)
    })
  })
  
  # --- TTE: NPS ---
  
  observeEvent(input$run_tte_nps, {
    req(nps_model_inputs())
    d <- nps_model_inputs()
    tte_nps_debug(make_model_debug(
      "TTE",
      "Uploaded field data",
      "Uploaded field data only",
      "Current build passes total animal events per camera and camera-days."
    ))
    update_model_debug(
      tte_nps_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "TTE only runs on uploaded field data in this app.",
      raw_error = NULL,
      context = summarize_tte_context(d),
      log_entry = "Uploaded-data TTE run requested."
    )
    validate(
      need(all(d$camera_days > 0),
           "Some cameras have zero camera-days; check deployment dates.")
    )
    
    # Reset stop flag
    stop_tte_nps(FALSE)
    tte_nps_running(TRUE)
    tte_nps_fit(NULL)  # Clear previous results
    
    showNotification(
      "Running TTE on uploaded data... This may take several minutes.",
      type = "message",
      duration = NULL,
      id = "tte_nps_status"
    )
    
    fit <- tryCatch(
      {
        # Check if stopped before starting
        if (stop_tte_nps()) {
          showNotification("TTE (uploaded data) run cancelled.", type = "warning")
          return(NULL)
        }
        
        status_callback <- make_encounter_status_callback(
          tte_nps_debug,
          "TTE",
          notification_id = "tte_nps_status"
        )
        status_callback(
          stage = "setup",
          detail = "Preparing TTE inputs from uploaded data and initial model build...",
          value = 0.1
        )
        
        # Check stop flag again before running
        if (stop_tte_nps()) {
          showNotification("TTE (uploaded data) run cancelled.", type = "warning")
          return(NULL)
        }
        
        run_TTE(
          y            = d$camera_counts,
          r_km         = d$out$`Detection Distance` / 1000,
          camera_days  = d$camera_days,
          theta_deg    = get_camera_angle_vector(d$out, fallback_deg = input$theta),
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = app_n_chains(),
          D_max        = input$D_max,
          log_v_mean   = input$log_v_mean,
          log_v_sd     = input$log_v_sd,
          sd_eps_shape = input$sd_eps_shape,
          sd_eps_rate  = input$sd_eps_rate,
          status_callback = status_callback
        )
      },
      error = function(e) {
        removeNotification("tte_nps_status")
        if (stop_tte_nps()) {
          update_model_debug(
            tte_nps_debug,
            status = "stopped",
            stage = "Stopped by user",
            finished_at = Sys.time(),
              guidance = "The uploaded-data TTE run was stopped manually before completion.",
              raw_error = e$message,
              log_entry = "Uploaded-data TTE run stopped by user."
          )
          showNotification("TTE (uploaded data) run stopped by user.", type = "warning")
        } else {
          update_model_debug(
            tte_nps_debug,
            status = "error",
            stage = "Failed",
            finished_at = Sys.time(),
            guidance = friendly_model_error("TTE", "uploaded field data", e$message),
            raw_error = e$message,
            log_entry = paste("Uploaded-data TTE failed:", e$message)
          )
          showNotification(
            paste("TTE (uploaded data) failed:", e$message),
            type = "error", duration = NULL
          )
        }
        return(NULL)
      },
      finally = {
        removeNotification("tte_nps_status")
        tte_nps_running(FALSE)
        stop_tte_nps(FALSE)  # Reset stop flag
      }
    )
    tte_nps_fit(fit)
    if (!is.null(fit) && !stop_tte_nps()) {
      removeNotification("tte_nps_status")
      update_model_debug(
        tte_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = "TTE completed successfully. Review the summary above.",
        log_entry = "Uploaded-data TTE run completed."
      )
      showNotification("TTE (uploaded data) complete!", type = "message")
    }
  })
  
  # ================================================================
  # MODEL TAB OUTPUTS (USCR / REM / TTE)
  # ================================================================
  
  # USCR summaries
  output$uscr_sim_text <- renderPrint({
    dbg <- uscr_sim_debug()
    if (uscr_sim_running()) {
      cat("⏳ USCR model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("USCR (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("USCR (simulated): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- uscr_sim_fit()
    if (is.null(fit)) {
      cat("USCR (simulated): not run yet. Click 'Run USCR on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      note                      = if (!is.null(sim()))
        sprintf("True simulated density = %.1f animals/km²", sim()$truth$D_per_km2)
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    out
  })
  
  output$uscr_nps_text <- renderPrint({
    dbg <- uscr_nps_debug()
    if (uscr_nps_running()) {
      cat("⏳ USCR model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("USCR (uploaded data): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("USCR (uploaded data): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- uscr_nps_fit()
    if (is.null(fit)) {
      cat("USCR (uploaded data): not run yet. Click 'Run USCR on uploaded data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Uploaded field data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    if (isTRUE(fit$provisional)) {
      out$Status <- "Provisional paused estimate"
      out$Warning <- fit$provisional_reason %||% "This estimate comes from the most recent completed USCR round and should not be interpreted as a final converged result."
    }
    out
  })
  
  # REM summaries
  output$rem_sim_text <- renderPrint({
    dbg <- rem_sim_debug()
    if (rem_sim_running()) {
      cat("⏳ REM model is running on the shared spatial simulator...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("REM (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    fit <- rem_sim_fit()
    if (is.null(fit)) {
      cat("REM (simulated): not run yet. Generate the shared simulated dataset first, then click 'Run REM on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    truth <- sim()
    list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(truth))
        sprintf(
          "%s true simulated density = %.1f animals/km²",
          "Shared spatial simulator",
          truth$truth$D_per_km2
        )
    )
  })
  
  output$rem_nps_text <- renderPrint({
    dbg <- rem_nps_debug()
    if (rem_nps_running()) {
      cat("⏳ REM model is running in the background...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("This session should stay responsive while the REM job runs.\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("REM (uploaded data): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("REM (uploaded data): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- rem_nps_fit()
    if (is.null(fit)) {
      cat("REM (uploaded data): not run yet. Click 'Run REM on uploaded data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Uploaded field data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    out
  })
  
  # TTE summaries
  output$tte_sim_text <- renderPrint({
    dbg <- tte_sim_debug()
    if (tte_sim_running()) {
      cat("⏳ TTE model is running on the shared spatial simulator...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("TTE (simulated): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    fit <- tte_sim_fit()
    if (is.null(fit)) {
      cat("TTE (simulated): not run yet. Generate the shared simulated dataset first, then click 'Run TTE on simulated data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    truth <- sim()
    list(
      dataset                   = "Shared spatial simulator",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2)),
      WAIC                      = round(s$waic, 2),
      note                      = if (!is.null(truth))
        sprintf(
          "%s true simulated density = %.1f animals/km²",
          "Shared spatial simulator",
          truth$truth$D_per_km2
        )
    )
  })
  
  output$tte_nps_text <- renderPrint({
    dbg <- tte_nps_debug()
    if (tte_nps_running()) {
      cat("⏳ TTE model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("TTE (uploaded data): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "stopped")) {
      cat("TTE (uploaded data): the last run was stopped before completion.\n")
      return(invisible(NULL))
    }
    fit <- tte_nps_fit()
    if (is.null(fit)) {
      cat("TTE (uploaded data): not run yet. Click 'Run TTE on uploaded data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Uploaded field data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    out
  })
  
  output$uscr_sim_debug <- renderText({
    state <- uscr_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, uscr_sim_fit())
  })
  
  output$uscr_nps_debug <- renderText({
    state <- uscr_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, uscr_nps_fit())
  })
  
  output$rem_sim_debug <- renderText({
    state <- rem_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, rem_sim_fit())
  })
  
  output$rem_nps_debug <- renderText({
    state <- rem_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, rem_nps_fit())
  })
  
  output$tte_sim_debug <- renderText({
    state <- tte_sim_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, tte_sim_fit())
  })
  
  output$tte_nps_debug <- renderText({
    state <- tte_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, tte_nps_fit())
  })
  
  # ================================================================
  # COMPARE & COMBINE (WAIC-based)
  # ================================================================
  
  shared_sim_fits <- reactive({
    current_id <- current_shared_sim_id()
    if (is.null(current_id)) return(list())

    fits <- list(
      REM = if (identical(rem_sim_dataset_id(), current_id)) rem_sim_fit() else NULL,
      TTE = if (identical(tte_sim_dataset_id(), current_id)) tte_sim_fit() else NULL,
      USCR = if (identical(uscr_sim_dataset_id(), current_id)) uscr_sim_fit() else NULL
    )
    fits[!vapply(fits, is.null, logical(1))]
  })

  sim_combo <- reactive({
    fits <- shared_sim_fits()
    if (!length(fits)) return(NULL)
    build_combo_table_from_fits(
      fits,
      density_threshold = input$combo_density_threshold
    )
  })

  sim_combo_interval_plot_obj <- reactive({
    fits <- shared_sim_fits()
    if (!length(fits)) return(NULL)
    build_combo_interval_plot(
      fits,
      title_text = "Shared spatial simulator: density estimates by model"
    )
  })

  sim_combo_density_plot_obj <- reactive({
    fits <- shared_sim_fits()
    if (!length(fits)) return(NULL)
    build_combo_density_overlay_plot(
      fits,
      title_text = "Shared spatial simulator: posterior density overlap"
    )
  })
  
  nps_combo <- reactive({
    build_combo_table_from_fits(
      list(
        REM = rem_nps_fit(),
        TTE = tte_nps_fit(),
        USCR = uscr_nps_fit()
      ),
      density_threshold = input$combo_density_threshold
    )
  })

  nps_combo_fits <- reactive({
    fits <- list(
      REM = rem_nps_fit(),
      TTE = tte_nps_fit(),
      USCR = uscr_nps_fit()
    )
    fits[!vapply(fits, is.null, logical(1))]
    fits <- fits[!vapply(fits, function(x) isTRUE(x$provisional), logical(1))]
    fits
  })

  nps_combo_interval_plot_obj <- reactive({
    fits <- nps_combo_fits()
    if (!length(fits)) return(NULL)
    build_combo_interval_plot(
      fits,
      title_text = "Uploaded field data: density estimates by model"
    )
  })

  nps_combo_density_plot_obj <- reactive({
    fits <- nps_combo_fits()
    if (!length(fits)) return(NULL)
    build_combo_density_overlay_plot(
      fits,
      title_text = "Uploaded field data: posterior density overlap"
    )
  })
  
  output$sim_combo_table <- renderDT({
    combo <- sim_combo()
    if (is.null(combo)) {
      return(DT::datatable(
        data.frame(Note = "Run one or more models from the shared spatial simulator first."),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      ))
    }
    df <- combo$table %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
    datatable(df, options = list(pageLength = 5, dom = "t", autoWidth = TRUE), rownames = FALSE)
  })

  output$sim_combo_interval_plot <- renderPlot({
    p <- sim_combo_interval_plot_obj()
    validate(need(!is.null(p), "Run one or more shared-simulation models first."))
    p
  }, res = 120)

  output$sim_combo_density_plot <- renderPlot({
    p <- sim_combo_density_plot_obj()
    validate(need(!is.null(p), "Run one or more shared-simulation models first."))
    p
  }, res = 120)
  
  output$dl_sim_uscr_csv <- downloadHandler(
    filename = function() {
      paste0("DEER_shared_sim_posterior_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fits <- shared_sim_fits()
      req(length(fits) > 0)
      out <- purrr::imap_dfr(fits, function(fit, model_name) {
        df <- posterior_summary_df(fit)
        if (is.null(df)) return(NULL)
        dplyr::mutate(df, model = model_name, .before = 1)
      })
      req(nrow(out) > 0)
      readr::write_csv(out, file)
    }
  )

  output$dl_sim_combo_interval_png <- downloadHandler(
    filename = function() {
      paste0("DEER_shared_sim_interval_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- sim_combo_interval_plot_obj()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )

  output$dl_sim_combo_density_png <- downloadHandler(
    filename = function() {
      paste0("DEER_shared_sim_density_overlay_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- sim_combo_density_plot_obj()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )
  
  output$dl_nps_all_csv <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_posterior_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fits <- nps_combo_fits()
      req(length(fits) > 0)
      df <- dplyr::bind_rows(
        lapply(names(fits), function(model_name) {
          dplyr::mutate(
            posterior_summary_df(fits[[model_name]]),
            model = model_name,
            .before = 1
          )
        })
      )
      readr::write_csv(df, file)
    }
  )
  
  output$nps_combo_table <- renderDT({
    combo <- nps_combo()
    if (is.null(combo)) {
      note_text <- if (isTRUE(uscr_nps_fit()$provisional)) {
        "A paused provisional USCR fit is currently saved, but provisional USCR fits are excluded from Compare & combine until you resume or rerun that model."
      } else {
        "Run at least one uploaded-data model from its tab first. The table will update as REM, TTE, and USCR finish."
      }
      return(DT::datatable(
        data.frame(
          Note = note_text
        ),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      ))
    }
    df <- combo$table %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
    datatable(df, options = list(pageLength = 5, dom = "t", autoWidth = TRUE), rownames = FALSE)
  })

  output$nps_combo_interval_plot <- renderPlot({
    p <- nps_combo_interval_plot_obj()
    validate(need(!is.null(p), "Run at least one uploaded-data model first."))
    p
  }, res = 120)

  output$nps_combo_density_plot <- renderPlot({
    p <- nps_combo_density_plot_obj()
    validate(need(!is.null(p), "Run at least one uploaded-data model first."))
    p
  }, res = 120)

  output$dl_nps_combo_interval_png <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_interval_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- nps_combo_interval_plot_obj()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )

  output$dl_nps_combo_density_png <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_density_overlay_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- nps_combo_density_plot_obj()
      req(!is.null(p))
      ggplot2::ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------

shinyApp(ui, server)
