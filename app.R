# Deer Density Lab – USCR vs REM vs TTE
# Models run via run_USCR(), run_REM(), run_TTE() in R/sim_and_models.R

# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------

needs <- c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr", "stringr", "data.table",
  "leaflet", "ggrepel",
  "nimble", "parallel", "MCMCvis", "lubridate", "sf", "tibble"
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
  library(data.table)
  library(leaflet)
  library(ggrepel)
  library(tibble)
  library(sf)
  library(nimble)
})

# -------------------------------------------------------------------
# Helper files
# -------------------------------------------------------------------

source("R/sim_and_models.R")   # run_USCR(), run_REM(), run_TTE()
source("R/data_checks.R")      # QC + summary helpers

# Pin the active USCR implementation to the current file contents so stale
# objects from prior sessions cannot route the app through an older wrapper.
uscr_model_env <- new.env(parent = globalenv())
sys.source("R/sim_and_models.R", envir = uscr_model_env)
run_USCR_app <- uscr_model_env$run_USCR

# -------------------------------------------------------------------
# Default priors (used for the Advanced inputs' starting values, and always
# used in Simple mode)
# -------------------------------------------------------------------
# REM/TTE upper limit on density under the uniform prior (animals/km^2). Not
# user-facing: the app doubles it if the posterior for D presses against it.
rem_tte_start_D_max <- 200

advanced_prior_defaults <- list(
  log_v_mean     = 1.130,
  log_v_sd       = 0.3372,
  sd_eps_shape   = 1,
  sd_eps_rate    = 1,
  log_sigma_mean = -1.5269,
  log_sigma_sd   = 0.1535,
  log_lam0_mean  = 0,
  log_lam0_sd    = 1
)

# -------------------------------------------------------------------
# Simulation helpers
# -------------------------------------------------------------------
# Simulation and model-input helpers now come directly from R/sim_and_models.R:
#   simulate_camera_counts(), sim_model_inputs(), build_nps_model_inputs()

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
  D_km2 <- fit$samples_all[, "D"]

  list(
    mean_km2  = mean(D_km2),
    median_km2 = stats::median(D_km2),
    sd_km2    = sd(D_km2),
    q2.5_km2  = quantile(D_km2, 0.025),
    q97.5_km2 = quantile(D_km2, 0.975),
    waic      = get_waic_value(fit)
  )
}

extract_density_draws_km2 <- function(fit) {
  if (is.null(fit) || is.null(fit$samples_all) || !"D" %in% colnames(fit$samples_all)) {
    return(numeric(0))
  }
  x <- as.numeric(fit$samples_all[, "D"])
  x[is.finite(x)]
}

generate_mixture_posterior <- function(draws_list, weights, n_draws = NULL) {

  if (is.null(n_draws)) {
    n_draws <- min(lengths(draws_list))
  }

  model_ids <- sample(
    seq_along(draws_list),
    size = n_draws,
    replace = TRUE,
    prob = weights
  )

  mixture_draws <- numeric(n_draws)

  for (i in seq_len(n_draws)) {

    m <- model_ids[i]

    mixture_draws[i] <- sample(
      draws_list[[m]],
      size = 1,
      replace = TRUE
    )
  }

  mixture_draws
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
  #For mixture distribution
  raw_density_draws <- lapply(
    fits,
    extract_density_draws_km2
  )
  draw_lengths <- lengths(density_draws)
  if (any(!is.finite(draw_lengths)) || any(draw_lengths < 1L)) return(NULL)

  # Quantile-based averaging (below) assumes every model contributes the same
  # number of posterior draws, so that sorted position s means the same thing
  # (the same quantile) in each model. In practice this always holds because
  # D is a bounded quantity for all three models (REM/TTE: D_max times a Beta
  # draw; USCR: N / area_km2, with N in [0, M] and area_km2 a fixed
  # constant), so it can't produce non-finite draws that get filtered out
  # unevenly across models. If that assumption is ever violated,
  # fail loudly rather than silently truncating every model down to the
  # shortest one, which would quietly drop the upper tail of the longer
  # model(s) instead of matching quantiles.
  if (length(unique(draw_lengths)) > 1L) {
    stop(
      "prepare_combo_density_data: models have unequal numbers of finite ",
      "posterior draws (",
      paste0(model_names, " = ", draw_lengths, collapse = ", "),
      "). Quantile-based model averaging assumes an equal number of draws ",
      "per model; this most likely means D contained non-finite values ",
      "in one model's samples that were filtered out. Investigate the ",
      "affected fit(s) instead of truncating to the shortest model.",
      call. = FALSE
    )
  }

  mat <- do.call(cbind, density_draws)
  colnames(mat) <- model_names

  density_est <- mat
  summary_rows <- character()
  if (ncol(mat) > 1L) {

    density_est <- cbind(
      mat,
      `unweighted mean` = rowMeans(mat)
    )

    summary_rows <- "unweighted mean"

    if (isTRUE(waic_weights_available)) {

      # weighted quantile averaging
      density_est <- cbind(
        density_est,
        `weighted quantile mean` =
          as.numeric(mat %*% waic_tbl$w)
      )

      # mixture distribution averaging
      mixture_draws <- generate_mixture_posterior(
        draws_list = raw_density_draws,
        weights = waic_tbl$w,
        n_draws = nrow(mat)
      )

      density_est <- cbind(
        density_est,
        `mixture distribution` =
          mixture_draws
      )

      summary_rows <- c(
        summary_rows,
        "weighted quantile mean",
        "mixture distribution"
      )
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
  medians <- apply(density_est, 2, stats::median)
  lower <- apply(density_est, 2, stats::quantile, probs = 0.025)
  upper <- apply(density_est, 2, stats::quantile, probs = 0.975)
  prob_threshold <- apply(density_est, 2, function(x) mean(x > density_threshold))
  # Posterior predictive check (Bayesian p-value) for each model; summary rows get NA.
  ppc_p <- vapply(combo_data$fits, function(f) {
    x <- ppc_summary(f)
    if (is.null(x)) NA_real_ else x$p
  }, numeric(1))
  threshold_label <- paste0(
    "P(D > ",
    format(density_threshold, trim = TRUE, scientific = FALSE),
    " animals/km²)"
  )

  table_out <- tibble::tibble(
    Method = c(model_names, summary_rows),
    `Mean density (animals/km²)` = as.numeric(means),
    `Median density (animals/km²)` = as.numeric(medians),
    `Lower 2.5%` = as.numeric(lower),
    `Upper 97.5%` = as.numeric(upper),
    `WAIC weight` = c(waic_tbl$w, rep(NA_real_, length(summary_rows))),
    `Bayesian p-value` = c(ppc_p, rep(NA_real_, length(summary_rows)))
  )
  table_out[[threshold_label]] <- as.numeric(prob_threshold)
  table_out <- table_out[, c(
    "Method",
    "Mean density (animals/km²)",
    "Median density (animals/km²)",
    "Lower 2.5%",
    "Upper 97.5%",
    threshold_label,
    "WAIC weight",
    "Bayesian p-value"
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
  "weighted quantile mean" = "#2F3E46",
  "mixture distribution" = "#8E44AD"
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
        Method == "weighted quantile mean" ~
          "Weighted quantile mean",

        Method == "mixture distribution" ~
          "Mixture distribution",
        TRUE ~ "Unweighted mean"
      ),
      Method = factor(Method, levels = rev(Method)),
      point_shape = dplyr::case_when(
        Summary == "Model" ~ 16,
        Summary == "Weighted quantile mean" ~ 15,
        Summary == "Mixture distribution" ~ 18,
        TRUE ~ 17
      )
    )

  interval_levels <- rev(colnames(combo_data$density_est))
  interval_colors <- combo_plot_colors[interval_levels]
  interval_colors <- interval_colors[!is.na(interval_colors)]

  ggplot(interval_df, aes(x = mean, y = Method, color = Method)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 1, alpha = 0.95) +
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

#' Posterior summary for CSV export (parameter names, mean, 95% CI)
posterior_summary_df <- function(fit) {
  if (is.null(fit) || is.null(fit$samples_all)) return(NULL)
  m <- as.data.frame(fit$samples_all)

  # Density: all three models monitor "D" in animals/km² (USCR derives it
  # as N / area_km2). Export it as "D_km2" so the units are explicit.
  if ("D" %in% names(m)) {
    m$D_km2 <- m$D
    m$D     <- NULL
  }

  # log_sigma/log_lam_0 (USCR) duplicate sigma/lam_0 on the log scale, and
  # sum_obs/sum_sim are just the intermediate discrepancy sums used to
  # compute bp (the Bayesian p-value) — drop all four rather than exporting
  # redundant or purely intermediate columns.
  drop_cols <- intersect(c("log_sigma", "log_lam_0", "sum_obs", "sum_sim"), names(m))
  if (length(drop_cols) > 0) m[drop_cols] <- NULL

  tibble::tibble(
    parameter = names(m),
    mean      = round(sapply(m, mean), 2),
    median    = round(sapply(m, stats::median), 2),
    q025      = round(sapply(m, function(x) stats::quantile(x, 0.025)), 2),
    q975      = round(sapply(m, function(x) stats::quantile(x, 0.975)), 2)
  )
}
  #' Shared legend describing posterior summary CSV columns, shown below the
  #' download buttons on the Compare & combine tab.
  posterior_param_legend <- function() {
    markdown(paste(
      "**Posterior summary column reference:**",
      "",
      "- `bp` — Bayesian p-value (posterior predictive check); values below 0.05 or above 0.95 mean the check failed.",
      "- `sd_eps` (σ<sub>ε</sub>) — standard deviation of the camera-level lognormal random effect in the Poisson-lognormal count model. Only included when camera-level noise is turned on (**Model settings > Advanced > Priors: Camera-level noise**); with it off, `sd_eps` is not estimated and does not appear here.",
      "- `v` — animal movement speed (km/day).",
      "- `D_km2` — density, in animals per km².",
      "- `N` — USCR only: abundance, the total number of animals estimated within the state space.",
      "- `lam_0` (λ₀) — USCR only: expected detections per deployed day if an animal's activity center were at the camera.",
      "- `psi` (ψ) — USCR only: the probability that a data-augmented individual is a real member of the population. Its prior sets the prior on abundance and density.",
      "- `sigma` (σ) — USCR only: spatial scale parameter (km), related to home-range size and how quickly detection falls with distance from an animal's activity center.",
      sep = "\n"
    ))
  }

# -------------------------------------------------------------------
# UI
# -------------------------------------------------------------------

ui <- page_fillable(
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
    # Shiny hides each file input's real <input> 99,999 px off the top of the
    # page. When a file dialog closes, Safari focuses that input and scrolls
    # to it, jumping to the top of the tab. Keeping the hidden input in place
    # (invisible) next to its Browse button avoids the jump. The custom
    # message clears a file input's display (used when an upload is reset).
    tags$script(HTML("
      function deerPinFileInputs() {
        document.querySelectorAll('input[type=\"file\"]').forEach(function (el) {
          ['top', 'left'].forEach(function (p) { el.style.setProperty(p, '0', 'important'); });
          el.style.setProperty('width', '1px', 'important');
          el.style.setProperty('height', '1px', 'important');
          el.style.setProperty('opacity', '0', 'important');
          el.style.setProperty('pointer-events', 'none', 'important');
        });
      }
      document.addEventListener('DOMContentLoaded', deerPinFileInputs);
      $(document).on('shiny:connected', function () {
        deerPinFileInputs();
        Shiny.addCustomMessageHandler('reset_file_input', function (id) {
          var el = document.getElementById(id);
          if (!el) return;
          el.value = '';
          var group = el.closest('.input-group');
          if (group) {
            var text = group.querySelector('input[type=\"text\"]');
            if (text) text.value = '';
          }
          var bar = document.getElementById(id + '_progress');
          if (bar) bar.style.visibility = 'hidden';
        });
      });
    ")),
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
      .prior-details {
        margin-bottom: 0.9rem;
      }
      .prior-details:last-of-type {
        margin-bottom: 3rem;
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
      .about-card h3 {
        font-size: 1.2rem;
        font-weight: 600;
        margin-top: 1rem;
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
      a {
        color: var(--rw1);
        text-decoration: none;
        border-bottom: 1px dotted rgba(0,0,0,.18);
      }
      a.text-link {
        color: var(--rw1);
        text-decoration: underline;
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
          style = "max-width: 980px; margin: 0 auto; padding: 0 16px; text-align: center;",
          tags$h1(style = "margin-top: 0.5rem;", "DEER App"),
          tags$p(
            class = "lead",
            style = "margin-top: 0.5rem;",
            "Density Estimation from Encounter Rates"
          ),
          tags$p(
            class = "small",
            "USCR · REM · TTE — three unmarked camera methods for estimating animal density, reported as animals/km²."
          )
        )
      ),
      tags$div(
        style = "max-width: 980px; margin: 0 auto; padding: 0 16px 42px;",

        # Welcome
        tags$section(
          tags$h2("Welcome to the DEER App"),
          tags$h3("What it does"),
          tags$p(
            class = "lead",
            "The DEER App estimates ",
            tags$strong("animal density at sites"),
            " from ",
            tags$strong("unmarked camera detections", .noWS = "after"),
            ". It fits ",
            tags$strong("three complementary models", .noWS = "after"),
            ", reports ",
            tags$strong("credible intervals", .noWS = "after"),
            ", and lets users ",
            tags$strong("compare, diagnose, and average results", .noWS = "after"),
            ". The upload workflow expects data on ",
            tags$strong("camera deployments"),
            " (e.g., where and when cameras were set and began recording) and ",
            tags$strong("images"),
            " (e.g., the date and time animals were detected at each camera). The broader goal is a ",
            tags$strong("consistent, reproducible workflow"),
            " for estimating ",
            tags$strong("animal density"),
            " from camera-trap data without requiring individually identifiable animals."
          )
        ),

        tags$div(class = "divider"),

        # How to use
        tags$section(
          tags$h2("How to use the DEER App"),
          tags$div(
            class = "about-card",

            tags$h3("Step 1: Upload field data"),
            tags$p(
              "The main analysis workflow uses uploaded field data. The app fits ",
              tags$strong("Unmarked Spatial Capture-Recapture (USCR)", .noWS = "after"),
              ", ",
              tags$strong("Random Encounter Model (REM)", .noWS = "after"),
              ", and ",
              tags$strong("Time-to-Event (TTE)"),
              " models to camera-trap data."
            ),
            tags$p("You will need two data files:"),
            tags$ul(
              tags$li(
                tags$strong("Deployment CSV"),
                " — camera deployment information such as locations, dates, and detection distances."
              ),
              tags$li(
                tags$strong("Images CSV"),
                " — detection records such as timestamps, species, and ",
                tags$strong("Cluster ID"),
                " values. ",
                tags$strong("Cluster ID"),
                " identifies which images belong to the same independent encounter event; multiple images (or species) from the same event share the same Cluster ID."
              )
            ),
            tags$p(
              "The upload pipeline expects the ",
              tags$strong("exact column names"),
              " described in the ",
              tags$strong("Add your data"),
              " tab, which lists column requirements and examples. The app will automatically:"
            ),
            tags$ul(
              tags$li("Check required columns and data types."),
              tags$li("Optionally trim each camera deployment to meet study design criteria.")
            ),

            tags$h3("Step 2: Adjust settings (optional)"),
            tags$p(
              "The default settings work well for many studies. To change them, open the ",
              tags$strong("Model settings"),
              " tab and choose ",
              tags$strong("Advanced", .noWS = "after"),
              ". There you can:"
            ),
            tags$ul(
              tags$li("Give all three models an informative prior on density, as a 95% interval."),
              tags$li("Add camera-level noise (a Poisson-lognormal model) to account for overdispersed counts."),
              tags$li("Adjust the prior distributions for every other parameter.")
            ),
            tags$p(
              "The default camera detection angle is set at the top of the ",
              tags$strong("Model settings"),
              " tab, or you can give camera-specific values in a ",
              tags$code("Camera Detection Angle"),
              " column of the deployment file."
            ),

            tags$h3("Step 3: Run the models"),
            tags$p(
              "Fit each model from its own tab: ",
              tags$strong("USCR model", .noWS = "after"),
              ", ",
              tags$strong("REM model", .noWS = "after"),
              ", and ",
              tags$strong("TTE model", .noWS = "after"),
              ". Each tab has a Run button, a note listing any settings that differ from the defaults, and a run-status and troubleshooting panel."
            ),
            tags$ul(
              tags$li("Click the Run button. Results appear below it when the model finishes."),
              tags$li(
                "Runs are adaptive: after each round, the app checks convergence and whether its settings are large enough, adjusts them, and reruns until every check passes or the adaptive-round limit is reached (6 rounds by default; change it in the ",
                tags$strong("Model settings"),
                " tab)."
              ),
              tags$li("Each run also reports a posterior predictive check, and flags a failed check.")
            ),
            tags$p(
              tags$em("Note:"),
              " model runs can take a long time, especially USCR at high-density sites or with large camera arrays. From the second adaptive round on, a pop-up estimates how long the current round will take. Plan to be patient and, if needed, leave the computer running. REM and TTE are usually much faster."
            ),
            tags$p(
              "For longer runs, you may prefer to download the app and run it locally. Download and setup instructions are in the ",
              tags$a(
                class = "text-link",
                href = "https://github.com/nationalparkservice/DEER_App",
                target = "_blank",
                "GitHub repository",
                .noWS = "after"
              ),
              "."
            ),

            tags$h3("Step 4: Compare results"),
            tags$p(
              "Use the ",
              tags$strong("Compare & combine"),
              " tab to:"
            ),
            tags$ul(
              tags$li("Compare models fit to uploaded field data using WAIC, and combine their density estimates (model averaging)."),
              tags$li("Check each model's posterior predictive check (Bayesian p-value)."),
              tags$li("Download posterior summaries with parameter names, means, medians, and 95% credible intervals.")
            ),

            tags$h3("Step 5: Optional design tools"),
            tags$p(
              "The ",
              tags$strong("Simulate data"),
              " tab can be used as a design tool. Default camera spacing and camera counts used in simulations are based on the ",
              tags$a(
                class = "text-link",
                href = "https://irma.nps.gov/DataStore/Reference/Profile/2307775",
                target = "_blank",
                "northeastern National Park Service deer monitoring program",
                .noWS = "after"
              ),
              "."
            ),
            tags$ul(
              tags$li("It works as a power analysis: choose a study design (camera grid, spacing, survey length, true density), simulate data, and see whether each model can recover the known true values."),
              tags$li("Each model is fit to data simulated from its own likelihood, so the three datasets differ even though they share the same design and true density."),
              tags$li("Results are shown as estimates with 95% credible intervals next to the true values, so you can judge whether a design is precise enough, or whether a model cannot recover the parameters under that design.")
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
                tags$strong("Unmarked Spatial Capture-Recapture (USCR)"),
                " — ",
                tags$em("Chandler & Royle, 2013", .noWS = "after"),
                ".",
                tags$br(),
                tags$strong("Spatial pattern → density."),
                " Where detections fall across the array informs density estimates inside the defined state space."
              ),
              tags$span(
                class = "ref",
                "Reference: ",
                tags$a(
                  class = "text-link",
                  href = "https://doi.org/10.1214/12-AOAS610",
                  target = "_blank",
                  "Chandler, R.B. & Royle, J.A. (2013)"
                )
              )
            ),
            tags$div(
              class = "about-card",
              tags$span(class = "tag tag-rem", "REM"),
              tags$p(
                tags$strong("Random Encounter Model (REM)"),
                " — ",
                tags$em("Rowcliffe et al., 2008", .noWS = "after"),
                ".",
                tags$br(),
                tags$strong("Events per unit time → density."),
                " Encounter rates are converted to density using animal speed and the camera viewshed."
              ),
              tags$span(
                class = "ref",
                "Reference: ",
                tags$a(
                  class = "text-link",
                  href = "https://doi.org/10.1111/j.1365-2664.2008.01473.x",
                  target = "_blank",
                  "Rowcliffe, J.M. et al. (2008)"
                )
              )
            ),
            tags$div(
              class = "about-card",
              tags$span(class = "tag tag-tte", "TTE"),
              tags$p(
                tags$strong("Time-to-Event (TTE)"),
                " — ",
                tags$em("Moeller et al., 2018", .noWS = "after"),
                ".",
                tags$br(),
                tags$strong("Waiting time → density."),
                " Shorter effective time between encounters implies higher density when speed and viewshed are held fixed."
              ),
              tags$span(
                class = "ref",
                "Reference: ",
                tags$a(
                  class = "text-link",
                  href = "https://doi.org/10.1002/ecs2.2331",
                  target = "_blank",
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
              tags$strong("animals per square kilometer (km²)"),
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
            tags$h3("Model development"),
            tags$p(
              "The underlying models and code were created by ",
              tags$strong("Dr. Amanda Van Buskirk"),
              " with guidance from ",
              tags$a(
                class = "text-link",
                href = "https://www.davis.wvu.edu/faculty-staff/directory/christopher-rota",
                target = "_blank",
                tags$strong("Dr. Christopher Rota")
              ),
              " at ",
              tags$strong("West Virginia University, School of Natural Resources and the Environment", .noWS = "after"),
              "."
            ),
            tags$p(
              "Scientific collaboration and feedback on model integration and quality control also came from ",
              tags$strong("Dr. Laura C. Gigliotti", .noWS = "after"),
              ", ",
              tags$strong("U.S. Geological Survey, West Virginia Cooperative Fish and Wildlife Research Unit, West Virginia University", .noWS = "after"),
              "."
            )
          ),
          tags$div(
            class = "about-card",
            tags$h3("Shiny app development"),
            tags$p(
              "This Shiny application was primarily developed by ",
              tags$strong("PhD Candidate Kacie Ring"),
              " as part of the ",
              tags$a(
                class = "text-link",
                href = "https://esa.org/programs/scip/",
                target = "_blank",
                tags$strong("SciComm in the Parks fellowship"),
                .noWS = "after"
              ),
              ". ",
              tags$strong("Dr. Christopher Rota"),
              " and ",
              tags$strong("Dr. Amanda Van Buskirk"),
              " also contributed edits to the app."
            ),
            tags$p(tags$strong("Fellowship support:")),
            tags$ul(
              tags$li(tags$strong("Dr. Brian Mitchell"), " (NPS) — Fellowship Liaison"),
              tags$li(tags$strong("Jasjeet Dhanota"), " (ESA) — Mentor"),
              tags$li(tags$strong("Mary Joy Mulumba"), " (ESA) — Program Assistant")
            )
          )
        ),

        tags$div(class = "divider")
      )
    ),

    # ---------------------- ADD YOUR DATA -------------------------
    nav_panel(
      "Add your data",
      markdown(paste(
        "Upload a **deployment CSV** (where and when cameras were set and recording) and an **images CSV** (timestamps, species, counts, and Cluster IDs).",
        "The upload checker expects the exact column names listed below.",
        sep = "\n"
      )),
      h3("The app will:"),
      markdown(paste(
        "1. Check required columns and data types.",
        "2. Flag image timestamps that fall outside each camera's deployment dates.",
        "3. Set each camera's analysis window from its deployment dates and times, optionally trimmed to a set number of days.",
        "4. Provide downloadable CSVs of the data with any corrections.",
        "",
        "After the data are uploaded and checked, review them in the **Data summary** tab.",
        "Then choose the species to analyze, and any non-default priors, in the **Model settings** tab,",
        "and run the models in the **USCR model**, **REM model**, and **TTE model** tabs.",
        "",
        "---",
        "",
        "## Required data columns",
        "",
        "### Deployment file columns (**bold = required**, others optional for quality control)",
        "",
        "- **`Site Name`**: Camera site identifier. Must be unique across all rows in the deployment file.",
        "- `Site`: Optional broader site identifier.",
        "- `Camera ID`: Optional camera identifier.",
        "- `SD Card ID`: Optional SD card identifier.",
        "- **`Start Date`**: Deployment start date (`MM/DD/YYYY`).",
        "- **`Start Time`**: Deployment start time (`HH:MM` or `HH:MM:SS`, 24-hour). If blank or unreadable, the camera's analysis window starts at midnight at the end of Start Date, so no time before the camera was set is counted.",
        "- **`End Date`**: Deployment end date (`MM/DD/YYYY`).",
        "- **`End Time`**: Deployment end time (`HH:MM` or `HH:MM:SS`, 24-hour). If blank or unreadable, the analysis window ends at midnight at the start of End Date, so no time after the camera was picked up is counted.",
        "- **`Latitude`**: Camera latitude (decimal degrees, WGS84, the standard datum used by GPS devices and phones). The app corrects the sign to match the hemisphere selected in **Model settings** (under Default camera geometry).",
        "- **`Longitude`**: Camera longitude (decimal degrees, WGS84). The app corrects the sign to match the hemisphere selected in **Model settings** (under Default camera geometry).",
        "- `Camera Model`: Optional camera make and model, for recordkeeping (for example, `Browning Strike Force Pro`).",
        "- `Camera Detection Angle`: Optional full detection angle in degrees for that camera (greater than 0 and at most 360). If left blank, the app uses the default angle from **Model settings**.",
        "- `Camera Height`: Optional camera height (meters).",
        "- `Camera Orientation`: Optional cardinal direction (`N`, `NE`, `E`, `SE`, `S`, `SW`, `W`, or `NW`) or 0-359 degrees.",
        "- **`Camera Functioning`**: Camera status; `Yes/No`, `TRUE/FALSE`, `T/F`, and `1/0` values are accepted.",
        "- **`Camera Malfunction Date`**: Leave blank for cameras that did not malfunction. If the camera failed (stopped working, or was knocked out of position but kept taking pictures), set Camera Functioning to `No` and enter the **first date whose data should be discarded** (`MM/DD/YYYY`). The camera's analysis window ends at midnight at the start of this date, so images from this date onward are excluded.",
        "- **`Detection Distance`**: Effective detection radius in meters.",
        "",
        "### Images file columns (all required unless marked optional)",
        "",
        "- **`Site Name`**: Camera site identifier. Must exactly match a Site Name in the deployment file.",
        "- **`Timestamp`**: Image date and time, written as month/day/year or year-month-day, with an optional 24-hour time (for example, `MM/DD/YYYY HH:MM:SS`, `YYYY-MM-DD HH:MM`, or `2/3/25 14:05`). Two-digit years are accepted.",
        "- **`Species`**: Species identifier. Pipe-delimited values such as `deer|squirrel` are allowed for multi-species rows; their order must match `Sighting Count`. A blank Species (for example, a trigger with no animal) is allowed, but that row isn't counted for any species.",
        "- **`Cluster ID`**: Identifier for an independent detection event (a \"cluster\"). Any format is accepted (text or numeric). Required on every row with a Species. Multiple rows can share a Cluster ID when they belong to the same event, for example burst photos, or a multi-species image split into one row per species.",
        "- **`Sighting Count`**: Number of individuals in the image. Pipe-delimited values such as `1|3` are allowed for multi-species rows, in the same order as `Species`. Required on every row with a Species.",
        "- `Image URL`: Optional image reference or link, stored as plain text, one value per row.",
        "",
        "### Notes",
        "",
        "- In the check logs, ❌ marks a problem that blocks the upload, ⚠️ marks a warning to review that doesn't block it, and ℹ️ marks a note about data the models will leave out.",
        "- **The deployment upload is blocked** if any of the following is true. These fields set each camera's analysis window and location, which every model depends on.",
        "  - A required column is missing.",
        "  - A Site Name is duplicated, or a row has deployment data but no Site Name.",
        "  - Start Date or End Date is missing or not in `MM/DD/YYYY` format, or End Date is before Start Date.",
        "  - A camera's analysis window (below) has no time in it.",
        "  - Latitude or Longitude is missing, non-numeric, or zero.",
        "  - Camera Functioning is `No` and Camera Malfunction Date is missing, not in `MM/DD/YYYY` format, on or before Start Date, or after End Date.",
        "- **The images upload is blocked** if any of the following is true.",
        "  - No deployment file has been uploaded and checked.",
        "  - A required column is missing.",
        "  - A Site Name is blank or doesn't exactly match a Site Name in the deployment file.",
        "  - A Timestamp is missing or can't be read.",
        "  - A row with a Species has a missing or non-numeric Sighting Count, or a missing Cluster ID.",
        "  - A row's Species and Sighting Count have different numbers of pipe-delimited values.",
        "- Other problems, such as a missing Camera Functioning or Detection Distance value, or image dates outside a camera's deployment dates, are reported as warnings.",
        "- **Analysis window:** each camera contributes the time from its Start Date and Start Time to its End Date and End Time (or to midnight at the start of its Camera Malfunction Date), shortened by the trim if it's on. A detection event counts only if its first image is inside the window; the images check log lists events left out this way (other than by the trim), since they can point to a wrong Start or End Time. Camera effort is the length of the window in days, minus the time spent within detection events (from first to last image of each event).",
        "- Uploading a new deployment file clears the images file; upload the images file again afterward.",
        "- Anything that changes the data the models are fit to (uploading a file, the hemisphere settings, the species to analyze, or the deployment trim) clears existing model results on uploaded data, so results from different data are never compared or combined. Changing priors or other model settings doesn't.",
        "- Changing the hemisphere settings in **Model settings** re-runs the checks on the uploaded files.",
        "- Latitude/Longitude sign correction only fixes the sign; a value of `0` is still blocked, since it can't be a real coordinate.",
        "- A Latitude/Longitude pair shared by different Site Names is flagged as a warning: two different sites shouldn't normally share the exact same coordinates.",
        "- A `Cluster ID` that appears at more than one `Site Name` is flagged as a warning: it usually means Cluster IDs aren't unique across cameras in your export.",
        "- Cross-year winter surveys (for example, `12/2025` to `01/2026`) are supported; the app uses the actual deployment dates and image timestamps, so no separate survey-year field is needed.",
        "- The app expects the exact column names listed above. Files from other tagging workflows can usually be renamed to match.",
        sep = "\n"
      )),

      h3("Step 1: Deployment file"),
      fileInput("deployment_csv", "Upload deployment CSV", accept = ".csv"),
      checkboxInput(
        "apply_trim",
        "Trim each camera's analysis window to a set number of days (24-hour periods) from its start date and time.",
        value = TRUE
      ),
      conditionalPanel(
        "input.apply_trim",
        numericInput("trim_days", "Number of days to keep", value = 56, min = 1, step = 1)
      ),
      h4("Deployment check log"),
      verbatimTextOutput("deployment_check_log"),
      h4("Preview of cleaned deployment data"),
      DTOutput("deployment_preview"),
      downloadButton("download_deployment_checked", "Download cleaned deployment CSV"),
      h4("Analysis window for each camera"),
      tags$p(
        class = "small",
        "The start and end date-times used in the analysis, after any missing-time, Camera Malfunction Date, and trimming adjustments. Days is the length of the window, before time spent within detection events is removed."
      ),
      DTOutput("deployment_windows"),

      hr(),

      h3("Step 2: Images file"),
      fileInput("images_csv", "Upload images CSV", accept = ".csv"),
      h4("Images check log"),
      verbatimTextOutput("images_check_log"),
      h4("Preview of cleaned images data"),
      DTOutput("images_preview"),
      downloadButton("download_images_checked", "Download cleaned images CSV"),

      hr()
    ),


    # ---------------------- DATA SUMMARY --------------------------
    nav_panel(
      "Data summary",
      p(
        class = "small",
        style = "color: var(--muted);",
        "This tab summarizes uploaded deployment and image data. Simulated datasets are shown in the ",
        tags$strong("Simulate data"),
        " tab."
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

    # ---------------------- MODEL SETTINGS -------------------------
    nav_panel(
      "Model settings",

      h4("Species to analyze"),
      "First, you must select the species to analyze.",
      selectInput(
        inputId = "species_to_analyze",
        label = "Select species to analyze",
        choices = character(0)
      ),
      tableOutput("species_count_summary"),
      hr(),

      markdown(paste(
        "### Run defaults",
        "",
        "**Current app defaults**:",
        "",
        "- **Number of chains**: 2 by default; the app enforces a minimum of 2 chains for model runs",
        "- **REM/TTE**: adaptive runs start from 6000 iterations, 1000 burn-in, and thin = 5",
        "- **USCR**: adaptive runs start from 6000 iterations, 1000 burn-in, thin = 5, and M = 100",
        "- **Convergence criteria**: R̂ < 1.1 for all monitored parameters",
        "- **Maximum adaptive rounds**: 6 by default. Only the settings tied to a failed check change between rounds: iterations, burn-in, and thinning double if R̂ fails; D_max (REM, TTE) or M (USCR) doubles if the posterior presses against it; and for USCR the state-space buffer grows if it is too small. If a model reaches this limit without passing every check, increase it under **Adaptive reruns** below and rerun.",
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
        "### Default camera geometry",
        "",
        "Default detection angle is 55° (based on Browning-style camera specifications). This value is used for uploaded cameras that do not include a `Camera Detection Angle` value. The **Simulate data** tab has its own detection angle setting.",
        "",
        "Longitude and Latitude sign are auto-corrected to match the selected hemispheres below during upload QC (for example, a positive Longitude gets flipped negative for Western sites, and a positive Latitude gets flipped negative for Southern sites).",
        "",
        sep = "\n"
      )),

      h4("Default camera geometry"),
      fluidRow(
        column(12,
          sliderInput("theta", "Default detection angle θ (degrees)",
                      min = 20, max = 80, value = 55, step = 1),
          radioButtons(
            "hemisphere", "Hemisphere (for Longitude sign auto-correction)",
            choices  = c("Western" = "Western", "Eastern" = "Eastern"),
            selected = "Western",
            inline   = TRUE
          ),
          radioButtons(
            "lat_hemisphere", "Hemisphere (for Latitude sign auto-correction)",
            choices  = c("Northern" = "Northern", "Southern" = "Southern"),
            selected = "Northern",
            inline   = TRUE
          )
        )
      ),

      hr(),

      h4("Adaptive reruns"),
      fluidRow(
        column(12,
          numericInput("max_adapt_rounds", "Maximum adaptive rounds per model run",
                       value = 6, min = 1, step = 1),
          tags$p(
            class = "small",
            style = "color: var(--muted);",
            "After each round, the app checks convergence (R\u0302 < 1.1), whether the upper limit on density or abundance is large enough (D_max for REM and TTE; M for USCR), and, for USCR, whether the state-space buffer is large enough. Only the settings tied to a failed check change before the next round: iterations, burn-in, and thinning double if R\u0302 fails; D_max or M doubles if it is too small; and the buffer grows to fit the posterior \u03c3 if it is too small. If a model reaches this limit without passing every check, the run-status panel will say so. Increase this number and run the model again to keep going."
          )
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
            tags$li("USCR starts from 6000 iterations, 1000 burn-in, thin = 5, and M = 100."),
            tags$li("The app checks convergence internally and can rerun models when needed."),
            tags$li("Most users only need to change the default camera angle above or switch to Advanced for species- or site-specific priors.")
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
            class = "prior-details",
            tags$summary(tags$strong("Priors: Density (all models)")),
            markdown(paste(
              "One prior on density is shared by USCR, REM, and TTE, for uploaded and simulated data.",
              "",
              "**Default (uniform) priors.** Each model's default prior is flat, but the models put it on different quantities, with different upper limits:",
              "",
              "- **REM and TTE** put it directly on density: D ~ Uniform(0, D<sub>max</sub>). D<sub>max</sub> starts at 200 animals/km², and the app doubles it and reruns if the posterior for density presses against it.",
              "- **USCR** puts it on abundance (N, the number of animals in the state space): N is uniform from 0 to M, and density is N divided by the state-space area. M starts at 100 (higher for shared-design simulations) and doubles if the posterior for N presses against it.",
              "",
              "So the default priors are all flat, but they are not the same prior.",
              "",
              "**Informative prior.** Uncheck the box to give all three models a prior with the same 95% interval for density. REM and TTE use D = D<sub>max</sub> × p, with p ~ Beta(a, b); USCR uses N ~ Binomial(M, ψ), with ψ ~ Beta(a, b). For each model, the app solves for the Beta parameters a and b so that the prior on density matches your interval, so the three priors nearly coincide. The prior's mean and median follow from the interval; for strongly skewed intervals (an upper bound far above the lower), the mean may differ from what you have in mind, so check the readout and plot before running.",
              sep = "\n"
            )),
            checkboxInput("D_prior_uniform", "Use each model's default uniform prior on density", value = TRUE),
            conditionalPanel(
              "!input.D_prior_uniform",
              fluidRow(
                column(6, numericInput("D_prior_lo", "Lower 2.5% (animals/km²)", value = NA, min = 0)),
                column(6, numericInput("D_prior_hi", "Upper 97.5% (animals/km²)", value = NA, min = 0))
              ),
              uiOutput("D_prior_readout"),
              plotOutput("D_prior_plot", height = "320px")
            )
          ),

          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: Camera-level noise (all models)")),
            markdown(paste(
              "By default, all three models assume camera counts are Poisson. Turning on camera-level noise makes them Poisson-lognormal: each camera's log expected count gets a random effect ε<sub>j</sub> ~ Normal(0, σ<sub>ε</sub>), which absorbs camera-to-camera variation in counts that the rest of the model doesn't explain (overdispersion).",
              "",
              "The setting applies to USCR, REM, and TTE, for uploaded and simulated data. Turning on camera-level noise can often fix a failed posterior predictive check (a Bayesian p-value near 0), which is reported after each run.",
              sep = "\n"
            )),
            checkboxInput("camera_noise", "Include camera-level noise (Poisson-lognormal counts)", value = FALSE),
            conditionalPanel(
              "input.camera_noise",
              markdown(paste(
                "- **σ<sub>ε</sub>** ~ Gamma(shape, rate): Prior on the standard deviation of the camera-level random effect. The default `Gamma(1, 1)` prior has mean 1 and variance 1.",
                sep = "\n"
              )),
              fluidRow(
                column(6,
                  numericInput("sd_eps_shape", HTML("σ<sub>ε</sub> gamma shape"),
                               value = advanced_prior_defaults$sd_eps_shape, min = 0.1, step = 0.1)
                ),
                column(6,
                  numericInput("sd_eps_rate", HTML("σ<sub>ε</sub> gamma rate"),
                               value = advanced_prior_defaults$sd_eps_rate, min = 0.1, step = 0.1)
                )
              )
            )
          ),

          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: REM & TTE")),
            markdown(paste(
              "- **log(v)** ~ Normal(mean, SD): Informative prior on daily movement speed (log scale). The current default corresponds to an average daily movement rate of `3.09 km/day` with an approximate 95% interval of `1.60` to `6.00 km/day` (Christensen et al. in prep).",
              "",
              "These defaults are not well calibrated for species other than white-tailed deer (*Odocoileus virginianus*). If you are working with another species, adjust the priors to match that species and study design. Model output is sensitive to choice of prior distribution.",
              "",
              sep = "\n"
            )),
            fluidRow(
              column(6,
                numericInput("log_v_mean", "log(v) mean (km/day)",
                             value = advanced_prior_defaults$log_v_mean, step = 0.1)
              ),
              column(6,
                numericInput("log_v_sd", "log(v) SD",
                             value = advanced_prior_defaults$log_v_sd, min = 0.01, step = 0.05)
              )
            )
          ),

          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: USCR")),
            markdown(paste(
              "- **log(σ)** ~ Normal(mean, SD): Informative prior on animal space use (log scale). The current default corresponds to an average 95% circular home-range size of `0.89 km²` with an approximate 95% interval of `0.48` to `1.62 km²` (Christensen et al. in prep).",
              "- **log(λ₀)** ~ Normal(mean, SD): Prior on expected detections per day when an animal's activity center is at the camera. The default mean of `0` corresponds to about one expected detection per deployed day at the detector because `exp(0) = 1`.",
              "- **State-space buffer**: Distance added around the camera array when defining the USCR state space. Animals with home range centroids outside the buffer should not be able to be detected within the state space. The default is set from the log(σ) prior so it comfortably covers a 99% home range, and updates automatically if you change that prior; the app also checks it against the posterior sigma after fitting and grows it (and reruns) if needed.",
              "",
              "These defaults are not well calibrated for species other than white-tailed deer (*Odocoileus virginianus*). If you are working with another species, adjust the priors to better match that species and study design. Model output is sensitive to choice of prior distribution.",
              "",
              sep = "\n"
            )),
            fluidRow(
              column(6,
                numericInput("log_sigma_mean", "log(σ) mean",
                             value = advanced_prior_defaults$log_sigma_mean, step = 0.1),
                numericInput("log_sigma_sd", "log(σ) SD",
                             value = advanced_prior_defaults$log_sigma_sd, min = 0.01, step = 0.05),
                numericInput("log_lam0_mean", "log(λ₀) mean",
                             value = advanced_prior_defaults$log_lam0_mean, step = 0.1)
              ),
              column(6,
                numericInput("uscr_buffer_m", "USCR state-space buffer (m)",
                             value = 1000, min = 50, step = 50),
                numericInput("log_lam0_sd", "log(λ₀) SD",
                             value = advanced_prior_defaults$log_lam0_sd, min = 0.1, step = 0.1)
              )
            )
          )
        )
      ),

      # Hidden run settings that aren't exposed in either Simple or
      # Advanced mode (chain count and MCMC iteration/burn-in/thin/M
      # starting points are managed automatically by the app).
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
          numericInput("M_uscr", NULL, value = 100)
        )
      )
    ),

    # ---------------------- USCR MODEL TAB ------------------------
    nav_panel(
      "USCR model",
      tags$div(
        id = "uscr-content",
        HTML('
        <h2>Model 1 — USCR (Unmarked Spatial Capture-Recapture)</h2>
        <p><strong>The gist:</strong> Estimates <strong>animal density in the site/state space</strong> by learning from <strong>where and when</strong> animals are detected <strong>across a camera array</strong>.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Per-camera observation model:</p>
            <p>$$y_j \\sim \\mathrm{Poisson}(\\mu_j), \\qquad \\log \\mu_j = \\log(\\mathrm{days}_j) + \\log(\\Lambda_j) + \\epsilon_j$$</p>
            <p>Expected encounter rate from augmented individuals:</p>
            <p>$$\\Lambda_j = \\sum_{i=1}^{M} z_i\\,\\lambda_0\\,\\exp\\!\\Big(-\\frac{d_{ij}^2}{2\\sigma^2}\\Big), \\qquad z_i \\sim \\mathrm{Bernoulli}(\\psi), \\qquad \\psi \\sim \\mathrm{Beta}(a, b)$$</p>
            <p>Activity centers and camera effects:</p>
            <p>$$s_i \\sim \\mathrm{Uniform}(S), \\qquad \\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>Here \\(S\\) is the buffered state space. Density is obtained from \\(N = \\sum_i z_i\\) divided by the study-area size. For both uploaded and simulated data, \\(S\\) is the union of circles, with radius equal to the buffer, around each camera, and its area is the area of that union. The state-space buffer can be adjusted in <strong>Model settings &gt; Advanced &gt; Priors: USCR</strong>. The app currently reports \\(D_{\\mathrm{km}^2}\\).</p>
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
              <li>\\(M\\), the upper bound on total abundance used for data augmentation. It starts at 100 (or about twice the upper bound on abundance implied by an informative density prior). If the posterior for \\(N\\) presses against \\(M\\), the app doubles it and reruns the model.</li>
              <li>\\(\\psi \\sim \\mathrm{Beta}(a, b)\\). The default \\(a = b = 1\\) (\\(\\psi \\sim \\mathrm{Uniform}(0, 1)\\)) gives a uniform prior on \\(N\\) from 0 to \\(M\\). An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) chooses \\(a\\) and \\(b\\) so that density \\(N / \\mathrm{area}\\) has the requested 95% prior interval.</li>
              <li>\\(\\log \\sigma \\sim \\mathcal{N}(-1.5269,\\,0.1535)\\), which corresponds to an average 95% circular home-range size of about 0.89 km² (Christensen et al. in prep). The state-space buffer defaults to 1.05&times; the 99% home-range radius implied by this prior, and the app grows it to 1.01&times; the required radius (and reruns) if the posterior sigma would need a bigger buffer.</li>
              <li>\\(\\log \\lambda_0 \\sim \\mathcal{N}(0,\\,1)\\), centered near one expected detection per deployed day when the activity center is at the camera.</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(1,1)\\), giving the camera-level random-effect standard deviation a prior mean of 1 and variance of 1. The random effect \\(\\epsilon_j\\) is included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise</strong> (off by default).</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
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
      uiOutput("run_note_uscr_nps"),
      div(
        actionButton("run_uscr_nps", "Run USCR on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      verbatimTextOutput("uscr_nps_text"),
      h5("Run status & troubleshooting"),
      verbatimTextOutput("uscr_nps_debug")
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
            <p>$$\\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>The app currently reports density as \\(D_{\\mathrm{km}^2}\\).</p>
            <h3>Variables</h3>
            <ul>
              <li>\\(y_j\\) — number of independent animal detection events at camera \\(j\\).</li>
              <li>\\(\\mathrm{days}_j\\) — deployed days for camera \\(j\\).</li>
              <li>\\(D\\) — animal density in animals/km².</li>
              <li>\\(v\\) — animal movement speed (km/day).</li>
              <li>\\(r_j\\) — effective detection radius for camera \\(j\\) (km).</li>
              <li>\\(\\theta_{\\mathrm{rad},j}\\) — full detection angle used for camera \\(j\\) (radians).</li>
              <li>\\(\\epsilon_j\\) — camera-level random effect for overdispersed counts, included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise</strong> (off by default).</li>
            </ul>
            <h3>Default priors/inputs used in the app</h3>
            <ul>
              <li>\\(D = D_{\\max}\\,p_D,\\; p_D \\sim \\mathrm{Beta}(a, b)\\). By default \\(a = b = 1\\), i.e. \\(D \\sim \\mathcal{U}(0, D_{\\max})\\) with \\(D_{\\max}\\) starting at 200 animals/km\u00b2; the app doubles \\(D_{\\max}\\) and reruns if the posterior for \\(D\\) presses against it. An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) sets \\(D_{\\max}\\) to about twice its upper bound and chooses \\(a\\) and \\(b\\) to match its 95% interval.</li>
              <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), which corresponds to an average daily movement rate of about 3.09 km/day (Christensen et al. in prep).</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(1, 1)\\), the default prior on the standard deviation of the camera-level lognormal random effect (used only when camera-level noise is on).</li>
              <li>Uploaded field data can provide a camera-specific `Camera Detection Angle`; otherwise the app uses the default angle from Model settings.</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
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
      uiOutput("run_note_rem_nps"),
      div(
        actionButton("run_rem_nps", "Run REM on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      verbatimTextOutput("rem_nps_text"),
      h5("Run status & troubleshooting"),
      verbatimTextOutput("rem_nps_debug")
    ),

    # ---------------------- TTE MODEL TAB -------------------------
    nav_panel(
      "TTE model",
      tags$div(
        id = "tte-content",
        HTML('
        <h2>Model 3 — TTE (Time-to-Event)</h2>
        <p><strong>The gist:</strong> Uses <strong>animal detection events per camera</strong>, scaled by <strong>movement-based time units</strong> and the <strong>viewshed area</strong>. Shorter effective time between encounters implies higher density.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Per-camera Poisson model:</p>
            <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j \\;=\\; \\log D \\;+\\; \\log(\\mathrm{days}_j)\\;+\\;\\log(U_j)\\;+\\;\\log(A_j)\\;+\\;\\epsilon_j$$</p>
            <p>Movement-based encounter multiplier, viewshed area, and camera effect:</p>
            <p>$$U_j \\;=\\; \\frac{v}{f(\\theta_j)\\,r_j}, \\qquad A_j \\;=\\; \\pi r_j^2 \\frac{\\theta_{\\mathrm{deg},j}}{360}, \\qquad \\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>with \\(f(\\theta_j) = 0.3324 + 0.005580\\,\\theta_j - 1.454\\times10^{-5}\\,\\theta_j^2\\), where \\(\\theta_j\\) is the full camera detection angle in degrees.</p>
            <p>The app currently reports density as \\(D_{\\mathrm{km}^2}\\).</p>
            <h3>Variables</h3>
            <ul>
              <li>\\(y_j\\) — number of animal detection events for camera \\(j\\).</li>
              <li>\\(\\mathrm{days}_j\\) — total deployed days for camera \\(j\\).</li>
              <li>\\(U_j\\) — the time needed to move across the average camera viewshed width \\(j\\).</li>
              <li>\\(v\\) — movement speed (km/day); prior on \\(\\log v\\) as below.</li>
              <li>\\(r_j\\) — effective detection radius for camera \\(j\\) (km).</li>
              <li>\\(\\theta_{\\mathrm{deg},j}\\) — full detection angle used for camera \\(j\\) (degrees).</li>
              <li>\\(A_j\\) — viewshed area for camera \\(j\\) (km²).</li>
              <li>\\(\\epsilon_j\\) — camera-level random effect for overdispersed counts, included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise</strong> (off by default).</li>
            </ul>
            <h3>Default priors/inputs used in the app</h3>
            <ul>
              <li>\\(D = D_{\\max}\\,p_D,\\; p_D \\sim \\mathrm{Beta}(a, b)\\). By default \\(a = b = 1\\), i.e. \\(D \\sim \\mathcal{U}(0, D_{\\max})\\) with \\(D_{\\max}\\) starting at 200 animals/km\u00b2; the app doubles \\(D_{\\max}\\) and reruns if the posterior for \\(D\\) presses against it. An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) sets \\(D_{\\max}\\) to about twice its upper bound and chooses \\(a\\) and \\(b\\) to match its 95% interval.</li>
              <li>\\(\\log v \\sim \\mathcal{N}(1.130,\\,0.3372)\\), which corresponds to an average daily movement rate of about 3.09 km/day (Christensen et al. in prep).</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(1, 1)\\), the default prior on the standard deviation of the camera-level lognormal random effect (used only when camera-level noise is on).</li>
              <li>Uploaded field data can provide a camera-specific `Camera Detection Angle`; otherwise the app uses the default angle from Model settings.</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
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
      uiOutput("run_note_tte_nps"),
      div(
        actionButton("run_tte_nps", "Run TTE on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      verbatimTextOutput("tte_nps_text"),
      h5("Run status & troubleshooting"),
      verbatimTextOutput("tte_nps_debug")
    ),

    # ---------------------- COMPARE & COMBINE --------------------
    nav_panel(
      "Compare & combine",
      markdown(paste(
        "**Uploaded field data:** the table below updates as REM, TTE, and USCR finish. If only some models have completed, the table will still summarize the completed fits.",
        "",
        "1. Compute ΔWAIC and WAIC weights when multiple models are available and WAIC is available for all completed fits;",
        "2. Combine posterior draws of density (animals/km²) across the completed fits;",
        "3. Report model-specific densities, and when possible also report unweighted, WAIC-weighted, and mixture distribution summaries plus the probability that density exceeds the user-set threshold. If WAIC is missing for one or more completed fits, only the unweighted summary is shown.",
        "",
        "Run the models from their tabs first, then check the tables.",
        sep = "\n"
      )),
      numericInput(
        "combo_density_threshold",
        "Density threshold for probability summary (animals/km²): Enter a density value below. The app will calculate the probability that the estimated density exceeds this threshold for each individual model and for the weighted and unweighted combined estimates.",
        value = 8,
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
      posterior_param_legend()
    ),

    # ---------------------- SIMULATE DATA -------------------------
    nav_panel(
      "Simulate data",
      markdown(paste(
        "Use this tab as a **power analysis** to see how study design affects the precision and accuracy of density estimates.",
        "Choose a study design (camera grid, spacing, survey length, true density) and model-specific parameter values, then click **Simulate shared design**.",
        "Each model block below simulates its own dataset under that model's assumptions: USCR from a spatial capture-recapture process, and REM and TTE from a movement and encounter process.",
        "The three datasets therefore differ, even though they share the same design and true density.",
        "",
        "Run each model separately in its block. Each block reports estimates with 95% credible intervals alongside the true values.",
        "If a model cannot recover the true values under a design, that design may be inadequate for that method.",
        "",
        "**Interpreting results:** each click simulates one dataset, so a single miss or hit is not conclusive. Try several random seeds before drawing conclusions.",
        "Model priors for σ, λ₀, and movement speed are centered on the true values you set here, so results show how well a design can recover parameters when priors are well calibrated.",
        "The prior SDs for these parameters are not set here: adjust them in the **Model settings** tab under **Advanced** (log(σ) SD and log(λ₀) SD in **Priors: USCR**; log(v) SD in **Priors: REM & TTE**).",
        "Priors for density, σ<sub>ε</sub>, and MCMC settings come from the **Model settings** tab (the density prior for all three models is under **Advanced** → **Priors: Density (all models)**).",
        sep = "\n"
      )),
      h3("Shared design parameters"),
      tags$p(
        class = "small",
        style = "color: var(--muted);",
        "The default values below match the app's default prior means (home-range size and σ, λ₀, movement speed, σ", tags$sub("ε", .noWS = c("before", "after")), ", and detection angle), so with default settings the simulated truth sits at the center of each model's priors."
      ),
      fluidRow(
        column(4,
               sliderInput("pshared_n_side", "Grid dimension (n x n cameras)",
                           min = 3, max = 8, value = 5, step = 1),
               sliderInput("pshared_spacing", "Camera spacing (m)",
                           min = 150, max = 500, value = 300, step = 10),
               sliderInput("pshared_days", "Survey length (days)",
                           min = 7, max = 60, value = 21, step = 1),
               sliderInput("pshared_Dtrue", "True density (animals/km²)",
                           min = 5, max = 80, value = 25, step = 1),
               conditionalPanel(
                 "input.camera_noise && input.mode == 'Advanced'",
                 sliderInput("pshared_sd_eps", HTML("Camera-level heterogeneity (σ<sub>ε</sub>)"),
                             min = 0.05, max = 2, value = 1, step = 0.05)
               ),
               conditionalPanel(
                 "!(input.camera_noise && input.mode == 'Advanced')",
                 tags$p(
                   class = "small",
                   style = "color: var(--muted);",
                   "Camera-level noise is off, so simulated counts are Poisson. Turn it on in ",
               tags$strong("Model settings > Advanced > Priors: Camera-level noise"),
               " to simulate and fit Poisson-lognormal counts."
                 )
               )
        ),
        column(4,
               h4("USCR inputs"),
               sliderInput("pshared_home_range_km2", "Mean 95% home-range size (km²)",
                           min = 0.10, max = 3.00, value = 0.89, step = 0.05),
               sliderInput("pshared_lambda0", "USCR detection rate λ₀ (expected detections/day at camera center)",
                           min = 0.05, max = 2, value = 1, step = 0.05),
               numericInput("pshared_buffer_m", "USCR State-space buffer around camera array (m)",
                            value = 1000, min = 50, step = 50),
               p(class = "small", style = "color: var(--muted);",
                 "Simulated activity centers are placed within this distance of at least one camera. By default it uses the same rule as uploaded data: 1.05 times the 99% home-range radius implied by the 99th percentile of the log(σ) prior, here centered on the σ implied by the home-range size above. It updates when you change the home-range size or the log(σ) SD, and you can override it.")
        ),
        column(4,
               h4("REM / TTE (movement) inputs"),
               sliderInput("pshared_v", "Movement speed (km/day)",
                           min = 0.5, max = 10, value = 3.09, step = 0.1),
               sliderInput("pshared_theta", "Detection angle θ (degrees)",
                           min = 20, max = 80, value = 55, step = 1),
               sliderInput("pshared_r", "Detection radius (m)",
                           min = 8, max = 25, value = 12, step = 1)
        )
      ),
      fluidRow(
        column(4, numericInput("pshared_seed", "Random seed", value = 1, min = 1)),
        column(8, br(), actionButton("pshared_run_sim", "Simulate shared design", class = "btn-primary"))
      ),
      tags$p(
        class = "small", style = "color: var(--muted);",
        "Note: USCR is simulated from a spatial (SECR) process; REM/TTE are simulated from a movement process. ",
        "Both use the design parameters and true density above, so they represent the same study design ",
        "under each model's own generative assumptions."
      ),
      tags$details(
        class = "about-card",
        style = "margin-top: 1rem;",
        tags$summary(
          style = "cursor: pointer; font-weight: 600;",
          "Simulation parameter guide"
        ),
        tags$div(style = "margin-top: 0.85rem;"),
        tags$p(tags$strong("Shared design (all models)")),
        tags$p(
          tags$strong("Grid dimension"),
          " sets the number of cameras along one side of the square grid. A 5 x 5 grid gives 25 cameras."
        ),
        tags$p(
          tags$strong("Camera spacing"),
          " is the distance between neighboring cameras. For USCR, cameras should be close enough that an individual can be detected at more than one camera."
        ),
        tags$p(
          tags$strong("Survey length"),
          " is the number of deployed days simulated for each camera."
        ),
        tags$p(
          tags$strong("True density"),
          " is the animal density used to generate the data, in animals/km². Estimates are reported in the same units."
        ),
        tags$p(
          tags$strong("Random seed"),
          " reproduces the same simulated datasets when you rerun with identical settings. Change it to see how results vary between datasets."
        ),
        tags$p(tags$strong("USCR inputs")),
        tags$p(
          tags$strong("Mean 95% home-range size"),
          " is converted in the background to the spatial scale parameter σ. Larger home ranges spread detections across more cameras."
        ),
        tags$p(
          tags$strong("Detection rate λ₀"),
          " is the expected number of detections per deployed day if an animal's activity center were at the camera. Larger values produce more detections everywhere in the array."
        ),
        tags$p(
          tags$strong("State-space buffer"),
          " is the maximum distance from the nearest camera at which simulated activity centers are placed, so the simulated state space is the union of circles of this radius around the cameras. The same buffer is used when fitting USCR. By default it follows the rule used for uploaded data (1.05 times the 99% home-range radius at the 99th percentile of the log(σ) prior, centered here on the true σ)."
        ),
        tags$p(tags$strong("REM and TTE inputs")),
        tags$p(
          tags$strong("Movement speed (km/day)"),
          " is the average daily animal movement rate. It is species-specific."
        ),
        tags$p(
          tags$strong("Detection radius"),
          " is how far a camera can detect an animal, in meters."
        ),
        tags$p(
          tags$strong("Detection angle θ"),
          " is the camera's detection angle in degrees. Manufacturer specifications typically list this value."
        ),
        tags$p(
          tags$strong(HTML("Camera-level heterogeneity (σ<sub>ε</sub>)")),
          " is the standard deviation of camera-to-camera variation in counts (on the log scale) that is not explained by deployment length, density, movement, or camera geometry. It is used only when camera-level noise is turned on in ", tags$strong("Model settings > Advanced > Priors: Camera-level noise", .noWS = "after"), "; then it is applied when simulating data for all three models, and all three models estimate ", HTML("σ<sub>ε</sub>"), " when fit."
        )
      ),
      hr(),
      h3("Fit each model to the shared design"),
      markdown(paste(
        "**Reading the results tables:**",
        "",
        "- **Truth**: the value used to simulate the data.",
        "- **Mean / Median**: posterior mean and median.",
        "- **Lower 2.5% / Upper 97.5%**: bounds of the 95% credible interval.",
        "- **Covered**: whether the 95% credible interval contains the true value. Across many simulated datasets, a well-calibrated model should cover the truth about 95% of the time, so a single \"No\" is not by itself a problem.",
        sep = "\n"
      )),
      uiOutput("pshared_noise_note"),

      # ---------------- USCR block ----------------
      tags$div(
        class = "about-card",
        h4("USCR"),
        h5("Simulated camera array"),
        plotOutput("pshared_uscr_grid_plot", height = "380px"),
        h5("Simulated camera data"),
        DTOutput("pshared_uscr_camera_table"),
        br(),
        uiOutput("run_note_uscr_pshared"),
        actionButton("pshared_run_uscr", "Run USCR", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_uscr_table"),
        h5("Run status & troubleshooting"),
        verbatimTextOutput("pshared_uscr_debug")
      ),

      # ---------------- REM block ----------------
      tags$div(
        class = "about-card",
        h4("REM"),
        h5("Simulated camera array"),
        plotOutput("pshared_rem_grid_plot", height = "380px"),
        h5("Simulated camera data"),
        DTOutput("pshared_rem_camera_table"),
        br(),
        uiOutput("run_note_rem_pshared"),
        actionButton("pshared_run_rem", "Run REM", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_rem_table"),
        h5("Run status & troubleshooting"),
        verbatimTextOutput("pshared_rem_debug")
      ),

      # ---------------- TTE block ----------------
      tags$div(
        class = "about-card",
        h4("TTE"),
        h5("Simulated camera array"),
        plotOutput("pshared_tte_grid_plot", height = "380px"),
        h5("Simulated camera data"),
        DTOutput("pshared_tte_camera_table"),
        br(),
        uiOutput("run_note_tte_pshared"),
        actionButton("pshared_run_tte", "Run TTE", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_tte_table"),
        h5("Run status & troubleshooting"),
        verbatimTextOutput("pshared_tte_debug")
      )
    )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------

  server <- function(input, output, session) {

    # make_model_debug() is used immediately below (reactiveVal() forces its
    # initial value right away, unlike observeEvent/render blocks), so it must
    # be defined here, before first use, rather than further down in the file.
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
        history = character()
      )
    }

  # ================================================================
  # POWER ANALYSIS: shared design comparison (USCR + REM + TTE)
  # ================================================================

  pshared_sim <- eventReactive(input$pshared_run_sim, {
    uscr_sim <- simulate_camera_counts(
      n_side = input$pshared_n_side,
      spacing_m = input$pshared_spacing,
      days = input$pshared_days,
      D_per_km2 = input$pshared_Dtrue,
      lambda0 = input$pshared_lambda0,
      home_range_km2 = input$pshared_home_range_km2,
      buffer_m = input$pshared_buffer_m,
      sd_eps = input$pshared_sd_eps,
      camera_noise = app_camera_noise(),
      seed = input$pshared_seed + 2  # offset so USCR, REM, and TTE use separate random streams
    )
    uscr_data <- sim_model_inputs(sim = uscr_sim, detection_radius_m = input$pshared_r)

    rem_sim <- simulate_teaching_counts(
      model = "REM",
      n_side = input$pshared_n_side, spacing_m = input$pshared_spacing, days = input$pshared_days,
      D_per_km2 = input$pshared_Dtrue, detection_radius_m = input$pshared_r,
      theta_deg = input$pshared_theta, v_km_day = input$pshared_v, sd_eps = input$pshared_sd_eps,
      camera_noise = app_camera_noise(),
      seed = input$pshared_seed
    )
    tte_sim <- simulate_teaching_counts(
      model = "TTE",
      n_side = input$pshared_n_side, spacing_m = input$pshared_spacing, days = input$pshared_days,
      D_per_km2 = input$pshared_Dtrue, detection_radius_m = input$pshared_r,
      theta_deg = input$pshared_theta, v_km_day = input$pshared_v, sd_eps = input$pshared_sd_eps,
      camera_noise = app_camera_noise(),
      seed = input$pshared_seed + 1  # offset so REM and TTE get independent camera noise
    )

    # Settings are saved with the simulated data so that the model fits
    # (priors, state-space buffer) and the truth-vs-estimate tables use the
    # values the data were simulated with, even if a slider is moved
    # afterwards without clicking 'Simulate shared design' again.
    list(uscr = list(sim = uscr_sim, data = uscr_data), rem = rem_sim, tte = tte_sim,
         D_true = input$pshared_Dtrue,
         settings = list(
           n_side   = input$pshared_n_side,
           spacing  = input$pshared_spacing,
           buffer_m = input$pshared_buffer_m,
           lambda0  = input$pshared_lambda0,
           v        = input$pshared_v,
           theta    = input$pshared_theta,
           sd_eps   = input$pshared_sd_eps,
           camera_noise = app_camera_noise(),
           # the automatic buffer at simulation time, to tell whether the
           # buffer used was a user override
           buffer_default_m = ceiling(uscr_default_buffer_m(
             log(home_range_km2_to_sigma_m(input$pshared_home_range_km2) / 1000),
             prior_input("log_sigma_sd")
           ) / 50) * 50
         ))
  }, ignoreInit = TRUE)

  # Keep the shared-design buffer default in sync with the home-range size,
  # using the same rule as the uploaded-data buffer: the 99th percentile of
  # the log(sigma) prior the shared-design fit uses (centered on the true
  # sigma, with the log(sigma) SD from Model settings). Users can override it.
  observe({
    req(is.numeric(input$pshared_home_range_km2))
    sigma_km <- home_range_km2_to_sigma_m(input$pshared_home_range_km2) / 1000
    log_sigma_sd <- prior_input("log_sigma_sd")
    req(is.numeric(log_sigma_sd), is.finite(log_sigma_sd))
    new_default <- uscr_default_buffer_m(log(sigma_km), log_sigma_sd)
    updateNumericInput(session, "pshared_buffer_m", value = ceiling(new_default / 50) * 50)
  })

  # Which noise setting the current simulated design (and so its fits) uses.
  output$pshared_noise_note <- renderUI({
    req(pshared_sim())
    sim_noise <- isTRUE(pshared_sim()$settings$camera_noise)
    toggle_noise <- app_camera_noise()
    tags$div(
      class = "about-card",
      style = "padding: 10px 14px;",
      tags$p(
        style = "margin: 0;",
        tags$strong(paste0("Current simulated design: camera-level noise ", if (sim_noise) "on" else "off", ". ")),
        "All three fits below use this setting. To change it, change the camera-level noise setting in ",
        tags$strong("Model settings > Advanced > Priors: Camera-level noise"),
        " and click ",
        tags$strong("Simulate shared design"),
        " again."
      ),
      if (sim_noise != toggle_noise) {
        tags$p(
          style = "margin: 0.4rem 0 0 0; color: firebrick;",
          paste0("The Model settings toggle is now ", if (toggle_noise) "on" else "off",
                 ". It will apply to the next simulation and to uploaded data, not to this design.")
        )
      }
    )
  })

  output$pshared_uscr_grid_plot <- renderPlot({
    req(pshared_sim())
    s <- pshared_sim()$uscr$sim
    cams_df <- s$out |> mutate(camera = row_number())
    ggplot() +
      geom_path(data = s$state_space_outline, aes(x = x, y = y, group = group),
                color = redwood_colors[5], linetype = "dashed") +
      geom_point(data = s$activity_centers, aes(x = x, y = y),
                 size = 1, alpha = 0.5, color = redwood_colors[2]) +
      geom_point(data = cams_df, aes(x = x, y = y), size = 3, color = redwood_colors[3]) +
      geom_label(data = cams_df, aes(x = x, y = y, label = camera), nudge_y = 18, size = 3.8,
                 color = redwood_colors[1], fill = "white", linewidth = 0.2) +
      coord_equal() + theme_minimal() +
      labs(x = "m", y = "m", title = sprintf(
        "%d cameras, spacing %dm, buffer %dm (true density: %.1f animals/km²)",
        s$truth$n_cams, s$truth$spacing_m, round(s$truth$buffer_m), s$truth$D_per_km2
      ), subtitle = sprintf(
        "%d simulated activity centers in a %.2f km² state space (dashed outline)",
        s$truth$N, s$truth$area_km2
      ))
  })

  output$pshared_uscr_camera_table <- renderDT({
    req(pshared_sim())
    s <- pshared_sim()$uscr$sim
    s$out |>
      mutate(camera_id = `Site Name`, total_counts = s$camera_counts) |>
      select(camera_id, x, y, total_counts) |>
      datatable(options = list(pageLength = 8))
  })

  output$pshared_rem_grid_plot <- renderPlot({
    req(pshared_sim())
    s <- pshared_sim()$rem
    df <- s$out |> mutate(total_detections = s$camera_counts)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(aes(size = total_detections), shape = 21, fill = redwood_colors[3],
                 color = redwood_colors[1], alpha = 0.75) +
      ggrepel::geom_label_repel(aes(label = `Site Name`), size = 3, max.overlaps = Inf) +
      scale_size_area(name = "Total detections", max_size = 14) +
      coord_equal() + theme_minimal() +
      labs(title = sprintf("REM simulated detections (true density: %.1f animals/km²)", s$truth$D_per_km2),
           x = "m", y = "m")
  })

  output$pshared_rem_camera_table <- renderDT({
    req(pshared_sim())
    s <- pshared_sim()$rem
    s$out |>
      mutate(total_detections = s$camera_counts, camera_days = s$camera_days) |>
      select(`Site Name`, total_detections, camera_days) |>
      datatable(options = list(pageLength = 8))
  })

  output$pshared_tte_grid_plot <- renderPlot({
    req(pshared_sim())
    s <- pshared_sim()$tte
    df <- s$out |> mutate(total_detections = s$camera_counts)
    ggplot(df, aes(x = x, y = y)) +
      geom_point(aes(size = total_detections), shape = 21, fill = redwood_colors[3],
                 color = redwood_colors[1], alpha = 0.75) +
      ggrepel::geom_label_repel(aes(label = `Site Name`), size = 3, max.overlaps = Inf) +
      scale_size_area(name = "Total detections", max_size = 14) +
      coord_equal() + theme_minimal() +
      labs(title = sprintf("TTE simulated detections (true density: %.1f animals/km²)", s$truth$D_per_km2),
           x = "m", y = "m")
  })

  output$pshared_tte_camera_table <- renderDT({
    req(pshared_sim())
    s <- pshared_sim()$tte
    s$out |>
      mutate(total_detections = s$camera_counts, camera_days = s$camera_days) |>
      select(`Site Name`, total_detections, camera_days) |>
      datatable(options = list(pageLength = 8))
  })

  pshared_uscr_fit <- reactiveVal(NULL)
  pshared_rem_fit  <- reactiveVal(NULL)
  pshared_tte_fit  <- reactiveVal(NULL)
  pshared_uscr_debug <- reactiveVal(make_model_debug("USCR", "Shared design", "Shared design comparison", ""))
  pshared_rem_debug  <- reactiveVal(make_model_debug("REM",  "Shared design", "Shared design comparison", ""))
  pshared_tte_debug  <- reactiveVal(make_model_debug("TTE",  "Shared design", "Shared design comparison", ""))

  observeEvent(input$pshared_run_sim, {
    pshared_uscr_fit(NULL); pshared_rem_fit(NULL); pshared_tte_fit(NULL)
  }, ignoreInit = TRUE)

  # ---- USCR ----
  observeEvent(input$pshared_run_uscr, {
    req(pshared_sim())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("USCR not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    s <- pshared_sim()
    d <- s$uscr$data
    sigma_km_true <- s$uscr$sim$truth$sigma_m / 1000

    update_model_debug(pshared_uscr_debug, status = "running", stage = "Running",
                       started_at = Sys.time(), guidance = "Fitting USCR on the shared design.",
                       context = c(
                         paste("Cameras:", nrow(d$out)),
                         paste("True density:", s$D_true, "animals/km^2"),
                         paste("True sigma:", round(sigma_km_true, 3), "km"),
                         paste("True lambda0:", s$settings$lambda0),
                         paste("Density prior:", describe_density_prior_spec(density_prior_spec())),
                         paste("Camera-level noise:", noise_label(s$settings$camera_noise))
                       ),
                       log_entry = "Shared-design USCR run requested.")
    showNotification("Running USCR on shared design...", type = "message", duration = NULL, id = "pshared_uscr_status")

    suggested_M <- pshared_suggested_M(s)

    fit <- tryCatch(
      run_USCR_app(
        out = d$out, camera_counts = d$camera_counts, camera_days = d$camera_days,
        iter = input$iter_uscr, burnin = input$burnin_uscr, thin = input$thin_uscr,
        n_chains = app_n_chains(), M = suggested_M,
        log_sigma_mean = log(sigma_km_true), log_sigma_sd = prior_input("log_sigma_sd"),
        log_lam0_mean = log(s$settings$lambda0), log_lam0_sd = prior_input("log_lam0_sd"),
        sd_eps_shape = prior_input("sd_eps_shape"), sd_eps_rate = prior_input("sd_eps_rate"),
        buffer_m = s$settings$buffer_m, adaptive = TRUE, compute_WAIC = FALSE,
        max_adapt_rounds = app_max_adapt_rounds(),
        density_prior = density_prior_spec(),
        camera_noise = s$settings$camera_noise,
        diagnostic_mode = FALSE, parallel_chains = isTRUE(app_n_chains() > 1),
        status_callback = make_model_status_callback(pshared_uscr_debug, "USCR", "pshared_uscr_status"),
        verbose = FALSE
      ),
      error = function(e) {
        update_model_debug(pshared_uscr_debug, status = "error", stage = "Failed", finished_at = Sys.time(),
                           raw_error = e$message, guidance = friendly_model_error("USCR", "shared design", e$message))
        showNotification(paste("USCR failed:", e$message), type = "error", duration = NULL)
        NULL
      },
      finally = removeNotification("pshared_uscr_status")
    )
    pshared_uscr_fit(fit)
    if (!is.null(fit)) {
      update_model_debug(pshared_uscr_debug, status = "success", stage = "Complete", finished_at = Sys.time(),
                         guidance = adaptive_limit_guidance(fit, "USCR", "Complete. Review the table above."))
      notify_adaptive_limit(fit, "USCR (shared design)")
      notify_ppc(fit, "USCR (shared design)")
    }
  })

  # Display version of a summarize_vs_truth() table: parameter symbols with
  # units instead of internal names, and Yes/No for coverage. The internal
  # names are still used to look up the posterior draws.
  pshared_truth_datatable <- function(tbl, fit = NULL) {
    param_labels <- c(
      D_km2  = "D (animals/km\u00b2)",
      sigma  = "\u03c3 (km)",
      lam_0  = "\u03bb\u2080 (detections/day)",
      v      = "v (km/day)",
      sd_eps = "σ<sub>ε</sub>"
    )
    tbl <- tbl |>
      mutate(
        Parameter = dplyr::coalesce(unname(param_labels[Parameter]), Parameter),
        Covered = ifelse(Covered, "Yes", "No"),
        across(where(is.numeric), ~round(.x, 3))
      )
    # Posterior predictive check under the table (red when it fails).
    ppc <- ppc_summary(fit)
    caption <- if (!is.null(ppc)) {
      htmltools::tags$caption(
        style = paste0("caption-side: bottom; text-align: left; color: ",
                       if (ppc$fail) "firebrick" else "var(--muted)", ";"),
        ppc_text(fit)
      )
    }
    # escape = FALSE so the <sub> tag in the sigma_epsilon label renders.
    datatable(tbl, caption = caption, options = list(dom = "t", paging = FALSE),
              rownames = FALSE, escape = FALSE)
  }

  output$pshared_uscr_table <- renderDT({
    fit <- pshared_uscr_fit()
    # req() rather than validate(): a validation message inside an empty
    # DT widget has no height and overlaps the heading below it.
    req(fit)
    s <- pshared_sim()
    sigma_km_true <- s$uscr$sim$truth$sigma_m / 1000
    tbl <- summarize_vs_truth(fit, list(
      D_km2 = s$D_true,
      sigma = sigma_km_true,
      lam_0 = s$settings$lambda0,
      sd_eps = if (isTRUE(s$settings$camera_noise)) s$settings$sd_eps  # row only with noise on
    ))
    pshared_truth_datatable(tbl, fit)
  })

  # ---- REM ----
  observeEvent(input$pshared_run_rem, {
    req(pshared_sim())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("REM not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    s <- pshared_sim()

    update_model_debug(pshared_rem_debug, status = "running", stage = "Running",
                       started_at = Sys.time(), guidance = "Fitting REM on the shared design.",
                       context = c(
                         paste("Cameras:", nrow(s$rem$out)),
                         paste("True density:", s$D_true, "animals/km^2"),
                         paste("True speed:", s$settings$v, "km/day"),
                         paste("Density prior:", describe_density_prior_spec(density_prior_spec())),
                         paste("Camera-level noise:", noise_label(s$settings$camera_noise))
                       ),
                       log_entry = "Shared-design REM run requested.")
    showNotification("Running REM on shared design...", type = "message", duration = NULL, id = "pshared_rem_status")

    fit <- tryCatch(
      run_REM(
        y = s$rem$camera_counts, r_km = s$rem$out$`Detection Distance` / 1000,
        camera_days = s$rem$camera_days, theta_deg = s$settings$theta,
        iter = input$iter_rem_tte, burnin = input$burnin_rem_tte, thin = input$thin_rem_tte,
        n_chains = app_n_chains(), D_max = rem_tte_start_D_max,
        log_v_mean = log(s$settings$v), log_v_sd = prior_input("log_v_sd"),
        sd_eps_shape = prior_input("sd_eps_shape"), sd_eps_rate = prior_input("sd_eps_rate"),
        max_adapt_rounds = app_max_adapt_rounds(),
        density_prior = density_prior_spec(),
        camera_noise = s$settings$camera_noise,
        status_callback = make_model_status_callback(pshared_rem_debug, "REM", "pshared_rem_status")
      ),
      error = function(e) {
        update_model_debug(pshared_rem_debug, status = "error", stage = "Failed", finished_at = Sys.time(),
                           raw_error = e$message, guidance = friendly_model_error("REM", "shared design", e$message))
        showNotification(paste("REM failed:", e$message), type = "error", duration = NULL)
        NULL
      },
      finally = removeNotification("pshared_rem_status")
    )
    pshared_rem_fit(fit)
    if (!is.null(fit)) {
      update_model_debug(pshared_rem_debug, status = "success", stage = "Complete", finished_at = Sys.time(),
                         guidance = adaptive_limit_guidance(fit, "REM", "Complete. Review the table above."))
      notify_adaptive_limit(fit, "REM (shared design)")
      notify_ppc(fit, "REM (shared design)")
    }
  })

  output$pshared_rem_table <- renderDT({
    fit <- pshared_rem_fit()
    # req() rather than validate(): a validation message inside an empty
    # DT widget has no height and overlaps the heading below it.
    req(fit)
    s <- pshared_sim()
    tbl <- summarize_vs_truth(fit, list(
      D_km2 = s$D_true, v = s$settings$v,
      sd_eps = if (isTRUE(s$settings$camera_noise)) s$settings$sd_eps  # row only with noise on
    ))
    pshared_truth_datatable(tbl, fit)
  })

  # ---- TTE ----
  observeEvent(input$pshared_run_tte, {
    req(pshared_sim())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("TTE not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    s <- pshared_sim()

    update_model_debug(pshared_tte_debug, status = "running", stage = "Running",
                       started_at = Sys.time(), guidance = "Fitting TTE on the shared design.",
                       context = c(
                         paste("Cameras:", nrow(s$tte$out)),
                         paste("True density:", s$D_true, "animals/km^2"),
                         paste("True speed:", s$settings$v, "km/day"),
                         paste("Density prior:", describe_density_prior_spec(density_prior_spec())),
                         paste("Camera-level noise:", noise_label(s$settings$camera_noise))
                       ),
                       log_entry = "Shared-design TTE run requested.")
    showNotification("Running TTE on shared design...", type = "message", duration = NULL, id = "pshared_tte_status")

    fit <- tryCatch(
      run_TTE(
        y = s$tte$camera_counts, r_km = s$tte$out$`Detection Distance` / 1000,
        camera_days = s$tte$camera_days, theta_deg = s$settings$theta,
        iter = input$iter_rem_tte, burnin = input$burnin_rem_tte, thin = input$thin_rem_tte,
        n_chains = app_n_chains(), D_max = rem_tte_start_D_max,
        log_v_mean = log(s$settings$v), log_v_sd = prior_input("log_v_sd"),
        sd_eps_shape = prior_input("sd_eps_shape"), sd_eps_rate = prior_input("sd_eps_rate"),
        max_adapt_rounds = app_max_adapt_rounds(),
        density_prior = density_prior_spec(),
        camera_noise = s$settings$camera_noise,
        status_callback = make_model_status_callback(pshared_tte_debug, "TTE", "pshared_tte_status")
      ),
      error = function(e) {
        update_model_debug(pshared_tte_debug, status = "error", stage = "Failed", finished_at = Sys.time(),
                           raw_error = e$message, guidance = friendly_model_error("TTE", "shared design", e$message))
        showNotification(paste("TTE failed:", e$message), type = "error", duration = NULL)
        NULL
      },
      finally = removeNotification("pshared_tte_status")
    )
    pshared_tte_fit(fit)
    if (!is.null(fit)) {
      update_model_debug(pshared_tte_debug, status = "success", stage = "Complete", finished_at = Sys.time(),
                         guidance = adaptive_limit_guidance(fit, "TTE", "Complete. Review the table above."))
      notify_adaptive_limit(fit, "TTE (shared design)")
      notify_ppc(fit, "TTE (shared design)")
    }
  })

  output$pshared_tte_table <- renderDT({
    fit <- pshared_tte_fit()
    # req() rather than validate(): a validation message inside an empty
    # DT widget has no height and overlaps the heading below it.
    req(fit)
    s <- pshared_sim()
    tbl <- summarize_vs_truth(fit, list(
      D_km2 = s$D_true, v = s$settings$v,
      sd_eps = if (isTRUE(s$settings$camera_noise)) s$settings$sd_eps  # row only with noise on
    ))
    pshared_truth_datatable(tbl, fit)
  })

  output$pshared_uscr_debug <- renderText({
    state <- pshared_uscr_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, pshared_uscr_fit())
  })

  output$pshared_rem_debug <- renderText({
    state <- pshared_rem_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, pshared_rem_fit())
  })

  output$pshared_tte_debug <- renderText({
    state <- pshared_tte_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, pshared_tte_fit())
  })

  # ================================================================
  # UPLOADED DATA: deployment + images QC (Add your data tab)
  # ================================================================

  # Run a QC step, collecting its messages and warnings for the check log.
  # Returns list(value, messages, warnings); value is NULL if the step stopped
  # with an error, and the error text is added to the warnings.
  run_qc_step <- function(expr) {
    msgs <- character()
    warns <- character()
    value <- tryCatch(
      withCallingHandlers(
        expr,
        message = function(m) {
          msgs <<- c(msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        warns <<- c(warns, paste0("ERROR: ", conditionMessage(e)))
        NULL
      }
    )
    list(value = value, messages = msgs, warnings = warns)
  }

  # Text for a check log from list(messages, warnings).
  format_check_log <- function(issues) {
    lines <- character()
    if (length(issues$messages)) {
      lines <- c(lines, "Messages:", paste0("  - ", trimws(issues$messages)), "")
    }
    if (length(issues$warnings)) {
      lines <- c(lines, "Warnings / issues:", paste0("  - ", trimws(issues$warnings)))
    }
    if (!length(lines)) lines <- "No issues reported."
    paste(lines, collapse = "\n")
  }

  # Deployment-length trim (days), or NULL when trimming is off.
  app_trim_days <- function() {
    if (!isTRUE(input$apply_trim)) return(NULL)
    trim_days <- suppressWarnings(as.integer(input$trim_days))
    if (length(trim_days) != 1 || is.na(trim_days)) trim_days <- 56L
    max(1L, trim_days)
  }

  # Cleaned (pre-QC) uploads are kept so the checks can be re-run when the
  # hemisphere settings change.
  deployment_raw <- reactiveVal(NULL)
  images_raw <- reactiveVal(NULL)

  deployment_checked <- reactiveVal(NULL)
  deployment_issues <- reactiveVal(NULL)
  images_checked <- reactiveVal(NULL)
  images_issues <- reactiveVal(NULL)

  run_deployment_check <- function() {
    req(deployment_raw())
    res <- run_qc_step(
      check_deployments(
        deployment_raw(),
        hemisphere = input$hemisphere,
        lat_hemisphere = input$lat_hemisphere
      )
    )
    deployment_checked(res$value)
    deployment_issues(res[c("messages", "warnings")])
  }

  run_images_check <- function() {
    req(images_raw())
    if (is.null(deployment_checked())) {
      images_checked(NULL)
      images_issues(list(
        messages = character(),
        warnings = "ERROR: ❌ Upload and check a deployment file before the images file."
      ))
      return(invisible(NULL))
    }
    res <- run_qc_step(check_images(images = images_raw(), deployments = deployment_checked()))
    images_checked(res$value)
    images_issues(res[c("messages", "warnings")])
  }

  # Model results on uploaded data are cleared whenever the data the models
  # were fit to change (uploads, hemisphere settings, species, or trim), so
  # the Compare & combine tab never mixes results from different data. A
  # pop-up explains what happened. Prior and other model settings don't
  # change the data, so they don't clear results.
  clear_uploaded_fits <- function(reason) {
    fits <- list(USCR = uscr_nps_fit(), REM = rem_nps_fit(), TTE = tte_nps_fit())
    cleared <- names(fits)[!vapply(fits, is.null, logical(1))]
    if (!length(cleared)) return(invisible(NULL))
    uscr_nps_fit(NULL)
    rem_nps_fit(NULL)
    tte_nps_fit(NULL)
    uscr_nps_debug(initial_nps_debug("USCR"))
    rem_nps_debug(initial_nps_debug("REM"))
    tte_nps_debug(initial_nps_debug("TTE"))
    showModal(modalDialog(
      title = "Model results cleared",
      tags$p(
        "You ", reason, ". The ", paste(cleared, collapse = ", "),
        if (length(cleared) == 1) " result was" else " results were",
        " fit to the previous data, so ",
        if (length(cleared) == 1) "it was" else "they were",
        " cleared. This keeps results from different data from being compared or combined on the ",
        tags$strong("Compare & combine"), " tab."
      ),
      tags$p(
        "Run the models again to get results for the current data. Changing priors or other model settings doesn't clear results; only changes to the data do (uploading files, the hemisphere settings, the species to analyze, or the deployment trim)."
      ),
      easyClose = TRUE,
      footer = modalButton("OK")
    ))
  }

  # A new deployment file clears the images, which were checked against the
  # previous deployment file.
  observeEvent(input$deployment_csv, {
    req(input$deployment_csv)
    raw <- readr::read_csv(input$deployment_csv$datapath, show_col_types = FALSE)
    deployment_raw(clean_deployment_import(raw))
    run_deployment_check()
    clear_uploaded_fits("uploaded a new deployment file")

    if (!is.null(images_raw()) || !is.null(images_checked())) {
      images_raw(NULL)
      images_checked(NULL)
      images_issues(list(
        messages = "A new deployment file was uploaded, so the images file was cleared. Upload the images file again.",
        warnings = character()
      ))
      session$sendCustomMessage("reset_file_input", "images_csv")
    }
  })

  observeEvent(input$images_csv, {
    req(input$images_csv)
    if (is.null(deployment_checked())) {
      showNotification(
        "Please upload and check a deployment CSV before the images CSV.",
        type = "error"
      )
      return(invisible(NULL))
    }
    raw <- readr::read_csv(
      input$images_csv$datapath,
      show_col_types = FALSE,
      col_types = readr::cols(
        `Sighting Count` = readr::col_character(),
        Species = readr::col_character()
      )
    )
    cleaned <- run_qc_step(clean_images_import(raw))
    if (is.null(cleaned$value)) {
      images_raw(NULL)
      images_checked(NULL)
      images_issues(cleaned[c("messages", "warnings")])
      return(invisible(NULL))
    }
    images_raw(cleaned$value)
    run_images_check()
    clear_uploaded_fits("uploaded a new images file")
  })

  # Changing a hemisphere setting changes the Latitude/Longitude sign
  # corrections, so re-run the checks on the uploaded data.
  observeEvent(list(input$hemisphere, input$lat_hemisphere), {
    if (is.null(deployment_raw())) return(invisible(NULL))
    run_deployment_check()
    if (!is.null(images_raw())) run_images_check()
    clear_uploaded_fits("changed a hemisphere setting, which changes the camera coordinates")
  }, ignoreInit = TRUE)

  # Species and trim changes alter the counts and camera-days the models use.
  # The trim is debounced so typing a number of days clears results once.
  observeEvent(input$species_to_analyze, {
    clear_uploaded_fits("changed the species to analyze")
  }, ignoreInit = TRUE)

  trim_setting <- debounce(reactive(list(input$apply_trim, input$trim_days)), 1000)
  observeEvent(trim_setting(), {
    clear_uploaded_fits("changed the deployment trim")
  }, ignoreInit = TRUE)

  output$deployment_check_log <- renderText({
    req(deployment_issues())
    format_check_log(deployment_issues())
  })

  output$deployment_preview <- DT::renderDT({
    req(deployment_checked())
    DT::datatable(deployment_checked(), options = list(pageLength = 10))
  })

  # The analysis window each camera contributes (see analysis_windows()).
  output$deployment_windows <- DT::renderDT({
    req(deployment_checked())
    win <- analysis_windows(deployment_checked(), max_days = app_trim_days())
    # What set each window's start and end; a trim overrides the end.
    flags <- cbind(
      "Start Time missing" = win$start_time_missing,
      "End Time missing" = win$end_time_missing & !win$trimmed,
      "Camera malfunction" = win$malfunction & !win$trimmed,
      "Trimmed" = win$trimmed
    )
    notes <- apply(flags, 1, function(f) paste(colnames(flags)[f], collapse = "; "))
    tbl <- data.frame(
      `Site Name` = win$`Site Name`,
      `Window start` = format(win$window_start, "%m/%d/%Y %H:%M"),
      `Window end` = format(win$window_end, "%m/%d/%Y %H:%M"),
      Days = round(win$days, 2),
      Notes = notes,
      check.names = FALSE
    )
    DT::datatable(tbl, options = list(pageLength = 10), rownames = FALSE)
  })

  output$download_deployment_checked <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$deployment_csv$name), "_CHECKED.csv")
    },
    content = function(file) {
      req(deployment_checked())
      readr::write_csv(deployment_checked(), file)
    }
  )

  output$images_check_log <- renderText({
    req(images_issues())
    format_check_log(images_issues())
  })

  output$images_preview <- DT::renderDT({
    req(images_checked())
    DT::datatable(images_checked(), options = list(pageLength = 10))
  })

  output$download_images_checked <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$images_csv$name), "_CHECKED.csv")
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

  image_summary <- reactive({
    summarize_images_by_species(images_checked())
  })

  available_species <- reactive({
    req(images_checked())
    species <- sort(unique(images_checked()$Species))
    species <- species[!is.na(species) & nzchar(species)]
    species
  })

  observe({
    species <- available_species()
    current <- isolate(input$summary_species)
    current_analysis <- isolate(input$species_to_analyze)

    updateSelectInput(
      session,
      "summary_species",
      choices = species,
      selected = species[1]
    )
    updateSelectInput(
      session,
      "species_to_analyze",
      choices = species
    )
  })

  species_objects <- reactive({
    req(images_checked(), input$summary_species)

    list(
      species_summary = species_summary_per_site(images_checked(), input$summary_species),
      species_counts  = species_counts_per_camera(
        images_checked(),
        input$summary_species,
        deployments = deployment_checked()
      ),
      daily_species   = species_daily_detections(images_checked(), input$summary_species)
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
  # MODEL FIT OBJECTS (reactiveVal) –  NPS PER MODEL
  # ================================================================

  uscr_nps_fit <- reactiveVal(NULL)
  rem_nps_fit  <- reactiveVal(NULL)
  tte_nps_fit  <- reactiveVal(NULL)

  update_model_debug <- function(rv,
                                 status = NULL,
                                 stage = NULL,
                                 started_at = NULL,
                                 finished_at = NULL,
                                 guidance = NULL,
                                 raw_error = NULL,
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

  get_camera_angle_vector <- function(out, default_deg = 55) {
    default_deg <- as.numeric(default_deg[[1]])
    if (!is.finite(default_deg)) default_deg <- 55
    default_deg <- max(1, min(360, default_deg))

    if (!"Camera Detection Angle" %in% names(out)) {
      return(rep(default_deg, nrow(out)))
    }

    theta_raw <- suppressWarnings(as.numeric(out$`Camera Detection Angle`))
    theta_use <- theta_raw
    theta_use[is.na(theta_use) | theta_use <= 0 | theta_use > 360] <- default_deg
    theta_use
  }

  summarize_camera_angle_source <- function(out, default_deg = 55) {
    theta_use <- get_camera_angle_vector(out, default_deg = default_deg)
    has_column <- "Camera Detection Angle" %in% names(out)
    theta_raw <- if (has_column) suppressWarnings(as.numeric(out$`Camera Detection Angle`)) else numeric(0)
    n_default <- if (has_column) {
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
      return(paste("Camera detection angle:", range_text, "(all cameras using default from Model settings)"))
    }

    if (n_default > 0) {
      return(paste(
        "Camera detection angle:",
        range_text,
        paste0("(", n_default, " camera(s) using the default angle from Model settings)")
      ))
    }

    paste("Camera detection angle:", range_text, "(from uploaded deployment data)")
  }

  summarize_uscr_context <- function(d, source_label = "Uploaded") {
    effective_iter <- as.integer(input$iter_uscr)
    effective_burnin <- as.integer(input$burnin_uscr)
    effective_thin <- as.integer(input$thin_uscr)
    effective_M <- as.integer(input$M_uscr)

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
      paste("USCR state-space buffer:", format_num(app_uscr_buffer_m(), 0), "m"),
      paste("USCR density prior:", describe_density_prior_spec(density_prior_spec())),
      paste("Camera-level noise:", noise_label(app_camera_noise())),
      "USCR reruns automatically at the same chain count shown above: iterations and burn-in double when Rhat is too high, M doubles when the posterior for abundance (N) presses against M, and the state-space buffer grows when the posterior sigma implies animals near the edge could still be detected.",
      if (source_label == "Uploaded") {
        "Supported sources here: uploaded field data."
      } else {
        "Supported sources here: uploaded field data."
      },
      "Preprocessing: total animal events per camera, camera-days, and buffered spatial state space."
    )

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
      summarize_camera_angle_source(d$out, default_deg = input$theta),
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
      paste("Density prior:", describe_density_prior_spec(density_prior_spec())),
      paste("Camera-level noise:", noise_label(app_camera_noise())),
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
      summarize_camera_angle_source(d$out, default_deg = input$theta),
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
      paste("Density prior:", describe_density_prior_spec(density_prior_spec())),
      paste("Camera-level noise:", noise_label(app_camera_noise())),
      "Preprocessing: TTE uses total animal events per camera and camera-days."
    )
  }

  summarize_shared_sim_context <- function(d, model_label) {
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
          ", default angle=", input$theta,
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
          ", default angle=", input$theta,
          " degrees"
        )
      ),
      NULL
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
      tips <- c(tips, "USCR doubles iterations and burn-in automatically when Rhat is too high, doubles M when the posterior for abundance presses against it, and grows the state-space buffer when the posterior sigma is larger than the buffer safely covers, so longer runs can be expected when any check fails. If it is too slow for testing, reduce chains, iterations, or M.")
    }
    if (identical(model, "REM")) {
      tips <- c(tips, "REM expects per-camera animal events, camera-days, and detection distances.")
      tips <- c(tips, "REM doubles iterations and burn-in automatically when Rhat is too high, and doubles D_max when the posterior for density presses against it, so longer runs can be expected when either check fails.")
    }
    if (identical(model, "TTE")) {
      tips <- c(tips, "TTE expects per-camera animal events, camera-days, and detection distances.")
      tips <- c(tips, "TTE doubles iterations and burn-in automatically when Rhat is too high, and doubles D_max when the posterior for density presses against it, so longer runs can be expected when either check fails.")
    }

    paste(tips, collapse = "\n")
  }

  # Single shared status callback used by all six run paths (USCR/REM/TTE x
  # simulated/uploaded). The model engine calls this once per MCMC round with
  # the round number and this round's settings; once a second round happens
  # (Rhat didn't hit the target, so the model reruns with more iterations),
  # it also gets the previous round's elapsed time so we can show a rough ETA.
  # USCR also passes M (its augmentation bound) and buffer_m (state-space
  # buffer); REM/TTE don't have either.
  make_model_status_callback <- function(rv, model_label, notification_id) {
    function(round, iter, burnin, M = NULL, buffer_m = NULL, prev_elapsed_min = NULL, prev_iter = NULL) {
      message <- paste0(model_label, ": ")
      if (!is.null(M)) message <- paste0(message, "M=", M, ", ")
      if (!is.null(buffer_m)) message <- paste0(message, "buffer=", round(buffer_m), "m, ")
      message <- paste0(message, iter, " iterations, ", burnin, " burn-in")
      if (round > 1L && !is.null(prev_elapsed_min) && !is.null(prev_iter) && prev_iter > 0) {
        est_min <- prev_elapsed_min * (iter / prev_iter)
        message <- paste0(message, " — approx. ", format(round(est_min, 1), nsmall = 1), " min")
      }
      showNotification(
        message,
        type = "message",
        duration = NULL,
        id = notification_id
      )
      update_model_debug(rv, stage = paste("Round", round), log_entry = message)
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
      if (!is.null(fit$settings$buffer_m)) {
        setting_lines <- c(setting_lines, paste0("- Final buffer: ", round(fit$settings$buffer_m), " m"))
      }
      if (!is.null(fit$settings$D_max)) {
        setting_lines <- c(setting_lines, paste0("- Final D_max: ", format_num(fit$settings$D_max, 1), " animals/km\u00b2"))
      }
      if (!is.null(fit$density_prior)) {
        prior_text <- if (!is.null(fit$density_prior$M)) {
          describe_uscr_density_prior(fit$density_prior)
        } else {
          describe_beta_density_prior(fit$density_prior)
        }
        setting_lines <- c(setting_lines, paste0("- Density prior: ", prior_text))
        if (isTRUE(fit$density_prior$off_target)) {
          setting_lines <- c(setting_lines, "  (a bound of the prior's 95% interval differs from the target by more than 10% of the interval's width)")
        }
      }
      if (!is.null(fit$settings$camera_noise)) {
        setting_lines <- c(setting_lines, paste0("- Camera-level noise: ", noise_label(fit$settings$camera_noise)))
      }
      lines <- c(lines, "", "Run settings used:", setting_lines)
      if (!is.null(fit$final_rhat_max) && is.finite(fit$final_rhat_max)) {
        lines <- c(lines, paste0("- Final max Rhat: ", format_num(fit$final_rhat_max, 3)))
      }
      ppc_line <- ppc_text(fit)
      if (!is.null(ppc_line)) {
        lines <- c(lines, "", ppc_line)
      }
    }

    if (!is.null(fit) && !is.null(fit$round_history) && nrow(fit$round_history) > 0) {
      round_df <- fit$round_history
      round_lines <- apply(
        round_df,
        1,
        function(row) {
          parts <- c(
            paste0("round ", row[["round"]]),
            paste0("iter=", row[["niter"]]),
            paste0("thin=", row[["thin"]]),
            paste0("chains=", row[["n_chains"]]),
            paste0("rhat=", format_num(as.numeric(row[["rhat_max"]]), 3))
          )
          if ("M" %in% names(round_df)) {
            parts <- append(parts, paste0("M=", row[["M"]]), after = 1)
          }
          if ("buffer_m" %in% names(round_df)) {
            parts <- append(parts, paste0("buffer=", round(as.numeric(row[["buffer_m"]])), "m"), after = 2)
          }
          if ("D_max" %in% names(round_df)) {
            parts <- append(parts, paste0("D_max=", format_num(as.numeric(row[["D_max"]]), 1)), after = 1)
          }
          if ("D_too_small" %in% names(round_df)) {
            parts <- c(parts, paste0("D_max_too_small=", row[["D_too_small"]] %in% c("TRUE", "1")))
          }
          if ("M_too_small" %in% names(round_df)) {
            parts <- c(parts, paste0("M_too_small=", row[["M_too_small"]] %in% c("TRUE", "1")))
          }
          if ("buffer_too_small" %in% names(round_df)) {
            parts <- c(parts, paste0("buffer_too_small=", row[["buffer_too_small"]] %in% c("TRUE", "1")))
          }
          paste(parts, collapse = ", ")
        }
      )
      lines <- c(lines, "", paste0(state$model, " run history:"), paste0("- ", round_lines))
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

  # Starting run-status panel for each model on uploaded data.
  initial_nps_debug <- function(model) {
    switch(
      model,
      USCR = make_model_debug(
        "USCR", "Uploaded field data", "Uploaded field data; simulated grid",
        "Total animal events per camera, camera-days, and buffered state space."
      ),
      REM = make_model_debug(
        "REM", "Uploaded field data", "Uploaded field data only",
        "REM uses total animal events per camera and camera-days."
      ),
      TTE = make_model_debug(
        "TTE", "Uploaded field data", "Uploaded field data only",
        "TTE uses total animal events per camera and camera-days."
      )
    )
  }
  uscr_nps_debug <- reactiveVal(initial_nps_debug("USCR"))
  rem_nps_debug <- reactiveVal(initial_nps_debug("REM"))
  tte_nps_debug <- reactiveVal(initial_nps_debug("TTE"))

  # Status tracking for model runs
  uscr_nps_running <- reactiveVal(FALSE)
  rem_nps_running <- reactiveVal(FALSE)
  tte_nps_running <- reactiveVal(FALSE)

  nps_model_inputs <- reactive({
    req(deployment_checked(), images_checked(), !is.null(input$species_to_analyze))
    build_nps_model_inputs(
      deployment_checked(),
      images_checked(),
      input$species_to_analyze,
      max_days = app_trim_days()
    )
  })

  app_n_chains <- function() {
    max(2L, as.integer(input$n_chains))
  }

  # ---- Effective prior settings ----
  # Simple mode always uses the default priors; Advanced mode uses the inputs.
  is_simple_mode <- function() !identical(input$mode, "Advanced")
  prior_input <- function(id) {
    if (is_simple_mode()) advanced_prior_defaults[[id]] else input[[id]]
  }
  default_uscr_buffer_m <- function() {
    ceiling(uscr_default_buffer_m(advanced_prior_defaults$log_sigma_mean,
                                  advanced_prior_defaults$log_sigma_sd) / 50) * 50
  }
  app_uscr_buffer_m <- function() {
    if (is_simple_mode()) default_uscr_buffer_m() else input$uscr_buffer_m
  }

  # Switching to Simple resets every Advanced input to its default.
  observeEvent(input$mode, {
    if (!is_simple_mode()) return()
    for (id in names(advanced_prior_defaults)) {
      updateNumericInput(session, id, value = advanced_prior_defaults[[id]])
    }
    updateNumericInput(session, "uscr_buffer_m", value = default_uscr_buffer_m())
    updateCheckboxInput(session, "D_prior_uniform", value = TRUE)
    updateCheckboxInput(session, "camera_noise", value = FALSE)
    for (id in c("D_prior_lo", "D_prior_hi")) {
      updateNumericInput(session, id, value = NA)
    }
  }, ignoreInit = TRUE)

  # Density prior (all models): uniform in Simple mode or when the box is checked;
  # otherwise the user's mean and 95% interval (animals/km^2).
  density_prior_spec <- function() {
    if (is_simple_mode() || !identical(input$D_prior_uniform, FALSE)) {
      return(list(uniform = TRUE))
    }
    list(
      uniform = FALSE,
      lo = input$D_prior_lo,
      hi = input$D_prior_hi
    )
  }
  density_prior_problem <- function(spec = density_prior_spec()) {
    if (isTRUE(spec$uniform)) return(NULL)
    check_density_prior(spec$lo, spec$hi)
  }
  describe_density_prior_spec <- function(spec) {
    if (isTRUE(spec$uniform)) return("uniform (default)")
    paste0("informative, 95% interval ", spec$lo, " to ", spec$hi, " animals/km\u00b2")
  }

  # Starting M for shared-design USCR runs: about twice the abundance expected
  # under the true density (the density prior can raise it further).
  pshared_suggested_M <- function(s) {
    area <- uscr_buffered_area_km2(
      cbind(s$uscr$data$out$utm_e, s$uscr$data$out$utm_n),
      s$settings$buffer_m / 1000
    )
    max(50, ceiling(2 * s$D_true * area / 50) * 50)
  }

  # ---- Density prior preview (Model settings > Advanced > Priors: USCR) ----
  # The prior on density depends on the state-space area (and M), so the
  # preview uses uploaded data if available, otherwise the shared design.
  uscr_prior_preview_basis <- reactive({
    d <- tryCatch(nps_model_inputs(), error = function(e) NULL)
    if (!is.null(d) && nrow(d$out) >= 2) {
      return(list(
        label = "your uploaded field data",
        out = d$out,
        buffer_m = app_uscr_buffer_m(),
        M_start = as.integer(input$M_uscr)
      ))
    }
    s <- tryCatch(pshared_sim(), error = function(e) NULL)
    if (!is.null(s)) {
      return(list(
        label = "the current shared-design simulation",
        out = s$uscr$data$out,
        buffer_m = s$settings$buffer_m,
        M_start = pshared_suggested_M(s)
      ))
    }
    NULL
  })

  # Debounced so the preview waits until typing pauses (solving the prior for
  # a large M takes up to a second or so).
  density_prior_spec_debounced <- debounce(reactive(density_prior_spec()), 600)

  density_prior_preview <- reactive({
    spec <- density_prior_spec_debounced()
    if (isTRUE(spec$uniform)) return(NULL)
    problem <- density_prior_problem(spec)
    if (!is.null(problem)) return(list(message = problem, kind = "error"))

    beta <- tryCatch(
      design_beta_density_prior(spec$lo, spec$hi, encounter_prior_start_D_max(spec$hi)),
      error = function(e) e
    )
    if (inherits(beta, "error")) {
      return(list(message = conditionMessage(beta), kind = "error"))
    }
    out <- list(beta = beta)

    # USCR's prior depends on the state-space area, so it needs camera data.
    basis <- uscr_prior_preview_basis()
    if (is.null(basis)) {
      out$uscr_message <- "its prior depends on the size of the state space, so it appears here once you upload field data or simulate a shared design."
      return(out)
    }
    area <- uscr_state_space_and_area(basis$out, buffer_m = basis$buffer_m)$area_km2
    M <- max(basis$M_start, uscr_prior_start_M(spec$hi, area))
    uscr <- tryCatch(design_uscr_density_prior(spec$lo, spec$hi, area, M), error = function(e) e)
    if (inherits(uscr, "error")) {
      out$uscr_message <- conditionMessage(uscr)
    } else {
      out$uscr <- uscr
      out$basis <- basis
    }
    out
  })

  # Note: pv[["uscr"]] rather than pv$uscr, because $ partially matches
  # "uscr_message" when there is no USCR prior yet.
  output$D_prior_readout <- renderUI({
    pv <- density_prior_preview()
    if (is.null(pv)) return(NULL)
    if (!is.null(pv$message)) {
      color <- if (identical(pv$kind, "error")) "firebrick" else "var(--muted)"
      return(tags$p(class = "small", style = paste0("color: ", color, ";"), pv$message))
    }
    fmt <- function(x) format(signif(x, 3), big.mark = ",")
    prior_line <- function(label, d, note) {
      tags$p(
        class = "small",
        tags$strong(label),
        "95% interval ", fmt(d$achieved[["lo"]]), " to ", fmt(d$achieved[["hi"]]),
        " animals/km\u00b2, mean ", fmt(d$achieved[["mean"]]),
        ", median ", fmt(d$achieved[["median"]]), ". ", note
      )
    }
    off_target_note <- function(label, d, why) {
      if (!isTRUE(d$off_target)) return(NULL)
      tags$p(
        class = "small",
        style = "color: firebrick;",
        tags$strong(paste0(label, ": interval not matched. ")),
        "At least one bound differs from your interval by more than ",
        round(100 * d$miss_tolerance), "% of its width. ", why,
        "Widen the interval, or use the prior shown in the plot if it is acceptable."
      )
    }
    tagList(
      prior_line("REM & TTE: ", pv$beta, paste0("(Starting D_max = ", fmt(pv$beta$D_max), " animals/km\u00b2.)")),
      off_target_note("REM & TTE", pv$beta, ""),
      if (!is.null(pv[["uscr"]])) {
        prior_line(
          "USCR: ", pv[["uscr"]],
          paste0("(Based on ", pv$basis$label, ": state space ", fmt(pv[["uscr"]]$area_km2),
                 " km\u00b2, starting M = ", pv[["uscr"]]$M, ".)")
        )
      } else {
        tags$p(class = "small", style = "color: var(--muted);", tags$strong("USCR: "), pv[["uscr_message"]])
      },
      if (!is.null(pv[["uscr"]])) {
        off_target_note(
          "USCR", pv[["uscr"]],
          "This usually means the interval is very narrow: abundance in the state space is a whole-number count with its own randomness, which limits how tight the prior on density can be. "
        )
      },
      tags$p(
        class = "small",
        style = "color: var(--muted);",
        "During fitting, the app may enlarge D_max, M, or the USCR state space. When it does, it re-solves the prior to keep your interval, so the prior used in the fit can differ slightly from this preview. The run-status panel reports the prior actually used."
      )
    )
  })

  output$D_prior_plot <- renderPlot({
    pv <- density_prior_preview()
    req(pv, pv$beta)
    plot_density_priors(pv$beta, pv[["uscr"]])
  }, res = 110)

  # Camera-level noise (Poisson-lognormal): off in Simple mode and by default.
  app_camera_noise <- function() {
    !is_simple_mode() && isTRUE(input$camera_noise)
  }
  noise_label <- function(on) if (isTRUE(on)) "on (Poisson-lognormal counts)" else "off (Poisson counts)"

  # Posterior predictive check text for a fit (see ppc_summary()).
  ppc_text <- function(fit) {
    ppc <- ppc_summary(fit)
    if (is.null(ppc)) return(NULL)
    base <- sprintf("Posterior predictive check: Bayesian p-value = %.2f", ppc$p)
    if (!ppc$fail) {
      return(paste0(base, " (passes; values from ", ppc$lower, " to ", ppc$upper, " pass)."))
    }
    paste0(
      base, ". FAILED (outside ", ppc$lower, " to ", ppc$upper, "): the observed counts vary ",
      if (ppc$p < ppc$lower) "more" else "less",
      " than the model expects, so its estimates and intervals may be unreliable.",
      if (!isTRUE(fit$settings$camera_noise)) {
        " Turning on camera-level noise (Model settings > Advanced > Priors: Camera-level noise) can often fix a failed check."
      } else {
        ""
      }
    )
  }
  notify_ppc <- function(fit, label) {
    ppc <- ppc_summary(fit)
    if (!is.null(ppc) && isTRUE(ppc$fail)) {
      showNotification(
        paste0(label, " failed its posterior predictive check (Bayesian p-value ",
               sprintf("%.2f", ppc$p), "). See the run-status panel."),
        type = "warning", duration = NULL
      )
    }
  }

  # ---- "Run with" notes next to each Run button ----
  # Lists what differs from the app defaults for the model that button runs.
  run_settings_diffs <- function(model, shared = FALSE) {
    fmt <- function(x) format(signif(as.numeric(x), 4))
    changed <- function(x, default, tol = 1e-8) {
      is.numeric(x) && length(x) == 1 && is.finite(x) && abs(x - default) > tol
    }
    # Priors centered on simulated values differ from the default means by
    # rounding at the default slider settings (e.g. log(3.09) vs 1.130), so
    # compare those on a coarser tolerance.
    sim_tol <- 0.01
    diffs <- character(0)
    add_prior <- function(id, label) {
      val <- prior_input(id)
      if (changed(val, advanced_prior_defaults[[id]])) {
        diffs <<- c(diffs, paste0(label, ": ", fmt(val), " (default ", fmt(advanced_prior_defaults[[id]]), ")"))
      }
    }
    s <- if (shared) tryCatch(pshared_sim(), error = function(e) NULL) else NULL

    spec <- density_prior_spec()
    if (!isTRUE(spec$uniform)) {
      diffs <- c(diffs, paste0("Density prior: 95% interval ", fmt(spec$lo), "\u2013", fmt(spec$hi), " animals/km\u00b2"))
    }

    noise_on <- if (shared) isTRUE(s$settings$camera_noise) else app_camera_noise()
    if (noise_on) {
      diffs <- c(diffs, paste0(
        "Camera-level noise: on, \u03c3<sub>\u03b5</sub> ~ Gamma(", fmt(prior_input("sd_eps_shape")), ", ",
        fmt(prior_input("sd_eps_rate")), ")", if (shared) " (from the simulated design)" else ""
      ))
    }

    if (identical(model, "USCR")) {
      if (shared && !is.null(s)) {
        sigma_km <- s$uscr$sim$truth$sigma_m / 1000
        if (changed(log(sigma_km), advanced_prior_defaults$log_sigma_mean, sim_tol)) {
          diffs <- c(diffs, paste0("log(\u03c3) mean: ", fmt(log(sigma_km)), " (centered on the simulated \u03c3; default ",
                                   fmt(advanced_prior_defaults$log_sigma_mean), ")"))
        }
        if (changed(log(s$settings$lambda0), advanced_prior_defaults$log_lam0_mean, sim_tol)) {
          diffs <- c(diffs, paste0("log(\u03bb\u2080) mean: ", fmt(log(s$settings$lambda0)),
                                   " (centered on the simulated \u03bb\u2080; default ", fmt(advanced_prior_defaults$log_lam0_mean), ")"))
        }
      } else if (!shared) {
        add_prior("log_sigma_mean", "log(\u03c3) mean")
        add_prior("log_lam0_mean", "log(\u03bb\u2080) mean")
      }
      add_prior("log_sigma_sd", "log(\u03c3) SD")
      add_prior("log_lam0_sd", "log(\u03bb\u2080) SD")
      # Buffer: flag only a user override of the automatic default.
      if (shared && !is.null(s)) {
        auto <- s$settings$buffer_default_m
        if (!is.null(auto) && changed(s$settings$buffer_m, auto)) {
          diffs <- c(diffs, paste0("State-space buffer: ", fmt(s$settings$buffer_m), " m (default ", fmt(auto), " m)"))
        }
      } else if (!shared) {
        auto <- ceiling(uscr_default_buffer_m(prior_input("log_sigma_mean"), prior_input("log_sigma_sd")) / 50) * 50
        if (changed(app_uscr_buffer_m(), auto)) {
          diffs <- c(diffs, paste0("State-space buffer: ", fmt(app_uscr_buffer_m()), " m (default ", fmt(auto), " m)"))
        }
      }
    } else {
      if (shared && !is.null(s)) {
        if (changed(log(s$settings$v), advanced_prior_defaults$log_v_mean, sim_tol)) {
          diffs <- c(diffs, paste0("log(v) mean: ", fmt(log(s$settings$v)), " (centered on the simulated speed; default ",
                                   fmt(advanced_prior_defaults$log_v_mean), ")"))
        }
      } else if (!shared) {
        add_prior("log_v_mean", "log(v) mean")
      }
      add_prior("log_v_sd", "log(v) SD")
    }

    if (app_max_adapt_rounds() != 6L) {
      diffs <- c(diffs, paste0("Maximum adaptive rounds: ", app_max_adapt_rounds(), " (default 6)"))
    }
    diffs
  }

  render_run_note <- function(model, shared = FALSE) {
    renderUI({
      if (shared) req(pshared_sim())
      diffs <- run_settings_diffs(model, shared)
      tags$div(
        class = "small",
        style = "margin: 0.4rem 0 0.5rem 0; color: var(--muted);",
        if (!length(diffs)) {
          tags$strong("Run with all default settings.")
        } else {
          tagList(
            tags$strong("The following differ from default:"),
            # HTML() so the <sub> in sigma_epsilon renders; diffs hold only app text and numbers.
            tags$ul(style = "margin: 0.2rem 0 0 0;", lapply(diffs, function(d) tags$li(HTML(d))))
          )
        }
      )
    })
  }
  output$run_note_uscr_nps     <- render_run_note("USCR")
  output$run_note_rem_nps      <- render_run_note("REM")
  output$run_note_tte_nps      <- render_run_note("TTE")
  output$run_note_uscr_pshared <- render_run_note("USCR", shared = TRUE)
  output$run_note_rem_pshared  <- render_run_note("REM", shared = TRUE)
  output$run_note_tte_pshared  <- render_run_note("TTE", shared = TRUE)

  app_max_adapt_rounds <- function() {
    x <- suppressWarnings(as.integer(input$max_adapt_rounds))
    if (length(x) != 1L || is.na(x) || x < 1L) 6L else x
  }

  # Run-status guidance for a completed fit. Flags fits that stopped at the
  # adaptive-round limit without passing every check.
  adaptive_limit_guidance <- function(fit, model, default_guidance) {
    if (!isTRUE(fit$hit_max_adapt_rounds)) return(default_guidance)
    paste0(
      model, " stopped after ", fit$settings$max_adapt_rounds,
      " adaptive rounds without passing these checks: ",
      paste(fit$unmet_checks, collapse = ", "),
      ". The results shown are from the last round and may not be reliable. ",
      "To keep going, increase 'Maximum adaptive rounds' in Model settings and run ",
      model, " again.",
      if (any(c("M", "D_max") %in% fit$unmet_checks) && isTRUE(fit$density_prior$uniform)) {
        paste0(
          " The posterior kept pressing against the upper limit on ",
          if ("M" %in% fit$unmet_checks) "abundance (M)" else "density (D_max)",
          " under the uniform density prior. An informative density prior ",
          "(Model settings > Advanced > Priors: Density (all models)) can rule out ",
          "implausibly high densities and stop the limit from growing."
        )
      }
    )
  }

  notify_adaptive_limit <- function(fit, label) {
    if (isTRUE(fit$hit_max_adapt_rounds)) {
      showNotification(
        paste0(label, " stopped at the adaptive-round limit without passing all checks. See the run-status panel."),
        type = "warning", duration = NULL
      )
    }
  }

  # Keep the USCR buffer default in sync with the log(sigma) prior: the
  # buffer should comfortably contain a 99% home range under that prior, so
  # if the prior mean/SD change, recompute and refresh the displayed default.
  observe({
    req(is.numeric(input$log_sigma_mean), is.numeric(input$log_sigma_sd))
    new_default <- uscr_default_buffer_m(input$log_sigma_mean, input$log_sigma_sd)
    updateNumericInput(session, "uscr_buffer_m", value = ceiling(new_default / 50) * 50)
  })

  uscr_run_args <- function(waic = TRUE, status_callback = NULL,
                            buffer_m = NULL, log_lam0_mean = NULL,
                            log_sigma_mean = NULL) {
    args <- list(
      iter = input$iter_uscr,
      burnin = input$burnin_uscr,
      thin = input$thin_uscr,
      n_chains = app_n_chains(),
      M = input$M_uscr,
      log_sigma_mean = log_sigma_mean %||% prior_input("log_sigma_mean"),
      log_sigma_sd = prior_input("log_sigma_sd"),
      buffer_m = buffer_m %||% app_uscr_buffer_m(),
      log_lam0_mean = log_lam0_mean %||% prior_input("log_lam0_mean"),
      log_lam0_sd =  prior_input("log_lam0_sd"),
      sd_eps_shape = prior_input("sd_eps_shape"),
      sd_eps_rate = prior_input("sd_eps_rate"),
      adaptive = TRUE,
      max_adapt_rounds = app_max_adapt_rounds(),
      density_prior = density_prior_spec(),
      camera_noise = app_camera_noise(),
      compute_WAIC = waic,
      diagnostic_mode = FALSE,
      parallel_chains = isTRUE(app_n_chains() > 1),
      status_callback = status_callback,
      verbose = FALSE
    )

    # Keep the app compatible with whichever run_USCR() definition is
    # currently loaded, including older variants without these extras.
    valid_names <- names(formals(run_USCR_app))
    args[names(args) %in% valid_names]
  }


    # --- USCR: NPS ---

  start_uscr_nps_run <- function() {
    req(nps_model_inputs())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("USCR not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    d <- nps_model_inputs()
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
      guidance = "USCR supports uploaded field data. Watch this panel for round-by-round Rhat, M, and buffer checks.",
      raw_error = NULL,
      context = summarize_uscr_context(d, source_label = "Uploaded"),
      log_entry = "Uploaded-data USCR run requested."
    )
    validate(
      need(nrow(d$out) >= 2,
           "At least 2 cameras are required to estimate density; only 1 camera was found in this dataset."),
      need(all(d$camera_days > 0),
           "Some cameras have zero or negative camera-days; check deployment dates, or check for overlapping/duplicate Cluster IDs at the same site (overlapping detection events can double-count time and drive camera-days to zero or below).")
    )

    uscr_nps_running(TRUE)
    uscr_nps_fit(NULL)  # Clear previous results

    showNotification(
      "Running USCR on uploaded data.",
      type = "message",
      duration = NULL,
      id = "uscr_nps_status"
    )

    fit <- tryCatch(
      {
        update_model_debug(
          uscr_nps_debug,
          stage = "Preparing run",
          log_entry = "Input checks passed. Starting USCR setup."
        )

        status_callback <- make_model_status_callback(uscr_nps_debug, "USCR", "uscr_nps_status")

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
              status_callback = status_callback
            )
          )
        )
      },
      error = function(e) {
        removeNotification("uscr_nps_status")
        update_model_debug(
          uscr_nps_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("USCR", "uploaded field data", e$message),
          raw_error = e$message,
          log_entry = paste("Uploaded-data USCR failed:", e$message)
        )
        showNotification(
          paste("USCR (uploaded data) failed:", e$message),
          type = "error", duration = NULL
        )
        return(NULL)
      },
      finally = {
        removeNotification("uscr_nps_status")
        uscr_nps_running(FALSE)
      }
    )
    uscr_nps_fit(fit)
    if (!is.null(fit)) {
      removeNotification("uscr_nps_status")
      update_model_debug(
        uscr_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = adaptive_limit_guidance(fit, "USCR", "Uploaded-data USCR completed successfully. Review the summary above and the run history below."),
        log_entry = "Uploaded-data USCR run completed."
      )
      showNotification("USCR (uploaded data) complete!", type = "message")
      notify_adaptive_limit(fit, "USCR (uploaded data)")
      notify_ppc(fit, "USCR (uploaded data)")
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
    req(nps_model_inputs())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("REM not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    d <- nps_model_inputs()
    rem_nps_debug(make_model_debug(
      "REM",
      "Uploaded field data",
      "Uploaded field data only",
      "REM uses total animal events per camera and camera-days."
    ))
    update_model_debug(
      rem_nps_debug,
      status = "running",
      stage = "Preflight checks",
      started_at = Sys.time(),
      guidance = "REM only runs on uploaded field data in this app.",
      raw_error = NULL,
      context = summarize_rem_context(d),
      log_entry = "Uploaded-data REM run requested."
    )
    validate(
      need(nrow(d$out) >= 2,
           "At least 2 cameras are required to estimate density; only 1 camera was found in this dataset."),
      need(all(d$camera_days > 0),
           "Some cameras have zero or negative camera-days; check deployment dates, or check for overlapping/duplicate Cluster IDs at the same site (overlapping detection events can double-count time and drive camera-days to zero or below).")
    )

    rem_nps_running(TRUE)
    rem_nps_fit(NULL)  # Clear previous results

    showNotification(
      "Running REM on uploaded data... This may take several minutes.",
      type = "message",
      duration = NULL,
      id = "rem_nps_status"
    )

    fit <- tryCatch(
      {
        status_callback <- make_model_status_callback(rem_nps_debug, "REM", "rem_nps_status")

        run_REM(
          y            = d$camera_counts,
          r_km         = d$out$`Detection Distance` / 1000,
          camera_days  = d$camera_days,
          theta_deg    = get_camera_angle_vector(d$out, default_deg = input$theta),
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = app_n_chains(),
          D_max        = rem_tte_start_D_max,
          log_v_mean   = prior_input("log_v_mean"),
          log_v_sd     = prior_input("log_v_sd"),
          sd_eps_shape = prior_input("sd_eps_shape"),
          sd_eps_rate  = prior_input("sd_eps_rate"),
          max_adapt_rounds = app_max_adapt_rounds(),
          density_prior = density_prior_spec(),
          camera_noise = app_camera_noise(),
          status_callback = status_callback
        )
      },
      error = function(e) {
        removeNotification("rem_nps_status")
        update_model_debug(
          rem_nps_debug,
          status = "error",
          stage = "Failed",
          finished_at = Sys.time(),
          guidance = friendly_model_error("REM", "uploaded field data", e$message),
          raw_error = e$message,
          log_entry = paste("Uploaded-data REM failed:", e$message)
        )
        showNotification(
          paste("REM (uploaded data) failed:", e$message),
          type = "error", duration = NULL
        )
        return(NULL)
      },
      finally = {
        removeNotification("rem_nps_status")
        rem_nps_running(FALSE)
      }
    )
    rem_nps_fit(fit)
    if (!is.null(fit)) {
      removeNotification("rem_nps_status")
      update_model_debug(
        rem_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = adaptive_limit_guidance(fit, "REM", "REM completed successfully. Review the summary above and the run history below."),
        log_entry = "Uploaded-data REM run completed."
      )
      showNotification("REM (uploaded data) complete!", type = "message")
      notify_adaptive_limit(fit, "REM (uploaded data)")
      notify_ppc(fit, "REM (uploaded data)")
    }
  })

  # --- TTE: NPS ---

  observeEvent(input$run_tte_nps, {
    req(nps_model_inputs())
    prior_problem <- density_prior_problem()
    if (!is.null(prior_problem)) {
      showNotification(paste("TTE not run. Density prior:", prior_problem), type = "error", duration = NULL)
      return(invisible(NULL))
    }
    d <- nps_model_inputs()
    tte_nps_debug(make_model_debug(
      "TTE",
      "Uploaded field data",
      "Uploaded field data only",
      "TTE uses total animal events per camera and camera-days."
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
      need(nrow(d$out) >= 2,
           "At least 2 cameras are required to estimate density; only 1 camera was found in this dataset."),
      need(all(d$camera_days > 0),
           "Some cameras have zero or negative camera-days; check deployment dates, or check for overlapping/duplicate Cluster IDs at the same site (overlapping detection events can double-count time and drive camera-days to zero or below).")
    )

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
        status_callback <- make_model_status_callback(tte_nps_debug, "TTE", "tte_nps_status")

        run_TTE(
          y            = d$camera_counts,
          r_km         = d$out$`Detection Distance` / 1000,
          camera_days  = d$camera_days,
          theta_deg    = get_camera_angle_vector(d$out, default_deg = input$theta),
          iter         = input$iter_rem_tte,
          burnin       = input$burnin_rem_tte,
          thin         = input$thin_rem_tte,
          n_chains     = app_n_chains(),
          D_max        = rem_tte_start_D_max,
          log_v_mean   = prior_input("log_v_mean"),
          log_v_sd     = prior_input("log_v_sd"),
          sd_eps_shape = prior_input("sd_eps_shape"),
          sd_eps_rate  = prior_input("sd_eps_rate"),
          max_adapt_rounds = app_max_adapt_rounds(),
          density_prior = density_prior_spec(),
          camera_noise = app_camera_noise(),
          status_callback = status_callback
        )
      },
      error = function(e) {
        removeNotification("tte_nps_status")
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
        return(NULL)
      },
      finally = {
        removeNotification("tte_nps_status")
        tte_nps_running(FALSE)
      }
    )
    tte_nps_fit(fit)
    if (!is.null(fit)) {
      removeNotification("tte_nps_status")
      update_model_debug(
        tte_nps_debug,
        status = "success",
        stage = "Complete",
        finished_at = Sys.time(),
        guidance = adaptive_limit_guidance(fit, "TTE", "TTE completed successfully. Review the summary above and the run history below."),
        log_entry = "Uploaded-data TTE run completed."
      )
      showNotification("TTE (uploaded data) complete!", type = "message")
      notify_adaptive_limit(fit, "TTE (uploaded data)")
      notify_ppc(fit, "TTE (uploaded data)")
    }
  })

  # ================================================================
  # MODEL TAB OUTPUTS (USCR / REM / TTE)
  # ================================================================

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
    fit <- uscr_nps_fit()
    if (is.null(fit)) {
      cat("USCR (uploaded data): not run yet. Click 'Run USCR on uploaded data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Uploaded field data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      median_density_animals_per_km2 = round(s$median_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    ppc <- ppc_summary(fit)
    if (!is.null(ppc)) {
      out$bayesian_p_value <- round(ppc$p, 2)
      out$posterior_predictive_check <- if (ppc$fail) "FAILED (see run status below)" else "passed"
    }
    out
  })

  output$rem_nps_text <- renderPrint({
    dbg <- rem_nps_debug()
    if (rem_nps_running()) {
      cat("⏳ REM model is running...\n")
      cat("Current stage:", dbg$stage, "\n")
      cat("Open 'Run status & troubleshooting' below for more detail.\n")
      return(invisible(NULL))
    }
    if (identical(dbg$status, "error")) {
      cat("REM (uploaded data): the last run failed.\n")
      cat("See 'Run status & troubleshooting' below for the raw error and guidance.\n")
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
      median_density_animals_per_km2 = round(s$median_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    ppc <- ppc_summary(fit)
    if (!is.null(ppc)) {
      out$bayesian_p_value <- round(ppc$p, 2)
      out$posterior_predictive_check <- if (ppc$fail) "FAILED (see run status below)" else "passed"
    }
    out
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
    fit <- tte_nps_fit()
    if (is.null(fit)) {
      cat("TTE (uploaded data): not run yet. Click 'Run TTE on uploaded data'.")
      return(invisible(NULL))
    }
    s <- summarize_method(fit)
    out <- list(
      dataset                   = "Uploaded field data",
      mean_density_animals_per_km2 = round(s$mean_km2, 2),
      median_density_animals_per_km2 = round(s$median_km2, 2),
      CI95_km2                  = c(round(s$q2.5_km2, 2), round(s$q97.5_km2, 2))
    )
    if (is.finite(s$waic)) out$WAIC <- round(s$waic, 2)
    ppc <- ppc_summary(fit)
    if (!is.null(ppc)) {
      out$bayesian_p_value <- round(ppc$p, 2)
      out$posterior_predictive_check <- if (ppc$fail) "FAILED (see run status below)" else "passed"
    }
    out
  })

   output$uscr_nps_debug <- renderText({
    state <- uscr_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, uscr_nps_fit())
  })


  output$rem_nps_debug <- renderText({
    state <- rem_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, rem_nps_fit())
  })

    output$tte_nps_debug <- renderText({
    state <- tte_nps_debug()
    if (identical(state$status, "running")) invalidateLater(1000, session)
    format_model_debug(state, tte_nps_fit())
  })

  # ================================================================
  # COMPARE & COMBINE (WAIC-based)
  # ================================================================
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
      note_text <- "Run at least one uploaded-data model from its tab first. The table will update as REM, TTE, and USCR finish."
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
    # Bayesian p-values outside 0.05-0.95 (failed posterior predictive check) in red.
    datatable(df, options = list(paging = FALSE, ordering = FALSE, dom = "t", autoWidth = TRUE), rownames = FALSE) |>
      DT::formatStyle(
        "Bayesian p-value",
        color = DT::styleInterval(c(0.05, 0.95), c("firebrick", "inherit", "firebrick")),
        fontWeight = DT::styleInterval(c(0.05, 0.95), c("bold", "normal", "bold"))
      )
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

  output$species_count_summary <- renderTable({
    tibble(
      category = c("Total detections", "Camera trap days"),
      result = c(
        sum(nps_model_inputs()$camera_counts),
        sum(nps_model_inputs()$camera_days)
        )
    )
  },
  colnames = F)
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------

shinyApp(ui, server)
