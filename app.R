# -------------------------------------------------------------------
# Packages
# -------------------------------------------------------------------

needs <- c(
  "shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr",
  "readr", "purrr",
  "leaflet", "ggrepel",
  "nimble", "parallel", "MCMCvis", "lubridate", "sf", "tibble"
)

missing <- needs[!needs %in% rownames(installed.packages())]
if (length(missing)) {
  message("DEER App: installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
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

# -------------------------------------------------------------------
# Colors
# -------------------------------------------------------------------

# Redwood-inspired palette (hex only — no extra color package)
redwood_colors <- c(
  "#303018", "#604830", "#609048", "#90A860", "#786048", "#B8C4A8"
)

# -------------------------------------------------------------------
# Default MCMC settings (not user-facing)
# -------------------------------------------------------------------
# Starting values for every model run. Adaptive reruns double iterations,
# burn-in, and thinning if R-hat fails; double M (USCR) or D_max (REM, TTE)
# if the posterior presses against it; and, for USCR, grow the state-space
# buffer if it is too small. Users can change only the maximum number of
# rounds.
mcmc_defaults <- list(
  n_chains         = 2L,
  iter             = 6000L,
  burnin           = 1000L,
  thin             = 5L,
  M_uscr           = 100L,
  max_adapt_rounds = 6L
)

# Checks after each round of an adaptive run. A round passes when the largest
# R-hat is at most rhat_max, and when the posterior probability that
# abundance (USCR) or density (REM, TTE) exceeds limit_fraction of its upper
# limit (M or D_max) is at most limit_prob.
adapt_checks <- list(
  rhat_max       = 1.1,
  limit_fraction = 0.9,
  limit_prob     = 0.01
)

# -------------------------------------------------------------------
# Defaults
# -------------------------------------------------------------------
# REM/TTE upper limit on density under the uniform prior (animals/km^2). Not
# user-facing: the app doubles it if the posterior for D presses against it.
rem_tte_start_D_max <- 200

# Default priors: the Advanced inputs' starting values, always used in
# Default mode.
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

# Default detection angle (degrees) for uploaded cameras without a valid
# Camera Detection Angle value.
default_theta_deg <- 55

# Number of days to keep when the deployment trim is turned on (the trim is
# off by default).
default_trim_days <- 56

# Bayesian p-values (posterior predictive check) inside this range pass; a
# model whose p-value falls outside it fails the check.
ppc_pass_range <- c(0.05, 0.95)

# Decimal places for model results.
result_digits <- 3

# Starting density threshold (animals/km²) for P(D > threshold) on the
# Compare & combine tab.
default_density_threshold <- 8

# R-hat for MathJax inside markdown(): markdown turns "\\(" into "\(", which
# MathJax then renders.
rhat_md <- "\\\\(\\hat{R}\\\\)"

# Format a number with a fixed number of decimal places ("NA" if it is
# missing or not finite).
format_fixed <- function(x, digits = 2) {
  if (length(x) != 1 || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

# Show every number in a results table (a datatable() of df) with
# result_digits decimal places. Unlike round(), this keeps trailing zeros
# (31.470, not 31.47); missing values stay blank.
format_result_columns <- function(dt, df) {
  formatRound(dt, names(df)[vapply(df, is.numeric, logical(1))], digits = result_digits)
}

# Round a USCR state-space buffer (m) up to the nearest 50 m.
round_buffer_m <- function(buffer_m) ceiling(buffer_m / 50) * 50

# Default USCR buffer (m) under the default log(sigma) prior.
default_uscr_buffer_m <- round_buffer_m(uscr_default_buffer_m(
  advanced_prior_defaults$log_sigma_mean,
  advanced_prior_defaults$log_sigma_sd
))

# The rule behind the default USCR buffer, for help text.
uscr_buffer_rule_text <- paste0(
  uscr_buffer_multiplier, " times the ", 100 * uscr_buffer_p,
  "% home-range radius at the ", 100 * uscr_buffer_p,
  "th percentile of the log(σ) prior, rounded up to the nearest 50 m"
)

# 95% circular home-range area (km^2) at the default prior mean of
# log(sigma), and its approximate 95% prior interval, quoted in the model
# descriptions.
default_hr95_km2 <- pi * uscr_home_range_radius_km(
  exp(advanced_prior_defaults$log_sigma_mean),
  p = 0.95
)^2
default_hr95_interval_km2 <- pi * uscr_home_range_radius_km(
  exp(
    advanced_prior_defaults$log_sigma_mean +
      c(-1, 1) * qnorm(0.975) * advanced_prior_defaults$log_sigma_sd
  ),
  p = 0.95
)^2

# Daily movement rate v (km/day) at the default prior mean of log(v), and its
# approximate 95% prior interval, quoted in the model descriptions.
default_v_km_day <- exp(advanced_prior_defaults$log_v_mean)
default_v_interval_km_day <- exp(
  advanced_prior_defaults$log_v_mean +
    c(-1, 1) * qnorm(0.975) * advanced_prior_defaults$log_v_sd
)

# -------------------------------------------------------------------
# Model summaries and model averaging
# -------------------------------------------------------------------

# Mixture of the models' posteriors: each draw comes from model (column) k
# with probability weights[k], as a random posterior draw from that model.
sample_mixture <- function(draws, weights) {
  n_draws <- nrow(draws)
  draws[cbind(
    sample.int(n_draws, n_draws, replace = TRUE),
    sample.int(ncol(draws), n_draws, replace = TRUE, prob = weights)
  )]
}

# Combines posterior draws of density across the models in `fits`, a named
# list in which models not run yet are NULL. Returns NULL if no model has been
# run; otherwise the fits that have been, the species they were fit to
# (results are cleared when the species changes, so the fits share one), their
# WAIC weights (NA if any model's WAIC is missing), and a matrix of density
# draws with one column per model followed by one column per model average.
prepare_combo_density_data <- function(fits) {
  fits <- compact(fits)
  if (!length(fits)) return(NULL)

  # WAIC weights: w_k = exp(-delta_k / 2) / sum_j exp(-delta_j / 2), where
  # delta_k = WAIC_k - min(WAIC).
  waic <- vapply(fits, \(fit) fit$waic, numeric(1))
  weights <- if (length(fits) == 1L) {
    1
  } else if (all(is.finite(waic))) {
    rel_lik <- exp(-0.5 * (waic - min(waic)))
    unname(rel_lik / sum(rel_lik))
  } else {
    rep(NA_real_, length(fits))
  }

  # Row s holds each model's s-th smallest draw, so averaging across a row
  # averages the models' quantiles. That needs the same number of draws from
  # every model, which holds because all three models start from the same
  # MCMC settings, and adaptive reruns double iterations, burn-in, and
  # thinning together, leaving (iterations - burn-in) / thinning unchanged.
  draws <- lapply(fits, \(fit) sort(fit$samples_all[, "D"]))
  stopifnot(length(unique(lengths(draws))) == 1L)
  draws <- do.call(cbind, draws)

  density_est <- draws
  if (ncol(draws) > 1L) {
    density_est <- cbind(
      density_est,
      `Unweighted quantile average` = rowMeans(draws)
    )
    if (!anyNA(weights)) {
      density_est <- cbind(
        density_est,
        `WAIC-weighted quantile average` = as.numeric(draws %*% weights),
        `WAIC-weighted mixture` = sample_mixture(draws, weights)
      )
    }
  }

  list(
    fits = fits,
    species = fits[[1]]$species,
    weights = weights,
    density_est = density_est
  )
}

# Density estimates, P(D > density_threshold), WAIC weights, and Bayesian
# p-values for each model and model average. The failed_check column marks
# models that failed their posterior predictive check (hidden in the app,
# where it sets the p-value's highlight).
build_combo_table <- function(combo_data, density_threshold, ppc_pass_range) {
  density_est <- combo_data$density_est
  n_averages <- ncol(density_est) - length(combo_data$fits)
  ppc <- lapply(combo_data$fits, ppc_summary, pass_range = ppc_pass_range)
  column_stat <- function(f, ...) unname(apply(density_est, 2, f, ...))
  threshold_label <- paste0(
    "P(D > ",
    format(density_threshold, trim = TRUE, scientific = FALSE),
    " animals/km²)"
  )

  tibble(
    Method = colnames(density_est),
    `Mean density (animals/km²)` = column_stat(mean),
    `Median density (animals/km²)` = column_stat(median),
    `Lower 2.5%` = column_stat(quantile, probs = 0.025),
    `Upper 97.5%` = column_stat(quantile, probs = 0.975),
    !!threshold_label := column_stat(\(x) mean(x > density_threshold)),
    `WAIC weight` = c(combo_data$weights, rep(NA_real_, n_averages)),
    `Bayesian p-value` = c(
      vapply(ppc, \(x) x$p, numeric(1), USE.NAMES = FALSE),
      rep(NA_real_, n_averages)
    ),
    failed_check = c(
      vapply(ppc, \(x) x$fail, logical(1), USE.NAMES = FALSE),
      rep(FALSE, n_averages)
    )
  )
}

# Colors and point shapes for each model and model average in the plots.
combo_plot_colors <- c(
  "USCR" = "#4F6FA3",
  "REM" = "#8C4F2D",
  "TTE" = "#3F7F5F",
  "Unweighted quantile average" = "#8A8677",
  "WAIC-weighted quantile average" = "#2F3E46",
  "WAIC-weighted mixture" = "#8E44AD"
)
combo_plot_shapes <- c(
  "USCR" = 16,
  "REM" = 16,
  "TTE" = 16,
  "Unweighted quantile average" = 17,
  "WAIC-weighted quantile average" = 15,
  "WAIC-weighted mixture" = 18
)

# Posterior mean and 95% credible interval of density for each model and
# model average.
build_combo_interval_plot <- function(combo_data, title_text) {
  density_est <- combo_data$density_est
  interval_df <- tibble(
    # Levels reversed so the first column is plotted at the top.
    Method = factor(colnames(density_est), levels = rev(colnames(density_est))),
    mean = unname(colMeans(density_est)),
    lower = unname(apply(density_est, 2, quantile, probs = 0.025)),
    upper = unname(apply(density_est, 2, quantile, probs = 0.975))
  )

  ggplot(
    interval_df,
    aes(x = mean, y = Method, color = Method, shape = Method)
  ) +
    geom_errorbar(
      aes(xmin = lower, xmax = upper),
      orientation = "y",
      width = 0.18,
      linewidth = 1,
      alpha = 0.95
    ) +
    geom_point(size = 3.2) +
    scale_color_manual(values = combo_plot_colors) +
    scale_shape_manual(values = combo_plot_shapes) +
    labs(title = title_text, x = "Density (animals/km²)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

# Posterior distributions of density for each model, overlaid.
build_combo_density_overlay_plot <- function(combo_data, title_text) {
  model_names <- names(combo_data$fits)
  density_df <- combo_data$density_est[, model_names, drop = FALSE] |>
    as_tibble() |>
    pivot_longer(
      everything(),
      names_to = "Method",
      values_to = "density_km2"
    ) |>
    mutate(Method = factor(Method, levels = model_names))

  ggplot(density_df, aes(x = density_km2, color = Method, fill = Method)) +
    geom_density(alpha = 0.16, linewidth = 1) +
    scale_color_manual(values = combo_plot_colors) +
    scale_fill_manual(values = combo_plot_colors) +
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

# Posterior summary of a fit's parameters for the CSV download, one row per
# parameter. D is renamed D_km2 so the units are explicit. log_sigma and
# log_lam_0 (USCR) duplicate sigma and lam_0 on the log scale, so they are
# dropped. bp is 0 or 1 in each draw, so only its mean (the Bayesian p-value)
# is reported.
posterior_summary_df <- function(fit) {
  draws <- as_tibble(fit$samples_all) |>
    rename(D_km2 = D) |>
    select(-any_of(c("log_sigma", "log_lam_0"))) |>
    relocate(any_of(
      c("D_km2", "v", "sd_eps", "N", "psi", "sigma", "lam_0", "bp")
    ))

  tibble(
    parameter = names(draws),
    mean = vapply(draws, mean, numeric(1)),
    median = vapply(draws, median, numeric(1)),
    q025 = vapply(draws, quantile, numeric(1), probs = 0.025),
    q975 = vapply(draws, quantile, numeric(1), probs = 0.975)
  ) |>
    mutate(
      across(
        c(median, q025, q975),
        \(x) if_else(parameter == "bp", NA_real_, x)
      ),
      across(where(is.numeric), \(x) round(x, result_digits))
    )
}

# Explains the posterior summary CSV, below its download button.
posterior_param_legend <- function() {
  markdown(paste(
    "**Posterior summary CSV:** one row per model and parameter, with the species analyzed and the posterior mean, median, and 2.5% and 97.5% quantiles (`q025`, `q975`).",
    "",
    "- `D_km2` — density (animals/km²).",
    "- `v` — REM and TTE only: daily movement rate (km/day).",
    "- `sd_eps` (σ<sub>ε</sub>) — standard deviation of the camera-level random effects. Only included when camera-level noise is turned on (**Model settings > Advanced > Priors: Camera-level noise (all models)**).",
    "- `N` — USCR only: abundance, the number of animals in the state space.",
    "- `psi` (ψ) — USCR only: probability that an augmented individual is part of the real population. Its prior sets the prior on abundance and density.",
    "- `sigma` (σ) — USCR only: spatial scale parameter (km), related to home-range size and how quickly detection falls with distance.",
    "- `lam_0` (λ₀) — USCR only: expected count of an individual per camera-day if its activity center were at the camera.",
    paste0(
      "- `bp` — posterior predictive check. Its mean is the Bayesian p-value; the check fails if it is below ",
      ppc_pass_range[1], " or above ", ppc_pass_range[2],
      ". `bp` is 0 or 1 in each draw, so its median and quantiles are left blank."
    ),
    sep = "\n"
  ))
}

# -------------------------------------------------------------------
# Upload checks (Add your data tab)
# -------------------------------------------------------------------

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

# Messages and warnings from several QC steps, in order, as one
# list(messages, warnings).
combine_issues <- function(...) {
  steps <- list(...)
  list(
    messages = c(character(), unlist(lapply(steps, `[[`, "messages"))),
    warnings = c(character(), unlist(lapply(steps, `[[`, "warnings")))
  )
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

# Legend breaks for whole-number counts: pretty() can return halves (for
# example 1, 1.5, 2), which don't occur in count data.
integer_breaks <- function(limits) {
  breaks <- pretty(limits)
  breaks[breaks == round(breaks)]
}

# -------------------------------------------------------------------
# Model inputs and settings
# -------------------------------------------------------------------

# Which uploaded cameras use the default detection angle (Model settings):
# those with no Camera Detection Angle, or one outside the range the upload
# check accepts (the check warns about these but doesn't block the upload).
camera_angle_uses_default <- function(out) {
  if (!"Camera Detection Angle" %in% names(out)) {
    return(rep(TRUE, nrow(out)))
  }
  theta <- suppressWarnings(as.numeric(out$`Camera Detection Angle`))
  is.na(theta) | theta < camera_angle_min_deg | theta > camera_angle_max_deg
}

# Detection angle (degrees) for each uploaded camera.
get_camera_angle_vector <- function(out, default_deg) {
  uses_default <- camera_angle_uses_default(out)
  theta <- rep(as.numeric(default_deg), nrow(out))
  # Read the column only when some cameras have a valid angle: it may not
  # exist, and a tibble warns about a missing column.
  if (!all(uses_default)) {
    theta[!uses_default] <- as.numeric(out$`Camera Detection Angle`[!uses_default])
  }
  theta
}

# Input-summary line describing the uploaded cameras' detection angles.
summarize_camera_angle_source <- function(out, default_deg) {
  theta <- get_camera_angle_vector(out, default_deg)
  n_default <- sum(camera_angle_uses_default(out))
  range_text <- paste0(
    format_fixed(min(theta), 1), " to ", format_fixed(max(theta), 1), " degrees"
  )
  source_text <- if (n_default == nrow(out)) {
    "(all cameras using the default angle from Model settings)"
  } else if (n_default > 0) {
    paste0("(", n_default, " camera(s) using the default angle from Model settings)")
  } else {
    "(from uploaded deployment data)"
  }
  paste("Camera detection angle:", range_text, source_text)
}

is_valid_buffer_m <- function(x) {
  is.numeric(x) && length(x) == 1 && is.finite(x) && x > 0
}
uscr_buffer_problem <- "The state-space buffer (Model settings > Advanced > Priors: USCR) must be a positive number of meters."

# Density prior specs come from Model settings: list(uniform = TRUE), or
# list(uniform = FALSE, lo = , hi = ), a 95% interval in animals/km^2.
density_prior_problem <- function(spec) {
  if (isTRUE(spec$uniform)) return(NULL)
  check_density_prior(spec$lo, spec$hi)
}
describe_density_prior_spec <- function(spec) {
  if (isTRUE(spec$uniform)) return("uniform (default)")
  paste0("informative, 95% interval ", spec$lo, " to ", spec$hi, " animals/km²")
}

# Starting upper limit on density (animals/km^2) for a REM or TTE run.
encounter_start_D_max <- function(density_prior) {
  if (isTRUE(density_prior$uniform)) {
    rem_tte_start_D_max
  } else {
    encounter_prior_start_D_max(density_prior$hi)
  }
}

# Starting M for a USCR run: mcmc_defaults$M_uscr, raised under an
# informative density prior so the prior fits below it (run_USCR() raises it
# the same way in later rounds if the state space grows).
uscr_start_M <- function(out, buffer_m, density_prior) {
  if (isTRUE(density_prior$uniform)) return(mcmc_defaults$M_uscr)
  area_km2 <- uscr_state_space_and_area(out, buffer_m)$area_km2
  max(mcmc_defaults$M_uscr, uscr_prior_start_M(density_prior$hi, area_km2))
}

noise_label <- function(on) {
  if (isTRUE(on)) "on (Poisson-lognormal counts)" else "off (Poisson counts)"
}

# Input-summary lines for a model run: the data, model_lines (the model's own
# settings), the starting MCMC settings (limit_setting is "M=..." for USCR or
# "D_max=..." for REM and TTE), and how the adaptive reruns work.
run_context_lines <- function(d, model, model_lines, limit_setting, max_adapt_rounds) {
  is_uscr <- identical(model, "USCR")
  c(
    paste("Cameras:", nrow(d$out)),
    paste("Animals counted:", sum(d$camera_counts, na.rm = TRUE)),
    paste("Total camera-days:", format_fixed(sum(d$camera_days, na.rm = TRUE), 1)),
    model_lines,
    paste0(
      "Initial MCMC settings: iter=", mcmc_defaults$iter,
      ", burn-in=", mcmc_defaults$burnin,
      ", thin=", mcmc_defaults$thin,
      ", chains=", mcmc_defaults$n_chains,
      ", ", limit_setting
    ),
    paste0(
      model, " reruns automatically, up to ", max_adapt_rounds,
      " rounds in all, at the same chain count shown above. ",
      "Each check that fails triggers its own change: iterations, burn-in, and ",
      "thinning double when the maximum R-hat is above ", adapt_checks$rhat_max, "; ",
      if (is_uscr) "M" else "D_max",
      " doubles when the posterior probability that ",
      if (is_uscr) "N" else "D",
      " exceeds ", 100 * adapt_checks$limit_fraction, "% of ",
      if (is_uscr) "M" else "D_max",
      " is above ", 100 * adapt_checks$limit_prob, "%",
      if (is_uscr) {
        paste0(
          "; and the state-space buffer grows when it is less than ",
          uscr_buffer_multiplier, "× the ", 100 * uscr_buffer_p,
          "% home-range radius at the largest posterior σ"
        )
      },
      "."
    )
  )
}

# -------------------------------------------------------------------
# Model runs and run-status panels
# -------------------------------------------------------------------

# Starting state of a model's run-status panel. data_label is "simulated
# data" or "uploaded data".
make_model_debug <- function(model, data_label) {
  list(
    model = model,
    data_source = data_label,
    started_at = NULL,
    finished_at = NULL,
    guidance = "No run has been started yet.",
    raw_error = NULL,
    context = character(),
    history = character()
  )
}

# Update the run-status panel held in the reactiveVal rv, adding log_entry
# (if given) to its log.
update_model_debug <- function(rv,
                               started_at = NULL,
                               finished_at = NULL,
                               guidance = NULL,
                               raw_error = NULL,
                               context = NULL,
                               log_entry = NULL) {
  state <- rv()
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

friendly_model_error <- function(model, data_label, error_text) {
  tips <- c(
    paste0(model, " on ", data_label, " did not finish."),
    "Read the raw error below, then try the model-specific guidance."
  )

  if (grepl("cannot allocate|vector memory|std::bad_alloc", error_text, ignore.case = TRUE)) {
    tips <- c(tips, "This looks like a memory issue: the computer running the app ran out of memory. For USCR, memory use grows with M, which doubles in each adaptive round where the posterior for abundance presses against it. An informative density prior (Model settings > Advanced > Priors: Density (all models)) can keep M from growing, and a lower 'Maximum adaptive rounds' in Model settings limits how many times it can double.")
  }
  if (grepl("there is no package called", error_text, ignore.case = TRUE)) {
    tips <- c(tips, "A required R package is missing in this session. Install the package named in the raw error, restart R, and rerun.")
  }
  if (grepl("\\b(NA|NaN|Inf)\\b", error_text, perl = TRUE)) {
    tips <- c(tips, "Missing or invalid numeric values reached the model. Check detection distances, camera-days, timestamps, and coordinates on the Add your data and Data summary tabs.")
  }

  if (identical(model, "USCR")) {
    tips <- c(tips, "USCR doubles iterations, burn-in, and thinning automatically when R-hat is too high, doubles M when the posterior for abundance presses against it, and grows the state-space buffer when the posterior σ is larger than the buffer safely covers, so longer runs can be expected when any check fails. To cap run time, lower 'Maximum adaptive rounds' in Model settings (results from a run that stops at the limit may not be reliable). An informative density prior (Model settings > Advanced > Priors: Density (all models)) can also keep M, and so run time, from growing.")
  }
  if (model %in% c("REM", "TTE")) {
    tips <- c(tips, paste0(model, " expects per-camera animals counted, camera-days, and detection distances."))
    tips <- c(tips, paste0(model, " doubles iterations, burn-in, and thinning automatically when R-hat is too high, and doubles D_max when the posterior for density presses against it, so longer runs can be expected when either check fails. To cap run time, lower 'Maximum adaptive rounds' in Model settings (results from a run that stops at the limit may not be reliable)."))
  }

  paste(tips, collapse = "\n")
}

# Status callback for a model run. The model engine calls it once per MCMC
# round with the round number and this round's settings. From the second
# round on (a rerun because a check failed), it also gets the previous
# round's elapsed time and iterations, so the pop-up can show a rough time
# estimate. USCR also passes M (its augmentation bound), the previous round's
# M, and buffer_m (state-space buffer); REM/TTE pass D_max (their upper limit
# on density) instead.
#
# The estimate scales the previous round's time by the change in
# iterations and, for USCR, by (M / previous M)^2: each update to one
# augmented individual recomputes every camera's expected count, a sum over
# all M individuals, so the work per iteration grows roughly with M^2.
make_model_status_callback <- function(rv, model_label, notification_id) {
  function(round, iter, burnin, M = NULL, buffer_m = NULL, D_max = NULL,
           prev_elapsed_min = NULL, prev_iter = NULL, prev_M = NULL) {
    message <- paste0(model_label, ": ")
    if (!is.null(M)) message <- paste0(message, "M=", M, ", ")
    if (!is.null(D_max)) message <- paste0(message, "D_max=", format_fixed(D_max, 1), " animals/km², ")
    if (!is.null(buffer_m)) message <- paste0(message, "buffer=", round(buffer_m), " m, ")
    message <- paste0(message, iter, " iterations, ", burnin, " burn-in")
    if (round > 1L && !is.null(prev_elapsed_min) && !is.null(prev_iter) && prev_iter > 0) {
      est_min <- prev_elapsed_min * (iter / prev_iter)
      if (!is.null(M) && !is.null(prev_M) && prev_M > 0) {
        est_min <- est_min * (M / prev_M)^2
      }
      message <- paste0(message, " — approx. ", format(round(est_min, 1), nsmall = 1), " min")
    }
    showNotification(
      message,
      type = "message",
      duration = NULL,
      id = notification_id
    )
    update_model_debug(rv, log_entry = message)
  }
}

# Posterior predictive check text for a fit (see ppc_summary()).
ppc_text <- function(fit) {
  ppc <- ppc_summary(fit, ppc_pass_range)
  if (is.null(ppc)) return(NULL)
  base <- paste0(
    "Posterior predictive check: Bayesian p-value = ",
    format_fixed(ppc$p, result_digits)
  )
  if (!ppc$fail) {
    return(paste0(base, " (passes; values from ", ppc$lower, " to ", ppc$upper, " pass)."))
  }
  # Camera-level noise adds variation between cameras, so suggest it only
  # when the observed counts vary more than the model expects.
  paste0(
    base, ". FAILED (outside ", ppc$lower, " to ", ppc$upper, "): the observed counts vary ",
    if (ppc$p < ppc$lower) "more" else "less",
    " than the model expects, so its estimates and intervals may be unreliable.",
    if (!isTRUE(fit$settings$camera_noise) && ppc$p < ppc$lower) {
      " Turning on camera-level noise (Model settings > Advanced > Priors: Camera-level noise (all models)) can often fix a failed check."
    } else {
      ""
    }
  )
}

notify_ppc <- function(fit, label) {
  ppc <- ppc_summary(fit, ppc_pass_range)
  if (!is.null(ppc) && isTRUE(ppc$fail)) {
    showNotification(
      paste0(
        label, " failed its posterior predictive check (Bayesian p-value ",
        format_fixed(ppc$p, result_digits), "). See the run-status panel."
      ),
      type = "warning", duration = NULL
    )
  }
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

# Run one model fit, shared by the six Run buttons (USCR, REM, and TTE on
# simulated and uploaded data). It clears the previous result in fit_rv,
# restarts the run-status panel in debug_rv with the input summary
# (context), shows progress pop-ups, and reports success or failure.
# fit_model(status_callback) runs the model and returns the fit.
# extra_warning(fit), if given, returns a warning to add to the guidance and
# show as a pop-up, or NULL.
run_model <- function(model,
                      data_label,
                      fit_rv,
                      debug_rv,
                      context,
                      fit_model,
                      extra_warning = NULL) {
  label <- paste0(model, " (", data_label, ")")
  notification_id <- paste0(tolower(model), "_", gsub(" ", "_", data_label), "_status")

  fit_rv(NULL)
  debug_rv(make_model_debug(model, data_label))
  update_model_debug(
    debug_rv,
    started_at = Sys.time(),
    context = context,
    log_entry = paste0("Input checks passed. Starting ", model, " MCMC.")
  )
  showNotification(
    paste0("Running ", model, " on ", data_label, "."),
    type = "message",
    duration = NULL,
    id = notification_id
  )

  fit <- tryCatch(
    fit_model(make_model_status_callback(debug_rv, model, notification_id)),
    error = function(e) {
      update_model_debug(
        debug_rv,
        finished_at = Sys.time(),
        guidance = friendly_model_error(model, data_label, conditionMessage(e)),
        raw_error = conditionMessage(e),
        log_entry = paste(label, "failed:", conditionMessage(e))
      )
      showNotification(
        paste(label, "failed:", conditionMessage(e)),
        type = "error", duration = NULL
      )
      NULL
    },
    finally = removeNotification(notification_id)
  )
  fit_rv(fit)
  if (is.null(fit)) return(invisible(NULL))

  warning_text <- if (!is.null(extra_warning)) extra_warning(fit)
  guidance <- adaptive_limit_guidance(
    fit, model,
    paste0(
      model, " on ", data_label,
      " completed successfully. Review the results above and the run history below."
    )
  )
  if (!is.null(warning_text)) guidance <- paste(guidance, warning_text)
  update_model_debug(
    debug_rv,
    finished_at = Sys.time(),
    guidance = guidance,
    log_entry = paste(label, "run completed.")
  )
  showNotification(paste0(label, " complete!"), type = "message")
  if (!is.null(warning_text)) {
    showNotification(paste0(label, ": ", warning_text), type = "warning", duration = NULL)
  }
  notify_adaptive_limit(fit, label)
  notify_ppc(fit, label)
  invisible(fit)
}

# Text for a run-status panel.
format_model_debug <- function(state, fit = NULL) {
  lines <- c(
    paste("Model:", state$model),
    paste("Data source:", state$data_source)
  )

  if (!is.null(state$started_at)) {
    lines <- c(lines, paste("Started:", format(state$started_at, "%Y-%m-%d %H:%M:%S")))
  }
  if (!is.null(state$finished_at)) {
    lines <- c(lines, paste("Finished:", format(state$finished_at, "%Y-%m-%d %H:%M:%S")))
  }
  if (!is.null(state$started_at) && !is.null(state$finished_at)) {
    elapsed_min <- as.numeric(difftime(state$finished_at, state$started_at, units = "mins"))
    lines <- c(lines, paste("Elapsed:", format_fixed(elapsed_min, 1), "minutes"))
  }

  lines <- c(lines, "", "Guidance:", state$guidance)

  if (length(state$context)) {
    lines <- c(lines, "", "Input summary:", paste0("- ", state$context))
  }

  if (!is.null(fit)) {
    setting_lines <- c(
      paste0("- Iterations: ", fit$settings$iter),
      paste0("- Burn-in: ", fit$settings$burnin),
      paste0("- Thinning: ", fit$settings$thin),
      paste0("- Chains: ", fit$settings$n_chains)
    )
    if (!is.null(fit$settings$M)) {
      setting_lines <- c(paste0("- M: ", fit$settings$M), setting_lines)
    }
    if (!is.null(fit$settings$buffer_m)) {
      setting_lines <- c(setting_lines, paste0("- State-space buffer: ", round(fit$settings$buffer_m), " m"))
    }
    if (!is.null(fit$settings$D_max)) {
      setting_lines <- c(setting_lines, paste0("- D_max: ", format_fixed(fit$settings$D_max, 1), " animals/km²"))
    }
    prior_text <- if (!is.null(fit$density_prior$M)) {
      describe_uscr_density_prior(fit$density_prior)
    } else {
      describe_beta_density_prior(fit$density_prior)
    }
    setting_lines <- c(setting_lines, paste0("- Density prior: ", prior_text))
    if (isTRUE(fit$density_prior$off_target)) {
      setting_lines <- c(setting_lines, paste0(
        "  (a bound of the prior's 95% interval differs from the target by more than ",
        round(100 * fit$density_prior$miss_tolerance), "% of the interval's width)"
      ))
    }
    setting_lines <- c(setting_lines, paste0("- Camera-level noise: ", noise_label(fit$settings$camera_noise)))
    if (is.finite(fit$final_rhat_max)) {
      setting_lines <- c(setting_lines, paste0("- Max R-hat: ", format_fixed(fit$final_rhat_max, result_digits)))
    }
    lines <- c(lines, "", "Final round:", setting_lines)
    ppc_line <- ppc_text(fit)
    if (!is.null(ppc_line)) {
      lines <- c(lines, "", ppc_line)
    }

    # One line per round. USCR rounds have M, buffer_m, M_too_small, and
    # buffer_too_small; REM and TTE rounds have D_max and D_too_small.
    round_df <- fit$round_history
    round_lines <- vapply(seq_len(nrow(round_df)), function(i) {
      row <- round_df[i, ]
      paste(c(
        paste0("round ", row$round),
        if (!is.null(row$M)) paste0("M=", row$M),
        if (!is.null(row$D_max)) paste0("D_max=", format_fixed(row$D_max, 1)),
        if (!is.null(row$buffer_m)) paste0("buffer=", round(row$buffer_m), " m"),
        paste0("iter=", row$niter),
        paste0("thin=", row$thin),
        paste0("chains=", row$n_chains),
        paste0("R-hat=", format_fixed(row$rhat_max, result_digits)),
        if (!is.null(row$D_too_small)) paste0("D_max_too_small=", row$D_too_small),
        if (!is.null(row$M_too_small)) paste0("M_too_small=", row$M_too_small),
        if (!is.null(row$buffer_too_small)) paste0("buffer_too_small=", row$buffer_too_small)
      ), collapse = ", ")
    }, character(1))
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

# Posterior summary table for a model fit to uploaded data: the parameters
# the simulated-data tables show (see sim_truth_values()), without the true
# values.
uploaded_summary_table <- function(fit) {
  params <- c(
    "D_km2",
    intersect(c("sigma", "lam_0", "v", "sd_eps"), colnames(fit$samples_all))
  )
  map(params, function(param) {
    draws <- extract_param_draws(fit$samples_all, param)
    q <- quantile(draws, c(0.025, 0.975), na.rm = TRUE)
    tibble(
      Parameter = param,
      Mean = mean(draws, na.rm = TRUE),
      Median = median(draws, na.rm = TRUE),
      `Lower 2.5%` = q[[1]],
      `Upper 97.5%` = q[[2]]
    )
  }) |>
    bind_rows()
}

# Results for a model fit to uploaded data (a model tab): the species and
# WAIC above the posterior summary table, or a note if the model hasn't been
# run or its last run failed.
uploaded_result_ui <- function(model, fit, debug) {
  if (!is.null(debug$raw_error)) {
    return(tags$p(
      model, " (uploaded data): the last run failed. See ",
      tags$strong("Run status & troubleshooting"),
      " below for the raw error and guidance."
    ))
  }
  if (is.null(fit)) {
    return(tags$p(
      class = "small",
      model, " (uploaded data): not run yet. Click ",
      tags$strong(paste("Run", model, "on uploaded data"), .noWS = "after"),
      "."
    ))
  }
  tagList(
    tags$p(
      tags$strong("Species:"), fit$species,
      if (is.finite(fit$waic)) {
        tagList(" · ", tags$strong("WAIC:"), format_fixed(fit$waic, result_digits))
      }
    ),
    results_datatable(uploaded_summary_table(fit), fit)
  )
}

# -------------------------------------------------------------------
# Simulate data tab
# -------------------------------------------------------------------

# Simulated counts by camera, with each camera's point sized by its count.
sim_count_plot <- function(sim, model, D_true) {
  df <- sim$out |> mutate(count = sim$camera_counts)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(
      aes(size = count),
      shape = 21, fill = redwood_colors[3], color = redwood_colors[1], alpha = 0.75
    ) +
    geom_label_repel(aes(label = `Site Name`), size = 3, max.overlaps = Inf) +
    scale_size_area(name = "Animals counted", max_size = 14) +
    coord_equal() +
    theme_minimal() +
    labs(
      title = sprintf("%s simulated counts (true density: %.1f animals/km²)", model, D_true),
      x = "m", y = "m"
    )
}

# Table of a simulated dataset's cameras and counts.
sim_camera_datatable <- function(sim) {
  sim$out |>
    transmute(
      Camera = `Site Name`,
      `x (m)` = x,
      `y (m)` = y,
      `Animals counted` = sim$camera_counts,
      `Camera-days` = sim$camera_days
    ) |>
    datatable(rownames = FALSE, options = list(pageLength = 8))
}

# Input-summary lines giving the true values a simulated dataset (s) was
# generated with.
sim_truth_context <- function(s, model) {
  st <- s$settings
  c(
    paste("True density:", st$D_true, "animals/km²"),
    if (identical(model, "USCR")) {
      c(
        paste("True σ:", format_fixed(st$sigma_km, 3), "km"),
        paste("True λ₀:", st$lambda0, "per camera-day"),
        paste("Simulated activity centers:", s$uscr$N),
        paste("Simulation buffer:", round(st$buffer_m), "m")
      )
    } else {
      paste("True daily movement rate v:", st$v, "km/day")
    },
    if (isTRUE(st$camera_noise)) paste("True camera-level noise SD (σε):", st$sd_eps)
  )
}

# True parameter values for a model's truth-vs-estimate table (see
# summarize_vs_truth()), from a simulation's settings. sd_eps is included
# only with camera-level noise on.
sim_truth_values <- function(st, model) {
  c(
    list(D_km2 = st$D_true),
    if (identical(model, "USCR")) {
      list(sigma = st$sigma_km, lam_0 = st$lambda0)
    } else {
      list(v = st$v)
    },
    if (isTRUE(st$camera_noise)) list(sd_eps = st$sd_eps)
  )
}

# Results table for a model fit, from summarize_vs_truth() (simulated data)
# or uploaded_summary_table() (uploaded data): parameter symbols with units
# instead of internal names, Yes/No for coverage (simulated data only), and
# the posterior predictive check under the table.
results_datatable <- function(tbl, fit) {
  param_labels <- c(
    D_km2 = "D (animals/km²)",
    sigma = "σ (km)",
    lam_0 = "λ₀ (per camera-day)",
    v = "v (km/day)",
    sd_eps = "σ<sub>ε</sub>"
  )
  tbl <- tbl |>
    mutate(
      Parameter = coalesce(unname(param_labels[Parameter]), Parameter),
      across(any_of("Covered"), \(x) if_else(x, "Yes", "No"))
    )
  # Posterior predictive check under the table (red when it fails).
  ppc <- ppc_summary(fit, ppc_pass_range)
  caption <- if (!is.null(ppc)) {
    tags$caption(
      style = paste0(
        "caption-side: bottom; text-align: left; color: ",
        if (ppc$fail) "firebrick" else "var(--muted)", ";"
      ),
      ppc_text(fit)
    )
  }
  # escape = FALSE so the <sub> tag in the sigma_epsilon label renders.
  datatable(
    tbl,
    caption = caption,
    options = list(dom = "t", paging = FALSE),
    rownames = FALSE,
    escape = FALSE
  ) |>
    format_result_columns(tbl)
}

# Warning for a USCR fit on simulated data whose state space reaches past
# the region the activity centers were simulated in.
sim_buffer_warning <- function(fit_buffer_m, sim_buffer_m) {
  if (!is_valid_buffer_m(fit_buffer_m) || fit_buffer_m <= sim_buffer_m) return(NULL)
  paste0(
    "The USCR state-space buffer (", round(fit_buffer_m), " m) is larger than the simulation buffer (",
    round(sim_buffer_m), " m), so part of the fit's state space has no simulated activity centers, ",
    "and density will likely be underestimated. Increase the simulation buffer and click Simulate data again."
  )
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
      h1, h2, h3, h4, h5 {
        color: var(--ink) !important;
      }
      /* Wrap long lines in text-output boxes instead of scrolling sideways. */
      pre.shiny-text-output {
        white-space: pre-wrap;
        overflow-wrap: break-word;
      }
      .hero h1, .hero h2, .hero h3, .hero h4, .hero h5 {
        color: #609048 !important;
      }
      /* Heading scale: each level is smaller than the one above. Tabs use h2
         for their main sections, h3 for subsections, and h4 for minor labels. */
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
        font-size: 1.25rem;
        margin: 1rem 0 0.4rem;
      }
      h4 {
        font-size: 1.1rem;
        margin: 0.9rem 0 0.35rem;
      }
      h5 {
        font-size: 1rem;
        font-weight: 600;
        margin: 0.8rem 0 0.3rem;
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
      .prior-details, .data-details {
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

            tags$h3("Step 1: Add your data"),
            tags$p(
              "The main analysis workflow uses uploaded data. The app fits ",
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
              tags$li("Optionally trim each camera's analysis window to a set number of days (off by default).")
            ),

            tags$h3("Step 2: Review the data summary"),
            tags$p(
              "The ",
              tags$strong("Data summary"),
              " tab summarizes the uploaded data: camera-days and a map of the cameras, detections by species, and counts by site and by day for a selected species. Check it before fitting the models."
            ),

            tags$h3("Step 3: Choose the species and review the priors"),
            tags$p(
              "On the ",
              tags$strong("Model settings"),
              " tab, choose the ",
              tags$strong("species to analyze", .noWS = "after"),
              ". Until you choose one, the app selects the first species in alphabetical order."
            ),
            tags$p(
              "Model output is sensitive to the priors. The default priors are calibrated for white-tailed deer (",
              tags$em("Odocoileus virginianus", .noWS = "outside"),
              "), so you may need to change them for other species or study sites. To see and change the priors, choose ",
              tags$strong("Advanced"),
              " on the ",
              tags$strong("Model settings"),
              " tab. There you can:"
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

            tags$h3("Step 4: Run the models"),
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
                paste0(
                  "Runs are adaptive: after each round, the app checks convergence and whether its settings are large enough, adjusts them, and reruns until every check passes or the adaptive-round limit is reached (",
                  mcmc_defaults$max_adapt_rounds,
                  " rounds by default; change it in the "
                ),
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

            tags$h3("Step 5: Compare results"),
            tags$p(
              "Use the ",
              tags$strong("Compare & combine"),
              " tab to:"
            ),
            tags$ul(
              tags$li("Compare models fit to uploaded data using WAIC, and combine their density estimates (model averaging)."),
              tags$li("Check each model's posterior predictive check (Bayesian p-value)."),
              tags$li("Download posterior summaries with parameter names, means, medians, and 95% credible intervals.")
            ),

            tags$h3("Step 6: Optional design tools"),
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
              tags$li("Choose a study design (camera grid, spacing, survey length, true density), simulate data, and see whether each model can recover the known true values."),
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
                " Encounter rates are converted to density using the animals' movement rate and the camera viewshed."
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
                " Shorter effective time between encounters implies higher density when movement rate and viewshed are held fixed."
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
      tags$p(tags$strong("The app will:")),
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
        "Expand each section below to read the detailed data requirements.",
        sep = "\n"
      )),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Deployment file columns")),
        markdown(paste(
          "**Bold** columns are required; the others are optional, for quality control.",
          "",
          "- **`Site Name`**: Camera site identifier. Must be unique across all rows in the deployment file.",
          "- `Site`: Optional broader site identifier.",
          "- `Camera ID`: Optional camera identifier.",
          "- `SD Card ID`: Optional SD card identifier.",
          "- **`Start Date`**: Deployment start date (`MM/DD/YYYY`).",
          "- **`Start Time`**: Deployment start time (`HH:MM` or `HH:MM:SS`, 24-hour). If blank or unreadable, the camera's analysis window starts at midnight at the end of Start Date, so no time before the camera was set is counted.",
          "- **`End Date`**: Deployment end date (`MM/DD/YYYY`).",
          "- **`End Time`**: Deployment end time (`HH:MM` or `HH:MM:SS`, 24-hour). If blank or unreadable, the analysis window ends at midnight at the start of End Date, so no time after the camera was picked up is counted.",
          "- **`Latitude`**: Camera latitude (decimal degrees, WGS84, the standard datum used by GPS devices and phones). The app corrects the sign to match the hemisphere selected in **Model settings** (under Camera geometry and coordinates).",
          "- **`Longitude`**: Camera longitude (decimal degrees, WGS84). The app corrects the sign to match the hemisphere selected in **Model settings** (under Camera geometry and coordinates).",
          "- `Camera Model`: Optional camera make and model, for recordkeeping (for example, `Browning Strike Force Pro`).",
          paste0("- `Camera Detection Angle`: Optional full detection angle in degrees for that camera (from ", camera_angle_min_deg, " to ", camera_angle_max_deg, "). If left blank or outside that range, the app uses the default angle from **Model settings**."),
          "- `Camera Height`: Optional camera height, in any unit, for recordkeeping (the models don't use it).",
          "- `Camera Orientation`: Optional cardinal direction (`N`, `NE`, `E`, `SE`, `S`, `SW`, `W`, or `NW`) or 0-359 degrees.",
          "- **`Camera Functioning`**: Camera status; `Yes/No`, `TRUE/FALSE`, `T/F`, and `1/0` values are accepted.",
          "- **`Camera Malfunction Date`**: Leave blank for cameras that did not malfunction (a date entered while Camera Functioning is `Yes` is ignored, with a warning). If the camera failed (stopped working, or was knocked out of position but kept taking pictures), set Camera Functioning to `No` and enter the **first date whose data should be discarded** (`MM/DD/YYYY`). The camera's analysis window ends at midnight at the start of this date, so images from this date onward are excluded.",
          "- **`Detection Distance`**: Effective detection radius in meters; a positive number for every camera. REM and TTE depend on it.",
          sep = "\n"
        ))
      ),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Images file columns")),
        markdown(paste(
          "All columns are required unless marked optional.",
          "",
          "- **`Site Name`**: Camera site identifier. Must exactly match a Site Name in the deployment file.",
          "- **`Timestamp`**: Image date and time, written as month/day/year or year-month-day, with an optional 24-hour time (for example, `MM/DD/YYYY HH:MM:SS`, `YYYY-MM-DD HH:MM`, or `2/3/25 14:05`). Two-digit years are accepted.",
          "- **`Species`**: Species identifier. Pipe-delimited values such as `deer|squirrel` are allowed for multi-species rows; their order must match `Sighting Count`. A blank Species (for example, a trigger with no animal) is allowed, but that row isn't counted for any species.",
          "- **`Cluster ID`**: Identifier for an independent detection event (a \"cluster\"). Any format is accepted (text or numeric). Required on every row with a Species. Multiple rows can share a Cluster ID when they belong to the same event, for example burst photos, or a multi-species image split into one row per species.",
          "- **`Sighting Count`**: Number of individuals in the image. Pipe-delimited values such as `1|3` are allowed for multi-species rows, in the same order as `Species`. Required on every row with a Species.",
          "- `Image URL`: Optional image reference or link, stored as plain text, one value per row.",
          sep = "\n"
        ))
      ),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Notes")),
        markdown(paste(
          "- In the check logs, ❌ marks a problem that blocks the upload, ⚠️ marks a warning to review that doesn't block it, ℹ️ marks a note about data the models will leave out, 🛠 marks a value the app corrected (a Latitude or Longitude sign), and ✅ means the file passed with no issues.",
          "- Row numbers in the check logs count data rows: the first row after the column names is row 1.",
          "- **The deployment upload is blocked** if any of the following is true. These fields set each camera's analysis window and location, which every model depends on.",
          "  - A required column is missing.",
          "  - A Site Name is duplicated, or a row has deployment data but no Site Name.",
          "  - Start Date or End Date is missing or not in `MM/DD/YYYY` format, or End Date is before Start Date.",
          "  - A camera's analysis window (below) has no time in it.",
          "  - Latitude or Longitude is missing, non-numeric, or zero.",
          "  - Detection Distance is missing, non-numeric, or not positive.",
          "  - Camera Functioning is `No` and Camera Malfunction Date is missing, not in `MM/DD/YYYY` format, on or before Start Date, or after End Date.",
          "- **The images upload is blocked** if any of the following is true.",
          "  - No deployment file has been uploaded and checked.",
          "  - A required column is missing.",
          "  - A Site Name is blank or doesn't exactly match a Site Name in the deployment file.",
          "  - A Timestamp is missing or can't be read.",
          "  - A row with a Species has a missing or non-numeric Sighting Count, or a missing Cluster ID.",
          "  - A row's Species and Sighting Count have different numbers of pipe-delimited values.",
          "- Other problems, such as a missing Camera Functioning value, a Camera Malfunction Date on a camera marked as functioning, or image dates outside a camera's deployment dates, are reported as warnings.",
          "- Completely blank rows (for example, rows of commas at the end of a CSV file) are removed from either file, and the check log notes which rows. Empty columns with no column name, which spreadsheets can leave at the right edge of a CSV file, are removed without a note.",
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
        ))
      ),

      hr(),

      h2("Deployment data upload"),
      fileInput("deployment_csv", "Upload deployment CSV", accept = ".csv"),
      checkboxInput(
        "apply_trim",
        "Trim each camera's analysis window to a set number of days (24-hour periods) from its start date and time.",
        value = FALSE
      ),
      conditionalPanel(
        "input.apply_trim",
        numericInput("trim_days", "Number of days to keep", value = default_trim_days, min = 1, step = 1),
        uiOutput("trim_days_problem")
      ),
      h3("Deployment check log"),
      verbatimTextOutput("deployment_check_log"),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Preview of cleaned deployment data")),
        DTOutput("deployment_preview")
      ),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Analysis window for each camera")),
        tags$p(
          class = "small",
          "Each camera's analysis window runs from Start Date and Start Time to End Date and End Time. It's shortened by a missing Start Time or End Time, a Camera Malfunction Date, or the trim. Days is the length of the window, before time spent within detection events is removed."
        ),
        DTOutput("deployment_windows")
      ),
      tags$p(
        class = "small",
        "Optional: download the cleaned deployment data for your records. Once the check passes, the data are already loaded in the app, so there's nothing to re-upload."
      ),
      downloadButton("download_deployment_checked", "Download cleaned deployment CSV"),

      hr(),

      h2("Image data upload"),
      fileInput("images_csv", "Upload images CSV", accept = ".csv"),
      h3("Images check log"),
      verbatimTextOutput("images_check_log"),
      tags$details(
        class = "data-details",
        tags$summary(tags$strong("Preview of cleaned images data")),
        DTOutput("images_preview")
      ),
      tags$p(
        class = "small",
        "Optional: download the cleaned images data for your records. Once both checks pass, your data are loaded in the app and the models are ready to run."
      ),
      downloadButton("download_images_checked", "Download cleaned images CSV"),

      hr()
    ),

    # ---------------------- DATA SUMMARY --------------------------
    nav_panel(
      "Data summary",
      p(
        class = "small",
        style = "color: var(--muted);",
        "This tab summarizes uploaded deployment and image data. Summaries appear once the files are uploaded and pass their checks on the ",
        tags$strong("Add your data", .noWS = "after"),
        " tab: the deployment summary needs a deployment file, and the species summaries also need an images file. Simulated datasets are shown in the ",
        tags$strong("Simulate data"),
        " tab."
      ),
      h2("Site deployment summary"),
      tags$p(
        class = "small",
        "Days come from each camera's analysis window (see the ",
        tags$strong("Add your data", .noWS = "after"),
        " tab), which runs from Start Date and Start Time to End Date and End Time. It's shortened by a missing Start Time or End Time, a Camera Malfunction Date, or the trim."
      ),
      DTOutput("deploy_summary_table"),
      leafletOutput("camera_map", height = "300px"),
      hr(),
      h2("Detections by species"),
      tags$p(
        class = "small",
        tags$strong("Images"),
        " is the number of rows for that species in the images file. ",
        tags$strong("Detection events"),
        " is the number of independent detections: one per Cluster ID at a camera, counted only if its first image is inside the camera's analysis window. ",
        tags$strong("Animals counted"),
        " adds up the animals in those detection events, taking each event's largest Sighting Count; this is the count the models use. For example, a burst of 5 photos showing 2 deer is 5 images, 1 detection event, and 2 animals counted. ",
        tags$strong("Cameras with events"),
        " is the number of cameras with at least one detection event."
      ),
      DTOutput("image_summary_table"),
      plotOutput("species_bar_plot", height = "300px"),
      hr(),
      h2("Selected species by site"),
      selectInput(
        "summary_species",
        "Species to summarize",
        choices = character(0)
      ),
      DTOutput("species_site_table"),
      plotOutput("species_site_map", height = "300px"),
      plotOutput("species_daily_plot", height = "300px")
    ),

    # ---------------------- MODEL SETTINGS -------------------------
    nav_panel(
      "Model settings",

      h2("Species to analyze"),
      tags$p(
        class = "small",
        "Choose the species to estimate density for. The table shows the animals counted and camera-days the models will use."
      ),
      selectInput(
        inputId = "species_to_analyze",
        label = "Select species to analyze",
        choices = character(0)
      ),
      tableOutput("species_count_summary"),
      hr(),

      h2("Camera geometry and coordinates"),
      markdown(paste(
        paste0("Default detection angle is ", default_theta_deg, "° (based on Browning-style camera specifications). This value is used for uploaded cameras that do not include a `Camera Detection Angle` value. The **Simulate data** tab has its own detection angle setting."),
        "",
        "Longitude and Latitude signs are corrected to match the hemispheres selected below when the deployment file is checked (for example, a positive Longitude is made negative for Western sites, and a positive Latitude is made negative for Southern sites).",
        sep = "\n"
      )),
      sliderInput(
        "theta", "Default detection angle θ (degrees)",
        min = camera_angle_min_deg, max = camera_angle_max_deg, value = default_theta_deg, step = 1
      ),
      radioButtons(
        "hemisphere", "Hemisphere (for Longitude sign auto-correction)",
        choices = c("Western", "Eastern"),
        selected = "Western",
        inline = TRUE
      ),
      radioButtons(
        "lat_hemisphere", "Hemisphere (for Latitude sign auto-correction)",
        choices = c("Northern", "Southern"),
        selected = "Northern",
        inline = TRUE
      ),
      hr(),

      h2("Run defaults"),
      markdown(paste(
        paste0("- **Number of chains**: ", mcmc_defaults$n_chains),
        paste0(
          "- **All models**: adaptive runs start from ", mcmc_defaults$iter,
          " iterations, ", mcmc_defaults$burnin, " burn-in, and thin = ",
          mcmc_defaults$thin
        ),
        paste0("- **Convergence criteria**: ", rhat_md, " ≤ ", adapt_checks$rhat_max, " for all monitored parameters"),
        paste0(
          "- **Maximum adaptive rounds**: ", mcmc_defaults$max_adapt_rounds,
          " by default. Only the settings tied to a failed check change between rounds: iterations, burn-in, and thinning double if ", rhat_md, " fails; D<sub>max</sub> (REM, TTE) or M (USCR) doubles if the posterior presses against it; and for USCR the state-space buffer grows if it is too small. If a model reaches this limit without passing every check, increase it below."
        ),
        sep = "\n"
      )),
      numericInput(
        "max_adapt_rounds", "Maximum adaptive rounds per model run",
        value = mcmc_defaults$max_adapt_rounds, min = 1, step = 1
      ),
      hr(),

      h2("Prior distributions"),
      markdown(paste(
        "Switch to **Advanced** to see and change the priors, grouped by model type.",
        "The current defaults are calibrated for white-tailed deer (*Odocoileus virginianus*) using home-range size and daily movement rates derived from the literature.",
        "These defaults may not be well calibrated for other species or field conditions.",
        "Model output can be sensitive to prior choice, so users working on other species should adjust the priors to better match their system.",
        sep = "\n"
      )),
      radioButtons(
        "mode", "Settings mode",
        choices = c("Default", "Advanced"),
        selected = "Default",
        inline = TRUE
      ),

      conditionalPanel(
        "input.mode == 'Default'",
        tags$div(
          class = "about-card",
          tags$p(
            tags$strong("Using the app defaults."),
            " The models use the default priors described under Advanced, a uniform prior on density, and no camera-level noise."
          )
        )
      ),

      conditionalPanel(
        "input.mode == 'Advanced'",
        tagList(
          tags$p(
            class = "small",
            style = "color: var(--muted);",
            "The advanced controls focus on priors. This is where you can tune the models for your species and study site."
          ),
          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: Density (all models)")),
            markdown(paste(
              "One prior on density is shared by USCR, REM, and TTE, for uploaded and simulated data.",
              "",
              "**Default (uniform) priors.** Each model's default prior is flat, but the models put it on different quantities, with different upper limits:",
              "",
              paste0(
                "- **REM and TTE** put it directly on density: D ~ Uniform(0, D<sub>max</sub>). D<sub>max</sub> starts at ",
                rem_tte_start_D_max,
                " animals/km², and the app doubles it and reruns if the posterior for density presses against it."
              ),
              paste0(
                "- **USCR** puts it on abundance (N, the number of animals in the state space): N is uniform from 0 to M, and density is N divided by the state-space area. M starts at ",
                mcmc_defaults$M_uscr,
                " and doubles if the posterior for N presses against it."
              ),
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
              markdown(paste0(
                "- **σ<sub>ε</sub>** ~ Gamma(shape, rate): Prior on the standard deviation of the camera-level random effect. The default `Gamma(",
                advanced_prior_defaults$sd_eps_shape, ", ", advanced_prior_defaults$sd_eps_rate, ")` prior has mean ",
                signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate, 3), " and variance ",
                signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate^2, 3), "."
              )),
              fluidRow(
                column(
                  6,
                  numericInput(
                    "sd_eps_shape", HTML("σ<sub>ε</sub> gamma shape"),
                    value = advanced_prior_defaults$sd_eps_shape, min = 0.1, step = 0.1
                  )
                ),
                column(
                  6,
                  numericInput(
                    "sd_eps_rate", HTML("σ<sub>ε</sub> gamma rate"),
                    value = advanced_prior_defaults$sd_eps_rate, min = 0.1, step = 0.1
                  )
                )
              )
            )
          ),

          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: REM & TTE")),
            markdown(paste0(
              "- **log(v)** ~ Normal(mean, SD): Informative prior on daily movement rate (log scale, v in km/day). The current default corresponds to a daily movement rate of `",
              format_fixed(default_v_km_day), " km/day` at the prior mean of log(v), with an approximate 95% interval of `",
              format_fixed(default_v_interval_km_day[1]), "` to `", format_fixed(default_v_interval_km_day[2]),
              " km/day` (Christensen et al. in prep)."
            )),
            fluidRow(
              column(
                6,
                numericInput(
                  "log_v_mean", "log(v) mean (v in km/day)",
                  value = advanced_prior_defaults$log_v_mean, step = 0.1
                )
              ),
              column(
                6,
                numericInput(
                  "log_v_sd", "log(v) SD",
                  value = advanced_prior_defaults$log_v_sd, min = 0.01, step = 0.05
                )
              )
            )
          ),

          tags$details(
            class = "prior-details",
            tags$summary(tags$strong("Priors: USCR")),
            markdown(paste(
              paste0(
                "- **log(σ)** ~ Normal(mean, SD): Informative prior on animal space use (log scale, σ in km). The current default corresponds to a 95% circular home-range size of `",
                format_fixed(default_hr95_km2), " km²` at the prior mean of log(σ), with an approximate 95% interval of `",
                format_fixed(default_hr95_interval_km2[1]), "` to `", format_fixed(default_hr95_interval_km2[2]),
                " km²` (Christensen et al. in prep)."
              ),
              paste0(
                "- **log(λ₀)** ~ Normal(mean, SD): Prior on the expected count per camera-day when an animal's activity center is at the camera. The default mean of `",
                advanced_prior_defaults$log_lam0_mean, "` corresponds to an expected count of `",
                signif(exp(advanced_prior_defaults$log_lam0_mean), 3), "` per camera-day at the camera (`exp(",
                advanced_prior_defaults$log_lam0_mean, ") = ", signif(exp(advanced_prior_defaults$log_lam0_mean), 3), "`)."
              ),
              paste0(
                "- **State-space buffer**: Distance added around the camera array when defining the USCR state space. The buffer is set wide enough that animals with activity centers outside it contribute almost nothing to the counts (under ", round(100 * (1 - uscr_buffer_p), 1), "% even at the high end of the σ prior). The default is ",
                uscr_buffer_rule_text,
                ", and updates automatically if you change that prior. After fitting, the app also checks it against the posterior σ: if it is less than ", uscr_buffer_multiplier, " times the ", 100 * uscr_buffer_p, "% home-range radius at the largest posterior σ, the app grows it to ", uscr_buffer_safety_factor, " times that minimum and reruns."
              ),
              sep = "\n"
            )),
            fluidRow(
              column(
                6,
                numericInput(
                  "log_sigma_mean", "log(σ) mean (σ in km)",
                  value = advanced_prior_defaults$log_sigma_mean, step = 0.1
                )
              ),
              column(
                6,
                numericInput(
                  "log_sigma_sd", "log(σ) SD",
                  value = advanced_prior_defaults$log_sigma_sd, min = 0.01, step = 0.05
                )
              )
            ),
            fluidRow(
              column(
                6,
                numericInput(
                  "log_lam0_mean", "log(λ₀) mean",
                  value = advanced_prior_defaults$log_lam0_mean, step = 0.1
                )
              ),
              column(
                6,
                numericInput(
                  "log_lam0_sd", "log(λ₀) SD",
                  value = advanced_prior_defaults$log_lam0_sd, min = 0.1, step = 0.1
                )
              )
            ),
            numericInput(
              "uscr_buffer_m", "USCR state-space buffer (m)",
              value = default_uscr_buffer_m, min = 50, step = 50
            )
          )
        )
      )
    ),

    # ---------------------- USCR MODEL TAB ------------------------
    nav_panel(
      "USCR model",
      HTML(paste0('
        <h2>Model 1 — USCR (Unmarked Spatial Capture-Recapture)</h2>
        <p><strong>The gist:</strong> Estimates <strong>animal density</strong> from <strong>where and how many</strong> animals are counted <strong>across a camera array</strong>.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Per-camera observation model:</p>
            <p>$$y_j \\sim \\mathrm{Poisson}(\\mu_j), \\qquad \\log \\mu_j = \\log(\\mathrm{days}_j) + \\log(\\Lambda_j) + \\epsilon_j$$</p>
            <p>Expected count per camera-day from the augmented individuals:</p>
            <p>$$\\Lambda_j = \\sum_{i=1}^{M} z_i\\,\\lambda_0\\,\\exp\\!\\Big(-\\frac{d_{ij}^2}{2\\sigma^2}\\Big), \\qquad z_i \\sim \\mathrm{Bernoulli}(\\psi), \\qquad \\psi \\sim \\mathrm{Beta}(a, b)$$</p>
            <p>Activity centers and camera effects:</p>
            <p>$$s_i \\sim \\mathrm{Uniform}(S), \\qquad \\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>Here \\(S\\) is the buffered state space: the union of circles, with radius equal to the state-space buffer, around each camera. Abundance is \\(N = \\sum_i z_i\\), and density is \\(N\\) divided by the area of \\(S\\), reported in animals/km². The state-space buffer can be adjusted in <strong>Model settings &gt; Advanced &gt; Priors: USCR</strong>. Throughout, \\(\\mathcal{N}(m, s)\\) is a normal distribution with mean \\(m\\) and standard deviation \\(s\\).</p>
            <h3>Variables</h3>
            <ul>
              <li>\\(y_j\\) — number of animals counted at camera \\(j\\): the sum, over detection events in the camera\'s analysis window, of each event\'s maximum Sighting Count.</li>
              <li>\\(\\mathrm{days}_j\\) — camera-days for camera \\(j\\): the length of its analysis window, minus time within detection events.</li>
              <li>\\(\\mu_j\\) — expected number of animals counted at camera \\(j\\).</li>
              <li>\\(\\Lambda_j\\) — expected number of animals counted per camera-day at camera \\(j\\), summed over the individuals in the population.</li>
              <li>\\(\\lambda_0\\) — expected count of an individual per camera-day if its activity center were at the camera.</li>
              <li>\\(\\sigma\\) — spatial scale parameter (km), related to home-range size and how quickly detection falls with distance.</li>
              <li>\\(s_i\\) — activity center of augmented individual \\(i\\), located within the state space \\(S\\).</li>
              <li>\\(d_{ij}\\) — projected distance (km) from augmented individual \\(i\\)\'s activity center to camera \\(j\\).</li>
              <li>\\(M\\) — number of augmented individuals, the upper limit of the prior distribution for abundance \\(N\\).</li>
              <li>\\(z_i\\) — indicator that augmented individual \\(i\\) is part of the real population.</li>
              <li>\\(\\psi\\) — probability that an augmented individual is part of the real population. Its prior sets the prior on abundance and density.</li>
              <li>\\(\\epsilon_j\\) — camera-level random effect (log scale) for camera \\(j\\), allowing counts to vary more between cameras than a Poisson distribution allows.</li>
              <li>\\(\\sigma_\\epsilon\\) — standard deviation of the camera-level random effects.</li>
            </ul>
            <h3>Default priors used in the app (can be changed in Model settings)</h3>
            <ul>
              <li>\\(M\\), the upper bound on total abundance used for data augmentation. It starts at ', mcmc_defaults$M_uscr, ' (or about twice the upper bound on abundance implied by an informative density prior). If the posterior for \\(N\\) presses against \\(M\\), the app doubles it and reruns the model.</li>
              <li>\\(\\psi \\sim \\mathrm{Beta}(a, b)\\). The default \\(a = b = 1\\) (\\(\\psi \\sim \\mathrm{Uniform}(0, 1)\\)) gives a uniform prior on \\(N\\) from 0 to \\(M\\). An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) chooses \\(a\\) and \\(b\\) so that density \\(N / \\mathrm{area}\\) has the requested 95% prior interval.</li>
              <li>\\(\\log \\sigma \\sim \\mathcal{N}(', advanced_prior_defaults$log_sigma_mean, ',\\,', advanced_prior_defaults$log_sigma_sd, ')\\), which corresponds to a 95% circular home-range size of about ', format(round(default_hr95_km2, 2)), ' km² at the prior mean of \\(\\log \\sigma\\) (Christensen et al. in prep).</li>
              <li>\\(\\log \\lambda_0 \\sim \\mathcal{N}(', advanced_prior_defaults$log_lam0_mean, ',\\,', advanced_prior_defaults$log_lam0_sd, ')\\), centered at \\(\\lambda_0 = ', format(signif(exp(advanced_prior_defaults$log_lam0_mean), 3)), '\\) expected count per camera-day when the activity center is at the camera.</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(', advanced_prior_defaults$sd_eps_shape, ',\\,', advanced_prior_defaults$sd_eps_rate, ')\\), giving the camera-level random-effect standard deviation a prior mean of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate, 3)), ' and variance of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate^2, 3)), '. The random effect \\(\\epsilon_j\\) is included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise (all models)</strong> (off by default).</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
            </ul>
          </div>
        </details>
      ')),
      hr(),
      h3("Uploaded data"),
      uiOutput("run_note_uscr_nps"),
      div(
        actionButton("run_uscr_nps", "Run USCR on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      uiOutput("uscr_nps_result"),
      h4("Run status & troubleshooting"),
      verbatimTextOutput("uscr_nps_debug")
    ),
    # ---------------------- REM MODEL TAB -------------------------
    nav_panel(
      "REM model",
      HTML(paste0('
        <h2>Model 2 — REM (Random Encounter Model)</h2>
        <p><strong>The gist:</strong> Converts <strong>how often animals pass a camera</strong> into density, correcting for <strong>movement rate</strong> and <strong>camera view geometry</strong>.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Per-camera observation model:</p>
            <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j = \\log D + \\log(\\mathrm{days}_j) + \\log v + \\log r_j + \\log\\!\\Big(\\frac{2 + \\theta_{\\mathrm{rad},j}}{\\pi}\\Big) + \\epsilon_j$$</p>
            <p>Density and camera effects:</p>
            <p>$$D = D_{\\max}\\,p_D, \\qquad p_D \\sim \\mathrm{Beta}(a, b), \\qquad \\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>Density is reported in animals/km². Throughout, \\(\\mathcal{N}(m, s)\\) is a normal distribution with mean \\(m\\) and standard deviation \\(s\\).</p>
            <h3>Variables</h3>
            <ul>
              <li>\\(y_j\\) — number of animals counted at camera \\(j\\): the sum, over detection events in the camera\'s analysis window, of each event\'s maximum Sighting Count.</li>
              <li>\\(\\mathrm{days}_j\\) — camera-days for camera \\(j\\): the length of its analysis window, minus time within detection events.</li>
              <li>\\(\\lambda_j\\) — expected number of animals counted at camera \\(j\\).</li>
              <li>\\(D\\) — animal density (animals/km²).</li>
              <li>\\(v\\) — daily movement rate (km/day).</li>
              <li>\\(r_j\\) — effective detection radius (km) for camera \\(j\\): its Detection Distance.</li>
              <li>\\(\\theta_{\\mathrm{rad},j}\\) — full detection angle (radians) for camera \\(j\\): its Camera Detection Angle, or the default detection angle in Model settings if that value is missing or outside the accepted range (from ', camera_angle_min_deg, ' to ', camera_angle_max_deg, ' degrees).</li>
              <li>\\(D_{\\max}\\) — upper limit of the prior distribution for density \\(D\\).</li>
              <li>\\(p_D\\) — density as a proportion of \\(D_{\\max}\\). Its prior sets the prior on density.</li>
              <li>\\(\\epsilon_j\\) — camera-level random effect (log scale) for camera \\(j\\), allowing counts to vary more between cameras than a Poisson distribution allows.</li>
              <li>\\(\\sigma_\\epsilon\\) — standard deviation of the camera-level random effects.</li>
            </ul>
            <h3>Default priors used in the app (can be changed in Model settings)</h3>
            <ul>
              <li>\\(D_{\\max}\\), the upper limit on density. It starts at ', rem_tte_start_D_max, ' animals/km² (or twice the upper bound of an informative density prior). If the posterior for \\(D\\) presses against \\(D_{\\max}\\), the app doubles it and reruns the model.</li>
              <li>\\(p_D \\sim \\mathrm{Beta}(a, b)\\). The default \\(a = b = 1\\) (\\(p_D \\sim \\mathrm{Uniform}(0, 1)\\)) gives a uniform prior on \\(D\\) from 0 to \\(D_{\\max}\\). An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) chooses \\(a\\) and \\(b\\) so that \\(D\\) has the requested 95% prior interval.</li>
              <li>\\(\\log v \\sim \\mathcal{N}(', advanced_prior_defaults$log_v_mean, ',\\,', advanced_prior_defaults$log_v_sd, ')\\), which corresponds to a daily movement rate of about ', format_fixed(default_v_km_day), ' km/day at the prior mean of \\(\\log v\\) (Christensen et al. in prep).</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(', advanced_prior_defaults$sd_eps_shape, ',\\,', advanced_prior_defaults$sd_eps_rate, ')\\), giving the camera-level random-effect standard deviation a prior mean of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate, 3)), ' and variance of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate^2, 3)), '. The random effect \\(\\epsilon_j\\) is included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise (all models)</strong> (off by default).</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
            </ul>
          </div>
        </details>
      ')),
      hr(),
      h3("Uploaded data"),
      uiOutput("run_note_rem_nps"),
      div(
        actionButton("run_rem_nps", "Run REM on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      uiOutput("rem_nps_result"),
      h4("Run status & troubleshooting"),
      verbatimTextOutput("rem_nps_debug")
    ),

    # ---------------------- TTE MODEL TAB -------------------------
    nav_panel(
      "TTE model",
      HTML(paste0('
        <h2>Model 3 — TTE (Time-to-Event)</h2>
        <p><strong>The gist:</strong> Converts <strong>how many animals each camera counts</strong> into density, using <strong>movement rate</strong> to measure camera time in <strong>viewshed-crossing time units</strong> and correcting for <strong>viewshed area</strong>. Shorter waits between detections imply higher density.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Per-camera observation model:</p>
            <p>$$y_j \\sim \\mathrm{Poisson}(\\lambda_j), \\qquad \\log \\lambda_j = \\log D + \\log(\\mathrm{days}_j) + \\log U_j + \\log A_j + \\epsilon_j$$</p>
            <p>Time units and viewshed area:</p>
            <p>$$U_j = \\frac{v}{f(\\theta_{\\mathrm{deg},j})\\,r_j}, \\qquad A_j = \\pi r_j^2 \\frac{\\theta_{\\mathrm{deg},j}}{360}$$</p>
            <p>$$f(\\theta_{\\mathrm{deg},j}) = 0.3324 + 0.005580\\,\\theta_{\\mathrm{deg},j} - 1.454\\times10^{-5}\\,\\theta_{\\mathrm{deg},j}^2$$</p>
            <p>Density and camera effects:</p>
            <p>$$D = D_{\\max}\\,p_D, \\qquad p_D \\sim \\mathrm{Beta}(a, b), \\qquad \\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)$$</p>
            <p>\\(f(\\theta_{\\mathrm{deg},j})\\,r_j\\) is the average length of a straight-line path across camera \\(j\\)\'s viewshed, a circular sector with radius \\(r_j\\) and full angle \\(\\theta_{\\mathrm{deg},j}\\); \\(f\\) was fitted to simulations of random straight-line paths across circular sectors. Density is reported in animals/km². Throughout, \\(\\mathcal{N}(m, s)\\) is a normal distribution with mean \\(m\\) and standard deviation \\(s\\).</p>
            <h3>Variables</h3>
            <ul>
              <li>\\(y_j\\) — number of animals counted at camera \\(j\\): the sum, over detection events in the camera\'s analysis window, of each event\'s maximum Sighting Count.</li>
              <li>\\(\\mathrm{days}_j\\) — camera-days for camera \\(j\\): the length of its analysis window, minus time within detection events.</li>
              <li>\\(\\lambda_j\\) — expected number of animals counted at camera \\(j\\).</li>
              <li>\\(D\\) — animal density (animals/km²).</li>
              <li>\\(v\\) — daily movement rate (km/day).</li>
              <li>\\(r_j\\) — effective detection radius (km) for camera \\(j\\): its Detection Distance.</li>
              <li>\\(\\theta_{\\mathrm{deg},j}\\) — full detection angle (degrees) for camera \\(j\\): its Camera Detection Angle, or the default detection angle in Model settings if that value is missing or outside the accepted range (from ', camera_angle_min_deg, ' to ', camera_angle_max_deg, ' degrees).</li>
              <li>\\(U_j\\) — number of time units per camera-day at camera \\(j\\), where a time unit is the time an animal with movement rate \\(v\\) takes to cross the camera\'s viewshed (average path length \\(f(\\theta_{\\mathrm{deg},j})\\,r_j\\)).</li>
              <li>\\(A_j\\) — viewshed area (km²) for camera \\(j\\).</li>
              <li>\\(D_{\\max}\\) — upper limit of the prior distribution for density \\(D\\).</li>
              <li>\\(p_D\\) — density as a proportion of \\(D_{\\max}\\). Its prior sets the prior on density.</li>
              <li>\\(\\epsilon_j\\) — camera-level random effect (log scale) for camera \\(j\\), allowing counts to vary more between cameras than a Poisson distribution allows.</li>
              <li>\\(\\sigma_\\epsilon\\) — standard deviation of the camera-level random effects.</li>
            </ul>
            <h3>Default priors used in the app (can be changed in Model settings)</h3>
            <ul>
              <li>\\(D_{\\max}\\), the upper limit on density. It starts at ', rem_tte_start_D_max, ' animals/km² (or twice the upper bound of an informative density prior). If the posterior for \\(D\\) presses against \\(D_{\\max}\\), the app doubles it and reruns the model.</li>
              <li>\\(p_D \\sim \\mathrm{Beta}(a, b)\\). The default \\(a = b = 1\\) (\\(p_D \\sim \\mathrm{Uniform}(0, 1)\\)) gives a uniform prior on \\(D\\) from 0 to \\(D_{\\max}\\). An informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) chooses \\(a\\) and \\(b\\) so that \\(D\\) has the requested 95% prior interval.</li>
              <li>\\(\\log v \\sim \\mathcal{N}(', advanced_prior_defaults$log_v_mean, ',\\,', advanced_prior_defaults$log_v_sd, ')\\), which corresponds to a daily movement rate of about ', format_fixed(default_v_km_day), ' km/day at the prior mean of \\(\\log v\\) (Christensen et al. in prep).</li>
              <li>\\(\\sigma_\\epsilon \\sim \\mathrm{Gamma}(', advanced_prior_defaults$sd_eps_shape, ',\\,', advanced_prior_defaults$sd_eps_rate, ')\\), giving the camera-level random-effect standard deviation a prior mean of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate, 3)), ' and variance of ', format(signif(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate^2, 3)), '. The random effect \\(\\epsilon_j\\) is included only when camera-level noise is turned on in <strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise (all models)</strong> (off by default).</li>
              <li><strong>Note:</strong> these defaults are not well calibrated for species other than white-tailed deer (<em>Odocoileus virginianus</em>), and model output is sensitive to choice of prior distribution.</li>
            </ul>
          </div>
        </details>
      ')),
      hr(),
      h3("Uploaded data"),
      uiOutput("run_note_tte_nps"),
      div(
        actionButton("run_tte_nps", "Run TTE on uploaded data", class = "btn-primary"),
        style = "margin-bottom: 10px;"
      ),
      br(),
      uiOutput("tte_nps_result"),
      h4("Run status & troubleshooting"),
      verbatimTextOutput("tte_nps_debug")
    ),

    # ---------------------- COMPARE & COMBINE --------------------
    nav_panel(
      "Compare & combine",
      HTML(paste0('
        <h2>Compare &amp; combine</h2>
        <p><strong>The gist:</strong> Compares the models run on <strong>uploaded data</strong> using <strong>WAIC</strong>, and combines their density estimates into <strong>model averages</strong>.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (equations)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>WAIC weights, where \\(\\mathrm{WAIC}_k\\) is model \\(k\\)\'s WAIC (averaged across its chains):</p>
            <p>$$w_k = \\frac{\\exp(-\\Delta_k / 2)}{\\sum_j \\exp(-\\Delta_j / 2)}, \\qquad \\Delta_k = \\mathrm{WAIC}_k - \\min_j \\mathrm{WAIC}_j$$</p>
            <p>Quantile averages, where \\(D_{k(s)}\\) is the \\(s\\)th smallest posterior draw of density from model \\(k\\):</p>
            <p>$$\\bar{D}_{(s)} = \\sum_k w_k\\,D_{k(s)}$$</p>
            <p>The WAIC-weighted quantile average uses the WAIC weights \\(w_k\\) above; the unweighted quantile average sets each \\(w_k = 1/K\\), where \\(K\\) is the number of models.</p>
            <p>WAIC-weighted mixture, where \\(p_k(D \\mid y)\\) is model \\(k\\)\'s posterior distribution of density:</p>
            <p>$$p_{\\mathrm{mix}}(D \\mid y) = \\sum_k w_k\\,p_k(D \\mid y)$$</p>
            <p>Each mixture draw comes from model \\(k\\) with probability \\(w_k\\), as a random posterior draw of density from that model.</p>
            <p>For each model and model average, \\(P(D > \\text{threshold})\\) is the proportion of draws of density above the threshold set below. Density is reported in animals/km².</p>
            <h3>Notes</h3>
            <ul>
              <li>The table summarizes whichever models have been run on uploaded data, and updates as each one finishes. Model averages need at least two models. The WAIC-weighted averages also need a WAIC for every model; otherwise only the unweighted quantile average is shown.</li>
              <li>A Bayesian p-value in red means that model failed its posterior predictive check (p-value below ', ppc_pass_range[1], ' or above ', ppc_pass_range[2], '), so its estimates and intervals may be unreliable.</li>
              <li>Results are cleared whenever the data the models were fit to change, so models are only compared and combined when fit to the same data.</li>
              <li>WAIC isn\'t computed for simulated data, so this tab covers uploaded data only.</li>
            </ul>
          </div>
        </details>
      ')),
      hr(),
      h3("Uploaded data"),
      numericInput(
        "combo_density_threshold",
        "Density threshold (animals/km²)",
        value = default_density_threshold,
        min = 0,
        step = 1
      ),
      DTOutput("nps_combo_table"),
      h4("Visual summaries"),
      plotOutput("nps_combo_interval_plot", height = "320px"),
      div(
        downloadButton("dl_nps_combo_interval_png", "Download interval plot (PNG)"),
        style = "margin-bottom: 12px;"
      ),
      plotOutput("nps_combo_density_plot", height = "320px"),
      div(
        downloadButton("dl_nps_combo_density_png", "Download posterior density plot (PNG)"),
        style = "margin-bottom: 12px;"
      ),
      h4("Posterior summaries"),
      downloadButton("dl_nps_all_csv", "Download posterior summaries (CSV)"),
      posterior_param_legend()
    ),

    # ---------------------- SIMULATE DATA -------------------------
    nav_panel(
      "Simulate data",
      HTML(paste0('
        <h2>Simulate data</h2>
        <p><strong>The gist:</strong> Simulates camera data from a <strong>study design you choose</strong>, then fits each model to see how <strong>precisely and accurately</strong> it recovers the <strong>true density</strong>.</p>
        <p>Choose a camera grid, spacing, survey length, and true density, plus the values each model needs, then click <strong>Simulate data</strong>. Then run each model in its own section below.</p>
        <details>
          <summary style="cursor: pointer; font-weight: 600; margin: 1rem 0; padding: 0.5rem; background: #f0f4e8; border-left: 3px solid #609048; border-radius: 6px;"><strong>Under the hood (how the simulations work)</strong></summary>
          <div style="margin: 1rem 0; padding-left: 1rem;">
            <p>Each model gets its own simulated dataset, generated from that model\'s own observation model (see its tab):</p>
            <ul>
              <li><strong>USCR:</strong> activity centers are placed uniformly at random within the simulation region (the union of circles, with radius equal to the simulation buffer, around each camera), and their number is Poisson, with mean equal to the true density times the region\'s area. Each camera\'s count is Poisson, with a mean that depends on how far each activity center is from the camera, on \\(\\sigma\\), and on \\(\\lambda_0\\). The fit uses its own state space, set by the state-space buffer in <strong>Model settings</strong>, as for uploaded data, and activity centers outside it still add to the counts, as in a field study.</li>
              <li><strong>REM and TTE:</strong> each camera\'s count is Poisson, with the mean the model expects from the true density, daily movement rate, detection radius, and detection angle.</li>
            </ul>
            <p>All three datasets use the same camera grid, survey length, and true density. When camera-level noise is on (<strong>Model settings &gt; Advanced &gt; Priors: Camera-level noise (all models)</strong>), each camera\'s mean is also multiplied by \\(\\exp(\\epsilon_j)\\), with \\(\\epsilon_j \\sim \\mathcal{N}(0, \\sigma_\\epsilon)\\), so the counts are Poisson-lognormal. Here \\(\\mathcal{N}(m, s)\\) is a normal distribution with mean \\(m\\) and standard deviation \\(s\\).</p>
            <h3>Notes</h3>
            <ul>
              <li>Each dataset is simulated from the same model that is fit to it, so results show how well a design can do when the model\'s assumptions hold.</li>
              <li>Every prior, including the prior means for \\(\\log \\sigma\\), \\(\\log \\lambda_0\\), and \\(\\log v\\), comes from <strong>Model settings</strong>, as for uploaded data, along with the USCR state-space buffer and the maximum number of adaptive rounds. The note next to each Run button lists any prior means that differ from the values the data were simulated with.</li>
              <li>You can also see how a fit behaves when its priors are far from the simulated truth: for example, a \\(\\log \\sigma\\) or \\(\\log v\\) prior mean far from the log of the true value, or an informative density prior (<strong>Model settings &gt; Advanced &gt; Priors: Density (all models)</strong>) whose 95% interval sits well away from the true density.</li>
              <li>Each click of <strong>Simulate data</strong> simulates one dataset per model, so a single hit or miss isn\'t conclusive. Try several random seeds before drawing conclusions. If a model often misses the truth under a design, that design may need more cameras, a longer survey, or different camera spacing.</li>
              <li>WAIC isn\'t computed for simulated data, so these fits don\'t appear on the <strong>Compare &amp; combine</strong> tab.</li>
            </ul>
          </div>
        </details>
      ')),
      hr(),
      h2("Simulation parameters"),
      tags$p(
        class = "small",
        style = "color: var(--muted);",
        "The default values match the app's default prior means and default detection angle, so with the default Model settings the simulated truth sits at the center of each model's priors."
      ),
      fluidRow(
        column(
          4,
          h3("All models"),
          sliderInput(
            "pshared_n_side", "Grid dimension (n × n cameras)",
            min = 3, max = 8, value = 5, step = 1
          ),
          sliderInput(
            "pshared_spacing", "Camera spacing (m)",
            min = 150, max = 500, value = 300, step = 10
          ),
          sliderInput(
            "pshared_days", "Survey length (days)",
            min = 7, max = 60, value = 21, step = 1
          ),
          sliderInput(
            "pshared_Dtrue", "True density (animals/km²)",
            min = 5, max = 80, value = 25, step = 1
          ),
          conditionalPanel(
            "input.camera_noise && input.mode == 'Advanced'",
            sliderInput(
              "pshared_sd_eps", HTML("Camera-level noise SD (σ<sub>ε</sub>)"),
              min = 0.05, max = 2, step = 0.01,
              value = round(advanced_prior_defaults$sd_eps_shape / advanced_prior_defaults$sd_eps_rate, 2)
            )
          ),
          conditionalPanel(
            "!(input.camera_noise && input.mode == 'Advanced')",
            tags$p(
              class = "small",
              style = "color: var(--muted);",
              "Camera-level noise is off, so simulated counts are Poisson. Turn it on in ",
              tags$strong("Model settings > Advanced > Priors: Camera-level noise (all models)"),
              " to simulate and fit Poisson-lognormal counts."
            )
          )
        ),
        column(
          4,
          h3("USCR inputs"),
          sliderInput(
            "pshared_home_range_km2", "95% home-range size (km²)",
            min = 0.10, max = 3.00, step = 0.01,
            value = round(default_hr95_km2, 2)
          ),
          sliderInput(
            "pshared_lambda0", "Detection rate λ₀ (expected count per camera-day at the camera)",
            min = 0.05, max = 2, step = 0.01,
            value = round(exp(advanced_prior_defaults$log_lam0_mean), 2)
          ),
          numericInput(
            "pshared_buffer_m", "Simulation buffer (m)",
            value = default_uscr_buffer_m, min = 50, step = 50
          ),
          tags$p(
            class = "small",
            style = "color: var(--muted);",
            paste0(
              "Simulated activity centers are placed within this distance of at least one camera. By default it is the larger of the USCR state-space buffer the fit starts with (Model settings) and ",
              uscr_buffer_multiplier, " times the ", 100 * uscr_buffer_p,
              "% home-range radius at the true σ, rounded up to the nearest 50 m. It updates when either changes, and you can override it."
            )
          )
        ),
        column(
          4,
          h3("REM and TTE inputs"),
          sliderInput(
            "pshared_v", "Daily movement rate v (km/day)",
            min = 0.5, max = 10, step = 0.1,
            value = round(default_v_km_day, 1)
          ),
          sliderInput(
            "pshared_theta", "Detection angle θ (degrees)",
            min = camera_angle_min_deg, max = camera_angle_max_deg, value = default_theta_deg, step = 1
          ),
          sliderInput(
            "pshared_r", "Detection radius (m)",
            min = 8, max = 25, value = 12, step = 1
          )
        )
      ),
      fluidRow(
        column(4, numericInput("pshared_seed", "Random seed", value = 1, min = 1, step = 1)),
        column(8, br(), actionButton("pshared_run_sim", "Simulate data", class = "btn-primary"))
      ),
      tags$p(
        class = "small",
        style = "color: var(--muted);",
        "USCR, REM, and TTE each get their own dataset, simulated from that model's own observation model with the design and true density above. The same seed and settings always give the same datasets."
      ),
      tags$details(
        class = "about-card",
        style = "margin-top: 1rem;",
        tags$summary(
          style = "cursor: pointer; font-weight: 600;",
          "Simulation parameter guide"
        ),
        tags$div(style = "margin-top: 0.85rem;"),
        tags$h4("All models"),
        tags$p(
          tags$strong("Grid dimension"),
          " sets the number of cameras along one side of the square grid. A 5 × 5 grid gives 25 cameras."
        ),
        tags$p(
          tags$strong("Camera spacing"),
          " is the distance between neighboring cameras. For USCR, cameras should be close enough that an individual can be detected at more than one camera."
        ),
        tags$p(
          tags$strong("Survey length"),
          " is the number of camera-days simulated for each camera."
        ),
        tags$p(
          tags$strong("True density"),
          " is the animal density used to generate the data, in animals/km². Estimates are reported in the same units."
        ),
        tags$p(
          tags$strong(HTML("Camera-level noise SD (σ<sub>ε</sub>)")),
          " is the standard deviation of camera-to-camera variation in counts (on the log scale) that is not explained by survey length, density, movement, or camera geometry. It is used only when camera-level noise is turned on in ",
          tags$strong("Model settings > Advanced > Priors: Camera-level noise (all models)", .noWS = "after"),
          "; then it is applied when simulating data for all three models, and all three models estimate ",
          HTML("σ<sub>ε</sub>"),
          " when fit."
        ),
        tags$p(
          tags$strong("Random seed"),
          " reproduces the same simulated datasets when you rerun with identical settings. Change it to see how results vary between datasets."
        ),
        tags$h4("USCR inputs"),
        tags$p(
          tags$strong("95% home-range size"),
          " is the area of a circular home range containing 95% of an animal's activity. It is converted to the spatial scale parameter σ, which every simulated animal shares. Larger home ranges spread each animal's detections across more cameras."
        ),
        tags$p(
          tags$strong("Detection rate λ₀"),
          " is the expected count of an individual per camera-day if its activity center were at the camera. Larger values produce larger counts everywhere in the array."
        ),
        tags$p(
          tags$strong("Simulation buffer"),
          paste0(
            " is the maximum distance from the nearest camera at which simulated activity centers are placed, so the simulated region is the union of circles of this radius around the cameras. The USCR fit uses its own state space, set by the state-space buffer in Model settings (by default ",
            uscr_buffer_rule_text,
            "), and grows it if it is less than ",
            uscr_buffer_multiplier, " times the ", 100 * uscr_buffer_p,
            "% home-range radius at the largest posterior σ. By default the simulation buffer is the larger of the fit's starting buffer and ",
            uscr_buffer_multiplier, " times the ", 100 * uscr_buffer_p,
            "% home-range radius at the true σ, so the simulated region covers the fit's state space and the cameras' surroundings. If the fit's buffer is ever larger than the simulation buffer, part of its state space has no simulated activity centers and density will likely be underestimated; the run note and run-status panel warn when this happens."
          )
        ),
        tags$h4("REM and TTE inputs"),
        tags$p(
          tags$strong("Daily movement rate v"),
          " is the average distance an animal moves per day, in km/day. It is species-specific."
        ),
        tags$p(
          tags$strong("Detection angle θ"),
          paste0(
            " is the camera's detection angle in degrees. Manufacturer specifications typically list this value. As for uploaded cameras, it is limited to ",
            camera_angle_min_deg, "–", camera_angle_max_deg, "°."
          )
        ),
        tags$p(
          tags$strong("Detection radius"),
          " is how far a camera can detect an animal, in meters."
        )
      ),
      hr(),
      h2("Fit each model to the simulated data"),
      markdown(paste(
        "**Reading the results tables:**",
        "",
        "- **Truth**: the value used to simulate the data.",
        "- **Mean / Median**: posterior mean and median.",
        "- **Lower 2.5% / Upper 97.5%**: bounds of the 95% credible interval.",
        "- **Covered**: whether the 95% credible interval contains the true value. Across many simulated datasets, a well-calibrated model should cover the truth about 95% of the time, so a single \"No\" is not by itself a problem.",
        paste0(
          "- The line under each table gives the model's posterior predictive check. It is red when the check fails (Bayesian p-value below ",
          ppc_pass_range[1], " or above ", ppc_pass_range[2], ")."
        ),
        sep = "\n"
      )),
      uiOutput("pshared_noise_note"),

      # ---------------- USCR block ----------------
      tags$div(
        class = "about-card",
        h3("USCR"),
        tags$p(
          class = "small",
          style = "color: var(--muted);",
          "Click the drop-downs below to see a summary of the simulated data: the camera array with the simulated activity centers, the counts at each camera, and a table of the cameras."
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated camera array and counts")),
          plotOutput("pshared_uscr_grid_plot", height = "380px"),
          plotOutput("pshared_uscr_count_plot", height = "380px")
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated camera data")),
          DTOutput("pshared_uscr_camera_table")
        ),
        uiOutput("run_note_uscr_pshared"),
        actionButton("pshared_run_uscr", "Run USCR", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_uscr_table"),
        h4("Run status & troubleshooting"),
        verbatimTextOutput("pshared_uscr_debug")
      ),

      # ---------------- REM block ----------------
      tags$div(
        class = "about-card",
        h3("REM"),
        tags$p(
          class = "small",
          style = "color: var(--muted);",
          "Click the drop-downs below to see a summary of the simulated data: the counts at each camera, and a table of the cameras."
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated counts")),
          plotOutput("pshared_rem_count_plot", height = "380px")
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated camera data")),
          DTOutput("pshared_rem_camera_table")
        ),
        uiOutput("run_note_rem_pshared"),
        actionButton("pshared_run_rem", "Run REM", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_rem_table"),
        h4("Run status & troubleshooting"),
        verbatimTextOutput("pshared_rem_debug")
      ),

      # ---------------- TTE block ----------------
      tags$div(
        class = "about-card",
        h3("TTE"),
        tags$p(
          class = "small",
          style = "color: var(--muted);",
          "Click the drop-downs below to see a summary of the simulated data: the counts at each camera, and a table of the cameras."
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated counts")),
          plotOutput("pshared_tte_count_plot", height = "380px")
        ),
        tags$details(
          class = "data-details",
          tags$summary(tags$strong("Simulated camera data")),
          DTOutput("pshared_tte_camera_table")
        ),
        uiOutput("run_note_tte_pshared"),
        actionButton("pshared_run_tte", "Run TTE", class = "btn-primary"),
        br(), br(),
        DTOutput("pshared_tte_table"),
        h4("Run status & troubleshooting"),
        verbatimTextOutput("pshared_tte_debug")
      )
    )
  )
)

# -------------------------------------------------------------------
# SERVER
# -------------------------------------------------------------------

server <- function(input, output, session) {

  # ================================================================
  # SIMULATE DATA: simulation and model fits (USCR + REM + TTE)
  # ================================================================

  pshared_sim <- reactiveVal(NULL)
  pshared_uscr_fit <- reactiveVal(NULL)
  pshared_rem_fit <- reactiveVal(NULL)
  pshared_tte_fit <- reactiveVal(NULL)
  pshared_uscr_debug <- reactiveVal(make_model_debug("USCR", "simulated data"))
  pshared_rem_debug <- reactiveVal(make_model_debug("REM", "simulated data"))
  pshared_tte_debug <- reactiveVal(make_model_debug("TTE", "simulated data"))

  # The simulation's sigma (km), from the 95% home-range size.
  pshared_sigma_km <- reactive({
    req(is.numeric(input$pshared_home_range_km2))
    uscr_sigma_km_from_home_range(input$pshared_home_range_km2, p = 0.95)
  })

  # Automatic simulation buffer (m): the larger of the USCR state-space
  # buffer the fit starts with (Model settings) and the post-fit buffer rule
  # at the true sigma, so the simulated activity centers cover the fit's
  # state space and every camera's neighborhood at the true sigma.
  pshared_default_buffer_m <- reactive({
    fit_buffer_m <- app_uscr_buffer_m()
    req(is_valid_buffer_m(fit_buffer_m))
    round_buffer_m(max(fit_buffer_m, uscr_buffer_from_sigma_km(pshared_sigma_km())))
  })

  # Keep the simulation buffer at the automatic value when the home-range
  # size or the USCR state-space buffer changes. Users can override it.
  observe({
    updateNumericInput(session, "pshared_buffer_m", value = pshared_default_buffer_m())
  })

  observeEvent(input$pshared_run_sim, {
    seed <- input$pshared_seed
    buffer_m <- input$pshared_buffer_m
    valid_seed <- is.numeric(seed) && length(seed) == 1 && is.finite(seed) &&
      seed == round(seed) && seed >= 1 && seed < .Machine$integer.max - 2
    problem <- if (!valid_seed) {
      "The random seed must be a positive whole number."
    } else if (!is_valid_buffer_m(buffer_m)) {
      "The simulation buffer (USCR inputs) must be a positive number of meters."
    }
    if (!is.null(problem)) {
      showNotification(paste("Data not simulated.", problem), type = "error", duration = NULL)
      return()
    }

    # Settings are saved with the simulated data so that the model fits
    # (camera-level noise), the truth-vs-estimate tables, and the run notes
    # use the values the data were simulated with, even if an input changes
    # afterward without clicking 'Simulate data' again.
    settings <- list(
      n_side = input$pshared_n_side,
      spacing_m = input$pshared_spacing,
      days = input$pshared_days,
      D_true = input$pshared_Dtrue,
      sigma_km = pshared_sigma_km(),
      lambda0 = input$pshared_lambda0,
      buffer_m = buffer_m,
      v = input$pshared_v,
      theta = input$pshared_theta,
      r = input$pshared_r,
      sd_eps = input$pshared_sd_eps,
      camera_noise = app_camera_noise()
    )

    # Seed offsets give REM, TTE, and USCR separate random streams.
    simulate_encounter <- function(model, seed) {
      simulate_encounter_counts(
        model = model,
        n_side = settings$n_side,
        spacing_m = settings$spacing_m,
        days = settings$days,
        D_per_km2 = settings$D_true,
        detection_radius_m = settings$r,
        theta_deg = settings$theta,
        v_km_day = settings$v,
        sd_eps = settings$sd_eps,
        camera_noise = settings$camera_noise,
        seed = seed
      )
    }
    pshared_sim(list(
      uscr = simulate_camera_counts(
        n_side = settings$n_side,
        spacing_m = settings$spacing_m,
        days = settings$days,
        D_per_km2 = settings$D_true,
        lambda0 = settings$lambda0,
        sigma_km = settings$sigma_km,
        buffer_m = settings$buffer_m,
        sd_eps = settings$sd_eps,
        camera_noise = settings$camera_noise,
        seed = seed + 2
      ),
      rem = simulate_encounter("REM", seed),
      tte = simulate_encounter("TTE", seed + 1),
      settings = settings
    ))

    # Results and run status from the previous simulation no longer apply.
    pshared_uscr_fit(NULL)
    pshared_rem_fit(NULL)
    pshared_tte_fit(NULL)
    pshared_uscr_debug(make_model_debug("USCR", "simulated data"))
    pshared_rem_debug(make_model_debug("REM", "simulated data"))
    pshared_tte_debug(make_model_debug("TTE", "simulated data"))
  })

  # Which noise setting the current simulated dataset (and so its fits) uses.
  output$pshared_noise_note <- renderUI({
    s <- pshared_sim()
    req(s)
    sim_noise <- isTRUE(s$settings$camera_noise)
    toggle_noise <- app_camera_noise()
    tags$div(
      class = "about-card",
      style = "padding: 10px 14px;",
      tags$p(
        style = "margin: 0;",
        tags$strong(paste0("Current simulated data: camera-level noise ", if (sim_noise) "on" else "off", ". ")),
        "All three fits below use this setting. To change it, change the camera-level noise setting in ",
        tags$strong("Model settings > Advanced > Priors: Camera-level noise (all models)"),
        " and click ",
        tags$strong("Simulate data"),
        " again."
      ),
      if (sim_noise != toggle_noise) {
        tags$p(
          style = "margin: 0.4rem 0 0 0; color: firebrick;",
          paste0(
            "The Model settings toggle is now ", if (toggle_noise) "on" else "off",
            ". It will apply to the next simulation and to uploaded data, not to the current simulated data."
          )
        )
      }
    )
  })

  output$pshared_uscr_grid_plot <- renderPlot({
    s <- pshared_sim()
    req(s)
    sim <- s$uscr
    fit_buffer_m <- app_uscr_buffer_m()
    fit_outline <- if (is_valid_buffer_m(fit_buffer_m)) {
      uscr_buffered_outline(cbind(sim$out$utm_e, sim$out$utm_n), fit_buffer_m / 1000, scale = 1000)
    }
    ggplot() +
      geom_path(
        data = sim$state_space_outline, aes(x = x, y = y, group = group),
        color = redwood_colors[5], linetype = "dashed"
      ) +
      geom_path(
        data = fit_outline, aes(x = x, y = y, group = group),
        color = redwood_colors[1]
      ) +
      geom_point(
        data = sim$activity_centers, aes(x = x, y = y),
        size = 1, alpha = 0.5, color = redwood_colors[2]
      ) +
      geom_point(data = sim$out, aes(x = x, y = y), size = 3, color = redwood_colors[3]) +
      geom_label(
        data = sim$out, aes(x = x, y = y, label = `Site Name`),
        nudge_y = 18, size = 3.8, color = redwood_colors[1], fill = "white", linewidth = 0.2
      ) +
      coord_equal() +
      theme_minimal() +
      labs(
        x = "m", y = "m",
        title = sprintf(
          "%d cameras, %d m spacing (true density: %.1f animals/km²)",
          nrow(sim$out), round(s$settings$spacing_m), s$settings$D_true
        ),
        subtitle = paste0(
          sprintf(
            "%d simulated activity centers within the %d m simulation buffer (dashed outline, %.2f km²)",
            sim$N, round(s$settings$buffer_m), sim$area_km2
          ),
          if (!is.null(fit_outline)) {
            sprintf("\nUSCR fit's starting state space: %d m buffer (solid outline)", round(fit_buffer_m))
          }
        )
      )
  })

  output$pshared_uscr_count_plot <- renderPlot({
    s <- pshared_sim()
    req(s)
    sim_count_plot(s$uscr, "USCR", s$settings$D_true)
  })

  output$pshared_rem_count_plot <- renderPlot({
    s <- pshared_sim()
    req(s)
    sim_count_plot(s$rem, "REM", s$settings$D_true)
  })

  output$pshared_tte_count_plot <- renderPlot({
    s <- pshared_sim()
    req(s)
    sim_count_plot(s$tte, "TTE", s$settings$D_true)
  })

  output$pshared_uscr_camera_table <- renderDT({
    s <- pshared_sim()
    req(s)
    sim_camera_datatable(s$uscr)
  })

  output$pshared_rem_camera_table <- renderDT({
    s <- pshared_sim()
    req(s)
    sim_camera_datatable(s$rem)
  })

  output$pshared_tte_camera_table <- renderDT({
    s <- pshared_sim()
    req(s)
    sim_camera_datatable(s$tte)
  })

  # Checks shared by the simulated-data runs. Returns the simulation, or NULL
  # after a pop-up explaining why the model can't be run.
  sim_run_inputs <- function(model) {
    problem <- if (is.null(pshared_sim())) {
      "Click Simulate data first."
    } else {
      model_settings_problem(model)
    }
    if (!is.null(problem)) {
      showNotification(paste(model, "not run.", problem), type = "error", duration = NULL)
      return(NULL)
    }
    pshared_sim()
  }

  # Run a model on the simulated data. Camera-level noise comes from the
  # simulation; there is no WAIC for simulated data.
  run_on_simulated <- function(model) {
    s <- sim_run_inputs(model)
    if (is.null(s)) return()
    camera_noise <- s$settings$camera_noise

    if (identical(model, "USCR")) {
      buffer_m <- app_uscr_buffer_m()
      M <- uscr_start_M(s$uscr$out, buffer_m, density_prior_spec())
      context <- summarize_uscr_context(s$uscr, buffer_m, M, camera_noise)
      fit_model <- uscr_fitter(s$uscr, M, camera_noise, compute_WAIC = FALSE)
      extra_warning <- function(fit) {
        sim_buffer_warning(fit$settings$buffer_m, s$settings$buffer_m)
      }
    } else {
      d <- s[[tolower(model)]]
      angle_summary <- paste("Camera detection angle:", s$settings$theta, "degrees")
      context <- summarize_encounter_context(d, model, angle_summary, camera_noise)
      fit_model <- encounter_fitter(model, d, s$settings$theta, camera_noise, compute_WAIC = FALSE)
      extra_warning <- NULL
    }

    run_model(
      model = model,
      data_label = "simulated data",
      fit_rv = list(USCR = pshared_uscr_fit, REM = pshared_rem_fit, TTE = pshared_tte_fit)[[model]],
      debug_rv = list(USCR = pshared_uscr_debug, REM = pshared_rem_debug, TTE = pshared_tte_debug)[[model]],
      context = c(sim_truth_context(s, model), context),
      fit_model = fit_model,
      extra_warning = extra_warning
    )
  }

  observeEvent(input$pshared_run_uscr, run_on_simulated("USCR"))
  observeEvent(input$pshared_run_rem, run_on_simulated("REM"))
  observeEvent(input$pshared_run_tte, run_on_simulated("TTE"))

  # Truth-vs-estimate table for a model fit to the simulated data.
  render_truth_table <- function(model, fit_rv) {
    renderDT({
      fit <- fit_rv()
      # req() rather than validate(): a validation message inside an empty
      # DT widget has no height and overlaps the heading below it.
      req(fit)
      tbl <- summarize_vs_truth(fit, sim_truth_values(pshared_sim()$settings, model))
      results_datatable(tbl, fit)
    })
  }
  output$pshared_uscr_table <- render_truth_table("USCR", pshared_uscr_fit)
  output$pshared_rem_table <- render_truth_table("REM", pshared_rem_fit)
  output$pshared_tte_table <- render_truth_table("TTE", pshared_tte_fit)

  output$pshared_uscr_debug <- renderText({
    format_model_debug(pshared_uscr_debug(), pshared_uscr_fit())
  })

  output$pshared_rem_debug <- renderText({
    format_model_debug(pshared_rem_debug(), pshared_rem_fit())
  })

  output$pshared_tte_debug <- renderText({
    format_model_debug(pshared_tte_debug(), pshared_tte_fit())
  })

  # ================================================================
  # UPLOADED DATA: deployment + images QC (Add your data tab)
  # ================================================================

  # Deployment trim problem message, or NULL when the trim is off or valid.
  trim_days_problem <- function() {
    if (!isTRUE(input$apply_trim)) return(NULL)
    x <- input$trim_days
    if (is.numeric(x) && length(x) == 1 && is.finite(x) && x >= 1 && x == round(x)) {
      return(NULL)
    }
    "The number of days to keep must be a whole number, 1 or more."
  }

  # Deployment trim (days), or NULL when trimming is off. While the trim
  # entry is invalid, outputs that depend on it show the problem instead.
  app_trim_days <- function() {
    problem <- trim_days_problem()
    validate(need(is.null(problem), problem))
    if (isTRUE(input$apply_trim)) as.integer(input$trim_days)
  }

  output$trim_days_problem <- renderUI({
    problem <- trim_days_problem()
    if (!is.null(problem)) tags$p(class = "small", style = "color: firebrick;", problem)
  })

  # Cleaned (pre-QC) uploads are kept, with the messages from reading and
  # cleaning them, so the checks can be re-run when the hemisphere settings
  # change.
  deployment_raw <- reactiveVal(NULL)
  deployment_import_issues <- reactiveVal(NULL)
  images_raw <- reactiveVal(NULL)
  images_import_issues <- reactiveVal(NULL)

  deployment_checked <- reactiveVal(NULL)
  deployment_issues <- reactiveVal(NULL)
  images_checked <- reactiveVal(NULL)
  images_issues <- reactiveVal(NULL)

  run_deployment_check <- function() {
    if (is.null(deployment_raw())) {
      deployment_checked(NULL)
      deployment_issues(deployment_import_issues())
      return(invisible(NULL))
    }
    res <- run_qc_step(
      check_deployments(
        deployment_raw(),
        hemisphere = input$hemisphere,
        lat_hemisphere = input$lat_hemisphere
      )
    )
    deployment_checked(res$value)
    deployment_issues(combine_issues(deployment_import_issues(), res))
  }

  run_images_check <- function() {
    if (is.null(images_raw())) {
      images_checked(NULL)
      images_issues(images_import_issues())
      return(invisible(NULL))
    }
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
    images_issues(combine_issues(images_import_issues(), res))
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
    uscr_nps_debug(make_model_debug("USCR", "uploaded data"))
    rem_nps_debug(make_model_debug("REM", "uploaded data"))
    tte_nps_debug(make_model_debug("TTE", "uploaded data"))
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
    imported <- run_qc_step(
      read_csv(
        input$deployment_csv$datapath,
        show_col_types = FALSE,
        name_repair = "unique_quiet"
      ) |>
        clean_deployment_import()
    )
    deployment_raw(imported$value)
    deployment_import_issues(imported[c("messages", "warnings")])
    run_deployment_check()
    clear_uploaded_fits("uploaded a new deployment file")

    if (!is.null(images_issues())) {
      images_raw(NULL)
      images_import_issues(NULL)
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
      session$sendCustomMessage("reset_file_input", "images_csv")
      return(invisible(NULL))
    }
    imported <- run_qc_step(
      read_csv(
        input$images_csv$datapath,
        show_col_types = FALSE,
        name_repair = "unique_quiet",
        col_types = cols(
          `Sighting Count` = col_character(),
          Species = col_character()
        )
      ) |>
        clean_images_import()
    )
    images_raw(imported$value)
    images_import_issues(imported[c("messages", "warnings")])
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
  observeEvent(input$species_to_analyze, {
    clear_uploaded_fits("changed the species to analyze")
  }, ignoreInit = TRUE)

  # Results are cleared only when the trim the models would use changes: not
  # while the entry is invalid, and not when the number of days changes while
  # the trim is off. Debounced so typing a number of days clears results once.
  trim_days_used <- debounce(reactive(list(app_trim_days())), 1000)
  last_trim_days <- NULL
  observeEvent(trim_days_used(), {
    trim <- trim_days_used()
    if (!is.null(last_trim_days) && !identical(trim, last_trim_days)) {
      clear_uploaded_fits("changed the deployment trim")
    }
    last_trim_days <<- trim
  })

  output$deployment_check_log <- renderText({
    req(deployment_issues())
    format_check_log(deployment_issues())
  })

  output$deployment_preview <- renderDT({
    req(deployment_checked())
    datatable(deployment_checked(), options = list(pageLength = 10))
  })

  # The analysis window each camera contributes (see analysis_windows()).
  output$deployment_windows <- renderDT({
    req(deployment_checked())
    win <- analysis_windows(deployment_checked(), max_days = app_trim_days())
    # What set each window's start and end; a trim overrides the end.
    flags <- cbind(
      "Start Time missing" = win$start_time_missing,
      "End Time missing" = win$end_time_missing & !win$trimmed,
      "Camera malfunction" = win$malfunction & !win$trimmed,
      "Trimmed" = win$trimmed
    )
    notes <- apply(flags, 1, \(f) paste(colnames(flags)[f], collapse = "; "))
    tbl <- data.frame(
      `Site Name` = win$`Site Name`,
      `Window start` = format(win$window_start, "%m/%d/%Y %H:%M"),
      `Window end` = format(win$window_end, "%m/%d/%Y %H:%M"),
      Days = round(win$days, 2),
      Notes = notes,
      check.names = FALSE
    )
    datatable(tbl, options = list(pageLength = 10), rownames = FALSE)
  })

  output$download_deployment_checked <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$deployment_csv$name), "_CHECKED.csv")
    },
    content = function(file) {
      req(deployment_checked())
      write_csv(deployment_checked(), file)
    }
  )

  output$images_check_log <- renderText({
    req(images_issues())
    format_check_log(images_issues())
  })

  output$images_preview <- renderDT({
    req(images_checked())
    datatable(images_checked(), options = list(pageLength = 10))
  })

  output$download_images_checked <- downloadHandler(
    filename = function() {
      paste0(tools::file_path_sans_ext(input$images_csv$name), "_CHECKED.csv")
    },
    content = function(file) {
      req(images_checked())
      write_csv(images_checked(), file)
    }
  )

  # ================================================================
  # DATA SUMMARY (uploaded data)
  # ================================================================

  deploy_summary <- reactive({
    req(deployment_checked())
    summarize_deployments(deployment_checked(), max_days = app_trim_days())
  })

  # Detection events in the uploaded images, counted the way the models count
  # them (see detection_events()).
  uploaded_events <- reactive({
    req(deployment_checked(), images_checked())
    detection_events(images_checked(), deployment_checked(), max_days = app_trim_days())
  })

  image_summary <- reactive({
    req(images_checked())
    summarize_species(images_checked(), uploaded_events())
  })

  # Species in the uploaded images; empty when no images are loaded (for
  # example, after a new deployment file clears them).
  available_species <- reactive({
    images <- images_checked()
    if (is.null(images)) return(character(0))
    species <- sort(unique(images$Species))
    species[!is.na(species) & nzchar(species)]
  })

  # Refill both species dropdowns, keeping each one's current choice when
  # that species is still in the list.
  observe({
    species <- available_species()
    keep_or_first <- function(current) {
      if (length(current) == 1 && current %in% species) current else species[1]
    }
    updateSelectInput(
      session,
      "summary_species",
      choices = species,
      selected = keep_or_first(isolate(input$summary_species))
    )
    updateSelectInput(
      session,
      "species_to_analyze",
      choices = species,
      selected = keep_or_first(isolate(input$species_to_analyze))
    )
  })

  # req() the uploads too: after a new deployment file clears the images (or
  # fails its checks), the species dropdown can still hold its old choice.
  species_objects <- reactive({
    req(deployment_checked(), images_checked())
    species_name <- req(input$summary_species)
    list(
      by_camera = species_by_camera(
        images_checked(), uploaded_events(), deployment_checked(), species_name
      ),
      by_day = species_by_day(uploaded_events(), species_name)
    )
  })

  output$deploy_summary_table <- renderDT({
    req(deploy_summary())
    summary_df <- deploy_summary() |>
      mutate(across(c(mean_days, total_camera_days), \(x) round(x, 1))) |>
      rename(
        `Number of cameras` = n_cameras,
        `Mean days per camera` = mean_days,
        `Total camera-days` = total_camera_days
      )
    datatable(
      summary_df,
      options = list(dom = "t", paging = FALSE, autoWidth = TRUE),
      rownames = FALSE
    )
  })

  output$camera_map <- renderLeaflet({
    req(deployment_checked())
    camera_icon <- makeAwesomeIcon(
      icon = "camera", markerColor = "green",
      library = "fa", iconColor = "black"
    )

    leaflet(deployment_checked()) |>
      addTiles() |>
      addAwesomeMarkers(
        lng   = ~Longitude,
        lat   = ~Latitude,
        icon  = camera_icon,
        popup = ~as.character(`Site Name`),
        label = ~as.character(`Site Name`)
      )
  })

  output$image_summary_table <- renderDT({
    req(image_summary())
    summary_df <- image_summary() |>
      rename(
        Images = n_images,
        `Detection events` = n_events,
        `Animals counted` = animals,
        `Cameras with events` = cameras
      )
    datatable(
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

    ggplot(df, aes(x = reorder(Species, -animals), y = animals)) +
      geom_col(fill = redwood_colors[3]) +
      labs(x = "Species", y = "Animals counted") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$species_site_table <- renderDT({
    species_df <- species_objects()$by_camera |>
      select(`Site Name`, n_images, n_events, animals) |>
      rename(
        Images = n_images,
        `Detection events` = n_events,
        `Animals counted` = animals
      )
    datatable(
      species_df,
      options = list(pageLength = 10, autoWidth = TRUE, dom = "tp"),
      rownames = FALSE
    )
  })

  output$species_site_map <- renderPlot({
    species_name <- req(input$summary_species)
    dc <- species_objects()$by_camera

    ggplot(dc, aes(x = Longitude, y = Latitude)) +
      geom_point(aes(size = animals),
                 color = redwood_colors[3], alpha = 0.7) +
      geom_text_repel(
        aes(label = `Site Name`),
        size = 3,
        max.overlaps = Inf
      ) +
      scale_size_continuous(
        name   = "Animals counted",
        range  = c(2, 14),
        limits = c(0, NA),
        breaks = integer_breaks
      ) +
      coord_quickmap() +
      labs(
        title = paste(species_name, "counted by site"),
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
  })

  output$species_daily_plot <- renderPlot({
    species_name <- req(input$summary_species)
    dd <- species_objects()$by_day
    validate(
      need(
        nrow(dd) > 0,
        paste0("No '", species_name, "' detection events inside the analysis windows.")
      )
    )

    ggplot(dd, aes(x = Date, y = `Site Name`, size = animals)) +
      geom_point(color = redwood_colors[3], alpha = 0.7) +
      scale_size_continuous(name = "Animals counted", range = c(2, 8), breaks = integer_breaks) +
      labs(
        title = paste(species_name, "counted per day by site"),
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
  # MODEL SETTINGS
  # ================================================================

  nps_model_inputs <- reactive({
    req(deployment_checked(), images_checked(), input$species_to_analyze)
    build_nps_model_inputs(
      deployment_checked(),
      images_checked(),
      input$species_to_analyze,
      max_days = app_trim_days()
    )
  })

  output$species_count_summary <- renderTable({
    d <- nps_model_inputs()
    tibble(
      category = c("Animals counted", "Total camera-days"),
      result = c(
        as.character(sum(d$camera_counts)),
        format_fixed(sum(d$camera_days), 1)
      )
    )
  }, colnames = FALSE)

  # ---- Effective settings ----
  # Default mode always uses the default priors; Advanced mode uses the inputs.
  is_default_mode <- function() !identical(input$mode, "Advanced")
  prior_input <- function(id) {
    if (is_default_mode()) advanced_prior_defaults[[id]] else input[[id]]
  }
  app_uscr_buffer_m <- function() {
    if (is_default_mode()) default_uscr_buffer_m else input$uscr_buffer_m
  }

  # Switching to Default resets every Advanced input to its default.
  observeEvent(input$mode, {
    if (!is_default_mode()) return()
    for (id in names(advanced_prior_defaults)) {
      updateNumericInput(session, id, value = advanced_prior_defaults[[id]])
    }
    updateNumericInput(session, "uscr_buffer_m", value = default_uscr_buffer_m)
    updateCheckboxInput(session, "D_prior_uniform", value = TRUE)
    updateCheckboxInput(session, "camera_noise", value = FALSE)
    for (id in c("D_prior_lo", "D_prior_hi")) {
      updateNumericInput(session, id, value = NA)
    }
  }, ignoreInit = TRUE)

  # Density prior (all models): uniform in Default mode or when the box is
  # checked; otherwise the user's 95% interval (animals/km^2).
  density_prior_spec <- function() {
    if (is_default_mode() || !identical(input$D_prior_uniform, FALSE)) {
      return(list(uniform = TRUE))
    }
    list(
      uniform = FALSE,
      lo = input$D_prior_lo,
      hi = input$D_prior_hi
    )
  }

  # Camera-level noise (Poisson-lognormal): off in Default mode and by default.
  app_camera_noise <- function() {
    !is_default_mode() && isTRUE(input$camera_noise)
  }

  app_max_adapt_rounds <- function() {
    x <- suppressWarnings(as.integer(input$max_adapt_rounds))
    if (length(x) != 1L || is.na(x) || x < 1L) mcmc_defaults$max_adapt_rounds else x
  }

  # Problem with the priors a model uses, or NULL. Means must be numbers;
  # SDs and the gamma shape and rate must be positive.
  prior_labels <- c(
    log_v_mean = "log(v) mean",
    log_v_sd = "log(v) SD",
    log_sigma_mean = "log(σ) mean",
    log_sigma_sd = "log(σ) SD",
    log_lam0_mean = "log(λ₀) mean",
    log_lam0_sd = "log(λ₀) SD",
    sd_eps_shape = "σε gamma shape",
    sd_eps_rate = "σε gamma rate"
  )
  prior_problem <- function(model) {
    ids <- c(
      if (identical(model, "USCR")) {
        c("log_sigma_mean", "log_sigma_sd", "log_lam0_mean", "log_lam0_sd")
      } else {
        c("log_v_mean", "log_v_sd")
      },
      "sd_eps_shape", "sd_eps_rate"
    )
    valid <- vapply(ids, function(id) {
      x <- prior_input(id)
      is.numeric(x) && length(x) == 1 && is.finite(x) && (grepl("_mean$", id) || x > 0)
    }, logical(1))
    if (all(valid)) return(NULL)
    paste0(
      "Check the priors in Model settings > Advanced: ",
      paste(prior_labels[ids[!valid]], collapse = ", "),
      if (sum(!valid) == 1) " isn't" else " aren't",
      " valid. Means must be numbers; SDs and the gamma shape and rate must be positive numbers."
    )
  }

  # Problem with the Model settings a model uses, or NULL.
  model_settings_problem <- function(model) {
    if (identical(model, "USCR") && !is_valid_buffer_m(app_uscr_buffer_m())) {
      return(uscr_buffer_problem)
    }
    problem <- prior_problem(model)
    if (!is.null(problem)) return(problem)
    problem <- density_prior_problem(density_prior_spec())
    if (!is.null(problem)) paste("Density prior:", problem)
  }

  # Input summary for USCR runs. buffer_m, M, and camera_noise are the
  # settings the run starts with.
  summarize_uscr_context <- function(d, buffer_m, M, camera_noise) {
    run_context_lines(
      d,
      "USCR",
      model_lines = c(
        paste("USCR state-space buffer:", format_fixed(buffer_m, 0), "m"),
        paste("USCR density prior:", describe_density_prior_spec(density_prior_spec())),
        paste("Camera-level noise:", noise_label(camera_noise))
      ),
      limit_setting = paste0("M=", M),
      max_adapt_rounds = app_max_adapt_rounds()
    )
  }

  # Input summary for REM and TTE runs. angle_summary describes the cameras'
  # detection angles, and camera_noise is the run's noise setting.
  summarize_encounter_context <- function(d, model, angle_summary, camera_noise) {
    density_prior <- density_prior_spec()
    distance_m <- d$out$`Detection Distance`
    run_context_lines(
      d,
      model,
      model_lines = c(
        paste0(
          "Detection distance range: ", format_fixed(min(distance_m), 1),
          " to ", format_fixed(max(distance_m), 1), " m"
        ),
        angle_summary,
        paste(model, "density prior:", describe_density_prior_spec(density_prior)),
        paste("Camera-level noise:", noise_label(camera_noise))
      ),
      limit_setting = paste0("D_max=", format_fixed(encounter_start_D_max(density_prior), 1)),
      max_adapt_rounds = app_max_adapt_rounds()
    )
  }

  # Model calls shared by simulated and uploaded data, each returned as a
  # function of the status callback for run_model(). d holds the cameras
  # (out), camera_counts, and camera_days.
  uscr_fitter <- function(d, M, camera_noise, compute_WAIC) {
    function(status_callback) {
      run_USCR(
        out = d$out,
        camera_counts = d$camera_counts,
        camera_days = d$camera_days,
        iter = mcmc_defaults$iter,
        burnin = mcmc_defaults$burnin,
        thin = mcmc_defaults$thin,
        n_chains = mcmc_defaults$n_chains,
        M = M,
        log_sigma_mean = prior_input("log_sigma_mean"),
        log_sigma_sd = prior_input("log_sigma_sd"),
        log_lam0_mean = prior_input("log_lam0_mean"),
        log_lam0_sd = prior_input("log_lam0_sd"),
        sd_eps_shape = prior_input("sd_eps_shape"),
        sd_eps_rate = prior_input("sd_eps_rate"),
        buffer_m = app_uscr_buffer_m(),
        density_prior = density_prior_spec(),
        camera_noise = camera_noise,
        max_adapt_rounds = app_max_adapt_rounds(),
        adapt_checks = adapt_checks,
        compute_WAIC = compute_WAIC,
        status_callback = status_callback
      )
    }
  }

  # theta_deg is one detection angle for every camera, or one per camera.
  encounter_fitter <- function(model, d, theta_deg, camera_noise, compute_WAIC) {
    run_fn <- if (identical(model, "REM")) run_REM else run_TTE
    function(status_callback) {
      run_fn(
        y = d$camera_counts,
        r_km = d$out$`Detection Distance` / 1000,
        camera_days = d$camera_days,
        theta_deg = theta_deg,
        iter = mcmc_defaults$iter,
        burnin = mcmc_defaults$burnin,
        thin = mcmc_defaults$thin,
        n_chains = mcmc_defaults$n_chains,
        D_max = encounter_start_D_max(density_prior_spec()),
        log_v_mean = prior_input("log_v_mean"),
        log_v_sd = prior_input("log_v_sd"),
        sd_eps_shape = prior_input("sd_eps_shape"),
        sd_eps_rate = prior_input("sd_eps_rate"),
        density_prior = density_prior_spec(),
        camera_noise = camera_noise,
        max_adapt_rounds = app_max_adapt_rounds(),
        adapt_checks = adapt_checks,
        compute_WAIC = compute_WAIC,
        status_callback = status_callback
      )
    }
  }

  # ---- Density prior preview (Model settings > Advanced > Priors: Density) ----
  # USCR's prior on density depends on the state-space area (and M), so its
  # preview uses the cameras in the uploaded data if available, otherwise
  # those in the current simulation.
  uscr_prior_preview_basis <- reactive({
    d <- tryCatch(nps_model_inputs(), error = \(e) NULL)
    if (!is.null(d) && nrow(d$out) >= 2) {
      return(list(label = "your uploaded data", out = d$out))
    }
    s <- pshared_sim()
    if (!is.null(s)) {
      return(list(label = "the current simulation", out = s$uscr$out))
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
    if (!is.null(problem)) return(list(message = problem))

    beta <- tryCatch(
      design_beta_density_prior(spec$lo, spec$hi, encounter_start_D_max(spec)),
      error = \(e) e
    )
    if (inherits(beta, "error")) {
      return(list(message = conditionMessage(beta)))
    }
    out <- list(beta = beta)

    # USCR's prior depends on the state-space area, so it needs camera data.
    basis <- uscr_prior_preview_basis()
    if (is.null(basis)) {
      out$uscr_message <- "its prior depends on the size of the state space, so it appears here once you upload data or simulate data."
      return(out)
    }
    buffer_m <- app_uscr_buffer_m()
    if (!is_valid_buffer_m(buffer_m)) {
      out$uscr_message <- uscr_buffer_problem
      return(out)
    }
    area <- uscr_state_space_and_area(basis$out, buffer_m = buffer_m)$area_km2
    M <- uscr_start_M(basis$out, buffer_m, spec)
    uscr <- tryCatch(design_uscr_density_prior(spec$lo, spec$hi, area, M), error = \(e) e)
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
      return(tags$p(class = "small", style = "color: firebrick;", pv$message))
    }
    fmt <- function(x) format(signif(x, 3), big.mark = ",")
    prior_line <- function(label, d, note) {
      tags$p(
        class = "small",
        tags$strong(label),
        paste0(
          "95% interval ", fmt(d$achieved[["lo"]]), " to ", fmt(d$achieved[["hi"]]),
          " animals/km², mean ", fmt(d$achieved[["mean"]]),
          ", median ", fmt(d$achieved[["median"]]), "."
        ),
        note
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
      prior_line(
        "REM & TTE: ", pv$beta,
        tagList(
          "(Starting D", tags$sub("max", .noWS = c("before", "after")),
          " = ", fmt(pv$beta$D_max), " animals/km².)"
        )
      ),
      off_target_note("REM & TTE", pv$beta, ""),
      if (!is.null(pv[["uscr"]])) {
        prior_line(
          "USCR: ", pv[["uscr"]],
          paste0("(Based on ", pv$basis$label, ": state space ", fmt(pv[["uscr"]]$area_km2),
                 " km², starting M = ", pv[["uscr"]]$M, ".)")
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
        "During fitting, the app may enlarge D", tags$sub("max", .noWS = c("before", "after")), ", M, or the USCR state space. When it does, it re-solves the prior to keep your interval, so the prior used in the fit can differ slightly from this preview. The run-status panel reports the prior actually used."
      )
    )
  })

  output$D_prior_plot <- renderPlot({
    pv <- density_prior_preview()
    req(pv, pv$beta)
    plot_density_priors(pv$beta, pv[["uscr"]])
  }, res = 110)

  # ---- "Run with" notes next to each Run button ----
  # Lists what differs from the app defaults for the model that button runs.
  # On simulated data, camera-level noise comes from the simulation.
  run_settings_diffs <- function(model, shared = FALSE) {
    fmt <- function(x) format(signif(as.numeric(x), 4))
    changed <- function(x, default, tol = 1e-8) {
      is.numeric(x) && length(x) == 1 && is.finite(x) && abs(x - default) > tol
    }
    diffs <- character(0)
    add_prior <- function(id) {
      val <- prior_input(id)
      if (changed(val, advanced_prior_defaults[[id]])) {
        diffs <<- c(diffs, paste0(
          prior_labels[[id]], ": ", fmt(val),
          " (default ", fmt(advanced_prior_defaults[[id]]), ")"
        ))
      }
    }

    spec <- density_prior_spec()
    if (!isTRUE(spec$uniform)) {
      diffs <- c(diffs, paste0("Density prior: 95% interval ", fmt(spec$lo), "–", fmt(spec$hi), " animals/km²"))
    }

    noise_on <- if (shared) isTRUE(pshared_sim()$settings$camera_noise) else app_camera_noise()
    if (noise_on) {
      diffs <- c(diffs, paste0(
        "Camera-level noise: on, σ<sub>ε</sub> ~ Gamma(", fmt(prior_input("sd_eps_shape")), ", ",
        fmt(prior_input("sd_eps_rate")), ")", if (shared) " (from the simulated data)" else ""
      ))
    }

    if (identical(model, "USCR")) {
      add_prior("log_sigma_mean")
      add_prior("log_sigma_sd")
      add_prior("log_lam0_mean")
      add_prior("log_lam0_sd")
      # Buffer: flag only a user override of the automatic default.
      auto <- round_buffer_m(uscr_default_buffer_m(prior_input("log_sigma_mean"), prior_input("log_sigma_sd")))
      if (changed(app_uscr_buffer_m(), auto)) {
        diffs <- c(diffs, paste0("State-space buffer: ", fmt(app_uscr_buffer_m()), " m (default ", fmt(auto), " m)"))
      }
    } else {
      add_prior("log_v_mean")
      add_prior("log_v_sd")
    }

    if (app_max_adapt_rounds() != mcmc_defaults$max_adapt_rounds) {
      diffs <- c(diffs, paste0(
        "Maximum adaptive rounds: ", app_max_adapt_rounds(),
        " (default ", mcmc_defaults$max_adapt_rounds, ")"
      ))
    }
    diffs
  }

  # Prior means (from Model settings) that differ from the logs of the values
  # the current simulated data were generated with. The slider defaults are
  # rounded to the slider steps, so at default settings the true values
  # differ slightly from the default prior means (e.g., log(3.1) = 1.131 vs
  # the default log(v) mean of 1.130); a coarser tolerance keeps those off
  # the list.
  sim_prior_mean_diffs <- function(model) {
    fmt <- function(x) format(signif(as.numeric(x), 4))
    st <- pshared_sim()$settings
    diffs <- character(0)
    add_mean <- function(id, label, true_value, true_label) {
      prior_mean <- prior_input(id)
      differs <- is.numeric(prior_mean) && length(prior_mean) == 1 &&
        is.finite(prior_mean) && abs(prior_mean - log(true_value)) > 0.01
      if (differs) {
        diffs <<- c(diffs, paste0(
          label, " mean: ", fmt(prior_mean), " (data-generating value ", fmt(log(true_value)),
          ", from ", true_label, ")"
        ))
      }
    }
    if (identical(model, "USCR")) {
      add_mean("log_sigma_mean", "log(σ)", st$sigma_km, paste0("σ = ", fmt(st$sigma_km), " km"))
      add_mean("log_lam0_mean", "log(λ₀)", st$lambda0, paste0("λ₀ = ", fmt(st$lambda0)))
    } else {
      add_mean("log_v_mean", "log(v)", st$v, paste0("v = ", fmt(st$v), " km/day"))
    }
    diffs
  }

  render_run_note <- function(model, shared = FALSE) {
    renderUI({
      if (shared) req(pshared_sim())
      diffs <- run_settings_diffs(model, shared)
      mean_diffs <- if (shared) sim_prior_mean_diffs(model)
      buffer_warning <- if (shared && identical(model, "USCR")) {
        sim_buffer_warning(app_uscr_buffer_m(), pshared_sim()$settings$buffer_m)
      }
      # HTML() so the <sub> in sigma_epsilon renders; the lists hold only app
      # text and numbers.
      note_list <- function(items) {
        tags$ul(style = "margin: 0.2rem 0 0 0;", lapply(items, \(d) tags$li(HTML(d))))
      }
      tags$div(
        class = "small",
        style = "margin: 0.4rem 0 0.5rem 0; color: var(--muted);",
        if (!length(diffs)) {
          tags$div(tags$strong("Run with all default settings."))
        } else {
          tags$div(tags$strong("The following differ from the defaults:"), note_list(diffs))
        },
        if (length(mean_diffs)) {
          tags$div(
            style = "margin-top: 0.3rem;",
            tags$strong("Prior means that differ from the data-generating values:"),
            note_list(mean_diffs)
          )
        },
        if (!is.null(buffer_warning)) {
          tags$div(style = "margin-top: 0.3rem; color: firebrick;", buffer_warning)
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

  # Keep the USCR state-space buffer at its automatic value (the rule in
  # uscr_buffer_rule_text) whenever the log(sigma) prior mean or SD changes.
  # A buffer the user types in stays until the log(sigma) prior changes
  # again.
  observe({
    req(
      is.finite(input$log_sigma_mean),
      is.finite(input$log_sigma_sd),
      input$log_sigma_sd > 0
    )
    new_default <- uscr_default_buffer_m(input$log_sigma_mean, input$log_sigma_sd)
    updateNumericInput(session, "uscr_buffer_m", value = round_buffer_m(new_default))
  })

  # ================================================================
  # MODEL FITS ON UPLOADED DATA (USCR, REM, and TTE tabs)
  # ================================================================

  uscr_nps_fit <- reactiveVal(NULL)
  rem_nps_fit <- reactiveVal(NULL)
  tte_nps_fit <- reactiveVal(NULL)
  uscr_nps_debug <- reactiveVal(make_model_debug("USCR", "uploaded data"))
  rem_nps_debug <- reactiveVal(make_model_debug("REM", "uploaded data"))
  tte_nps_debug <- reactiveVal(make_model_debug("TTE", "uploaded data"))

  # Checks shared by the uploaded-data runs. Returns the model inputs, or NULL
  # after a pop-up explaining why the model can't be run.
  nps_run_inputs <- function(model) {
    problem <- if (is.null(deployment_raw())) {
      "Upload a deployment file on the Add your data tab first."
    } else if (is.null(deployment_checked())) {
      "The deployment file didn't pass its checks. See the Deployment check log on the Add your data tab."
    } else if (is.null(images_raw())) {
      "Upload an images file on the Add your data tab first."
    } else if (is.null(images_checked())) {
      "The images file didn't pass its checks. See the Images check log on the Add your data tab."
    } else if (!is.null(trim_days_problem())) {
      paste(trim_days_problem(), "See the deployment trim on the Add your data tab.")
    } else if (!isTruthy(input$species_to_analyze)) {
      "Choose a species to analyze in Model settings."
    }
    if (!is.null(problem)) {
      showNotification(paste(model, "not run.", problem), type = "error", duration = NULL)
      return(NULL)
    }

    d <- nps_model_inputs()
    problem <- if (nrow(d$out) < 2) {
      "The models need data from at least 2 cameras."
    } else if (!all(d$camera_days > 0)) {
      "Some cameras have zero or negative camera-days; check deployment dates, or check for overlapping/duplicate Cluster IDs at the same site (overlapping detection events can double-count time and drive camera-days to zero or below)."
    } else {
      model_settings_problem(model)
    }
    if (!is.null(problem)) {
      showNotification(paste(model, "not run.", problem), type = "error", duration = NULL)
      return(NULL)
    }
    d
  }

  # Run a model on the uploaded data, with WAIC for the Compare & combine tab.
  run_on_uploaded <- function(model) {
    d <- nps_run_inputs(model)
    if (is.null(d)) return()
    species <- input$species_to_analyze
    camera_noise <- app_camera_noise()

    if (identical(model, "USCR")) {
      buffer_m <- app_uscr_buffer_m()
      M <- uscr_start_M(d$out, buffer_m, density_prior_spec())
      context <- summarize_uscr_context(d, buffer_m, M, camera_noise)
      fit_model <- uscr_fitter(d, M, camera_noise, compute_WAIC = TRUE)
    } else {
      angle_summary <- summarize_camera_angle_source(d$out, default_deg = input$theta)
      context <- summarize_encounter_context(d, model, angle_summary, camera_noise)
      theta_deg <- get_camera_angle_vector(d$out, default_deg = input$theta)
      fit_model <- encounter_fitter(model, d, theta_deg, camera_noise, compute_WAIC = TRUE)
    }

    # Record the species with the fit, for the results printout and the
    # Compare & combine tab.
    fit_species <- function(status_callback) {
      fit <- fit_model(status_callback)
      fit$species <- species
      fit
    }

    run_model(
      model = model,
      data_label = "uploaded data",
      fit_rv = list(USCR = uscr_nps_fit, REM = rem_nps_fit, TTE = tte_nps_fit)[[model]],
      debug_rv = list(USCR = uscr_nps_debug, REM = rem_nps_debug, TTE = tte_nps_debug)[[model]],
      context = c(paste("Species:", species), context),
      fit_model = fit_species
    )
  }

  observeEvent(input$run_uscr_nps, run_on_uploaded("USCR"))
  observeEvent(input$run_rem_nps, run_on_uploaded("REM"))
  observeEvent(input$run_tte_nps, run_on_uploaded("TTE"))

  output$uscr_nps_result <- renderUI({
    uploaded_result_ui("USCR", uscr_nps_fit(), uscr_nps_debug())
  })

  output$rem_nps_result <- renderUI({
    uploaded_result_ui("REM", rem_nps_fit(), rem_nps_debug())
  })

  output$tte_nps_result <- renderUI({
    uploaded_result_ui("TTE", tte_nps_fit(), tte_nps_debug())
  })

  output$uscr_nps_debug <- renderText({
    format_model_debug(uscr_nps_debug(), uscr_nps_fit())
  })

  output$rem_nps_debug <- renderText({
    format_model_debug(rem_nps_debug(), rem_nps_fit())
  })

  output$tte_nps_debug <- renderText({
    format_model_debug(tte_nps_debug(), tte_nps_fit())
  })

  # ================================================================
  # COMPARE & COMBINE
  # ================================================================
  # Combined once, so the table and plots use the same mixture draws.
  nps_combo_data <- reactive({
    prepare_combo_density_data(list(
      USCR = uscr_nps_fit(),
      REM = rem_nps_fit(),
      TTE = tte_nps_fit()
    ))
  })

  nps_combo_interval_plot_obj <- reactive({
    combo_data <- nps_combo_data()
    if (!is.null(combo_data)) {
      build_combo_interval_plot(
        combo_data,
        title_text = paste0(combo_data$species, " (uploaded data): density estimates")
      )
    }
  })

  nps_combo_density_plot_obj <- reactive({
    combo_data <- nps_combo_data()
    if (!is.null(combo_data)) {
      build_combo_density_overlay_plot(
        combo_data,
        title_text = paste0(combo_data$species, " (uploaded data): posterior distributions of density")
      )
    }
  })

  output$nps_combo_table <- renderDT({
    combo_data <- nps_combo_data()
    if (is.null(combo_data)) {
      return(datatable(
        data.frame(
          Note = "Run at least one model on uploaded data from its tab. The table updates as USCR, REM, and TTE finish."
        ),
        options = list(dom = "t", paging = FALSE),
        rownames = FALSE
      ))
    }
    df <- build_combo_table(
      combo_data,
      density_threshold = input$combo_density_threshold,
      ppc_pass_range = ppc_pass_range
    )
    # Bayesian p-values of models that failed their posterior predictive check
    # in bold red, set by the hidden failed_check column.
    datatable(
      df,
      caption = paste("Species:", combo_data$species),
      options = list(
        paging = FALSE,
        ordering = FALSE,
        dom = "t",
        autoWidth = TRUE,
        columnDefs = list(list(
          visible = FALSE,
          targets = which(names(df) == "failed_check") - 1L
        ))
      ),
      rownames = FALSE
    ) |>
      format_result_columns(df) |>
      formatStyle(
        "Bayesian p-value",
        valueColumns = "failed_check",
        color = styleEqual(TRUE, "firebrick"),
        fontWeight = styleEqual(TRUE, "bold")
      )
  })

  output$nps_combo_interval_plot <- renderPlot({
    p <- nps_combo_interval_plot_obj()
    validate(need(!is.null(p), "Run at least one model on uploaded data first."))
    p
  }, res = 120)

  output$nps_combo_density_plot <- renderPlot({
    p <- nps_combo_density_plot_obj()
    validate(need(!is.null(p), "Run at least one model on uploaded data first."))
    p
  }, res = 120)

  output$dl_nps_combo_interval_png <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_interval_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- nps_combo_interval_plot_obj()
      req(p)
      ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )

  output$dl_nps_combo_density_png <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_density_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- nps_combo_density_plot_obj()
      req(p)
      ggsave(file, plot = p, width = 8.5, height = 4.5, dpi = 300, bg = "white")
    }
  )

  output$dl_nps_all_csv <- downloadHandler(
    filename = function() {
      paste0("DEER_uploaded_data_posterior_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      combo_data <- nps_combo_data()
      req(combo_data)
      imap(combo_data$fits, function(fit, model) {
        mutate(posterior_summary_df(fit), species = fit$species, model = model, .before = 1)
      }) |>
        bind_rows() |>
        write_csv(file)
    }
  )
}

# -------------------------------------------------------------------
# Run the app
# -------------------------------------------------------------------

shinyApp(ui, server)
