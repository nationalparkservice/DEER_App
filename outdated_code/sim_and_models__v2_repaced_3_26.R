## R/sim_and_models.R
## Simulation helpers + NPS model inputs + REM / USCR / TTE models
## Models follow Camera_Unmarked_Analyses_updated.Rmd (REM/TTE/USCR: nimbleMCMC).

# -------------------------------------------------------------------
# Package checks (loaded lazily where needed)
# -------------------------------------------------------------------

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }
}

# -------------------------------------------------------------------
# Drop-in rewritten USCR implementation
# -------------------------------------------------------------------
# This repo includes a separate "rewritten_uscr_dropin.R" that provides
# a speed-aware adaptive USCR rerun logic (Rhat/psi/M) and a full USCR
# wrapper. We source it here and keep its USCR entry point as
# `run_USCR_dropin`, then define an app-compatible `run_USCR()` below.
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
script_dir <- if (is.na(this_file)) getwd() else dirname(this_file)

# rewritten_uscr_dropin.R lives in `R/` relative to app.R.
# Depending on how R is invoked, script_dir may be repo root or `R/`,
# so we check multiple candidates.
candidate_paths <- c(
  file.path(script_dir, "rewritten_uscr_dropin.R"),
  file.path(script_dir, "R", "rewritten_uscr_dropin.R"),
  file.path(getwd(), "R", "rewritten_uscr_dropin.R")
)
dropin_path <- candidate_paths[file.exists(candidate_paths)][1]
if (is.na(dropin_path) || is.null(dropin_path)) {
  stop(
    "Could not locate rewritten_uscr_dropin.R. Looked in: ",
    paste(candidate_paths, collapse = ", "),
    call. = FALSE
  )
}
source(dropin_path)
run_USCR_dropin <- run_USCR

# -------------------------------------------------------------------
# 1. Simulation helpers (for teaching/sim tab)
# -------------------------------------------------------------------
simulate_camera_counts <- function(n_side   = 5,
                                   spacing_m = 300,
                                   days      = 21,
                                   D_per_km2 = 25,
                                   lambda0   = 0.20,
                                   sigma_m   = 150,
                                   seed      = 1) {
  quiet_require("secr")
  
  set.seed(seed)
  
  traps_obj <- secr::make.grid(
    nx       = n_side,
    ny       = n_side,
    spacing  = spacing_m,
    detector = "count"
  )
  
  mask_obj <- secr::make.mask(
    traps   = traps_obj,
    buffer  = 4 * sigma_m,
    spacing = spacing_m / 3
  )
  
  # secr uses animals/ha when coords are in metres
  D_ha <- D_per_km2 / 100
  
  ch <- secr::sim.capthist(
    traps      = traps_obj,
    popn       = list(D = D_ha),
    detectfn   = "HHN",
    detectpar  = list(lambda0 = lambda0, sigma = sigma_m),
    noccasions = days,
    renumber   = FALSE
  )
  
  list(
    ch    = ch,
    traps = traps_obj,
    mask  = mask_obj,
    truth = list(
      D_per_km2 = D_per_km2,
      lambda0   = lambda0,
      sigma_m   = sigma_m,
      days      = days,
      spacing_m = spacing_m,
      n_cams    = nrow(traps_obj)
    )
  )
}

get_counts_matrix <- function(ch) {
  ch_dims <- dim(ch)
  if (length(ch_dims) == 3) {
    # [animals, occasions, traps] -> [traps, occasions]
    apply(ch, c(3, 2), sum, na.rm = TRUE)
  } else {
    as.matrix(ch)
  }
}

capthist_to_events <- function(ch) {
  quiet_require("dplyr")
  
  counts_mat <- get_counts_matrix(ch)
  n_occ      <- ncol(counts_mat)
  n_cams     <- nrow(counts_mat)
  
  events_list <- vector("list", length = sum(counts_mat > 0))
  idx <- 1L
  
  for (day in seq_len(n_occ)) {
    for (cam in seq_len(n_cams)) {
      count <- counts_mat[cam, day]
      if (count > 0) {
        for (i in seq_len(count)) {
          events_list[[idx]] <- data.frame(
            camera_id  = paste0("C", cam),
            day        = day,
            group_size = 1L
          )
          idx <- idx + 1L
        }
      }
    }
  }
  
  if (idx == 1L) {
    df <- data.frame(
      camera_id  = character(),
      day        = integer(),
      group_size = integer(),
      date_time  = as.POSIXct(character())
    )
  } else {
    df <- dplyr::bind_rows(events_list)
    origin <- as.POSIXct("2025-01-01 00:00:00", tz = "UTC")
    df <- df |>
      dplyr::mutate(
        date_time = origin + (day - 1) * 24 * 3600 + 12 * 3600
      )
  }
  
  df
}

# For simulations: create "NPS-like" inputs for models
sim_model_inputs <- function(sim,
                             detection_radius_m,
                             days_override = NULL) {
  quiet_require("dplyr")
  quiet_require("secr")
  
  ch       <- sim$ch
  traps_df <- as.data.frame(secr::traps(ch))
  counts   <- get_counts_matrix(ch)
  
  y <- rowSums(counts, na.rm = TRUE)
  days_used <- if (is.null(days_override)) ncol(counts) else days_override
  camera_days <- rep(days_used, length(y))
  
  out <- traps_df |>
    dplyr::mutate(
      Site              = paste0("C", dplyr::row_number()),
      utm_e             = x / 1000,   # treat metres as km for consistency
      utm_n             = y / 1000,
      `Detection Distance` = detection_radius_m,
      `Start Index`     = 1L,
      `End Index`       = days_used
    )
  
  list(
    out           = out,
    camera_counts = y,
    camera_days   = camera_days
  )
}

# -------------------------------------------------------------------
# 2. NPS inputs: format deployments + images as in the Rmd
# -------------------------------------------------------------------
build_nps_model_inputs <- function(deployments, images) {
  quiet_require("dplyr")
  quiet_require("tidyr")
  quiet_require("sf")
  quiet_require("lubridate")
  quiet_require("tibble")
  
  # This assumes format_deployments() is defined in data_checks.R
  if (!exists("format_deployments")) {
    stop("format_deployments() not found. Source data_checks.R first.",
         call. = FALSE)
  }
  
  deps <- format_deployments(deployments)
  
  deps <- deps |>
    dplyr::rename(Site = `Site Name`)
  
  # UTM conversion and centering (same logic as Rmd)
  mean_lon  <- mean(deps$Longitude, na.rm = TRUE)
  utm_zone  <- floor((mean_lon + 180) / 6) + 1
  epsg_code <- 26900 + utm_zone
  
  utm_coords <- deps |>
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4269) |>
    sf::st_transform(crs = epsg_code) |>
    sf::st_coordinates()
  
  deps <- deps |>
    dplyr::mutate(
      utm_e = (utm_coords[, "X"] -
                 mean(range(utm_coords[, "X"]))) / 1000,
      utm_n = (utm_coords[, "Y"] -
                 mean(range(utm_coords[, "Y"]))) / 1000
    )
  
  # Deer-only sequences
  seqs <- images |>
    dplyr::filter(Species == "Deer") |>
    dplyr::group_by(`Cluster ID`) |>
    dplyr::summarise(
      Site           = dplyr::first(`Site Name`),
      latitude       = dplyr::first(Latitude),
      longitude      = dplyr::first(Longitude),
      detection_date = lubridate::date(min(Timestamp)),
      deer_count     = max(as.numeric(`Sighting Count`), na.rm = TRUE),
      start          = min(Timestamp),
      end            = max(Timestamp),
      .groups        = "drop"
    ) |>
    dplyr::mutate(
      cluster_length_days = as.numeric(difftime(end, start, units = "days"))
    )
  
  # Ensure Start/End are Dates (try multiple formats like the Rmd)
  start_dates <- deps$`Start Date`
  end_dates   <- deps$`End Date`
  
  # Helper function to parse dates robustly (handles 2-digit years like "2/3/25")
  parse_date_robust <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXct")) return(as.Date(x))
    
    # Try lubridate first (handles 2-digit years automatically)
    parsed <- lubridate::parse_date_time(
      x,
      orders = c("mdy", "ymd", "m-d-y", "y-m-d", "mdY", "Ymd"),
      quiet = TRUE
    )
    
    # If still NA, try as.Date with common formats (including 2-digit year)
    if (any(is.na(parsed))) {
      parsed <- as.Date(x, tryFormats = c("%m/%d/%y", "%m/%d/%Y", "%Y-%m-%d", "%m-%d-%Y", "%d/%m/%Y"))
    }
    
    # Convert POSIXct to Date
    if (inherits(parsed, "POSIXct")) parsed <- as.Date(parsed)
    
    return(parsed)
  }
  
  start_dates <- parse_date_robust(start_dates)
  end_dates   <- parse_date_robust(end_dates)
  
  # Check for parsing failures
  if (all(is.na(start_dates)) || all(is.na(end_dates))) {
    # Show sample of what we're trying to parse
    cat("Sample Start Date values:", head(deps$`Start Date`, 3), "\n")
    cat("Sample End Date values:", head(deps$`End Date`, 3), "\n")
    stop("Failed to parse Start Date or End Date. Check date format in deployment file.",
         call. = FALSE)
  }
  
  # Warn about partial parsing failures
  if (any(is.na(start_dates)) || any(is.na(end_dates))) {
    n_missing_start <- sum(is.na(start_dates))
    n_missing_end <- sum(is.na(end_dates))
    warning("Could not parse ", n_missing_start, " Start Date(s) and ", 
            n_missing_end, " End Date(s). These rows will be excluded.")
  }
  
  # Filter out rows with missing dates
  valid_rows <- !is.na(start_dates) & !is.na(end_dates)
  if (!all(valid_rows)) {
    deps <- deps[valid_rows, ]
    start_dates <- start_dates[valid_rows]
    end_dates <- end_dates[valid_rows]
    warning("Excluded ", sum(!valid_rows), " deployment row(s) with missing dates.")
  }
  
  # Check that we have valid date ranges
  min_start <- min(start_dates, na.rm = TRUE)
  max_end   <- max(end_dates, na.rm = TRUE)
  
  if (!is.finite(min_start) || !is.finite(max_end)) {
    stop("No valid date range found. All dates failed to parse.", call. = FALSE)
  }
  
  # Zero-count rows across all sites × all days in deployment window
  zero_counts <- expand.grid(
    Site = unique(deps$Site),
    detection_date = seq(
      from = min_start,
      to   = max_end,
      by   = "day"
    )
  ) |>
    dplyr::mutate(
      deer_count         = 0,
      cluster_length_days = 0
    )
  
  counts_time <- seqs |>
    dplyr::bind_rows(zero_counts) |>
    dplyr::group_by(Site, detection_date) |>
    dplyr::summarise(
      detections       = sum(deer_count),
      camera_time_days = sum(cluster_length_days),
      .groups          = "drop"
    )
  
  counts <- counts_time |>
    dplyr::select(-camera_time_days) |>
    tidyr::pivot_wider(
      names_from  = detection_date,
      values_from = detections,
      values_fill = 0,
      names_sort  = TRUE
    )
  
  out <- counts |>
    dplyr::left_join(deps, by = "Site")
  
  detection_matrix <- counts |>
    tibble::column_to_rownames("Site") |>
    as.matrix()
  
  camera_time_matrix <- counts_time |>
    dplyr::select(-detections) |>
    tidyr::pivot_wider(
      names_from  = detection_date,
      values_from = camera_time_days,
      values_fill = 0,
      names_sort  = TRUE
    ) |>
    tibble::column_to_rownames("Site") |>
    as.matrix()
  
  # Create Start Index and End Index columns
  # These map each camera's deployment dates to column indices in the detection matrix
  all_dates <- sort(unique(counts_time$detection_date))
  
  # Parse Start/End dates from out (they should already be formatted by format_deployments)
  out_start_dates <- parse_date_robust(out$`Start Date`)
  out_end_dates   <- parse_date_robust(out$`End Date`)
  
  # Find column indices for each camera's deployment period
  out <- out |>
    dplyr::mutate(
      Start_Date_parsed = out_start_dates,
      End_Date_parsed   = out_end_dates
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      `Start Index` = {
        idx <- which(all_dates == Start_Date_parsed)
        if (length(idx) > 0) idx[1] else 1L
      },
      `End Index` = {
        idx <- which(all_dates == End_Date_parsed)
        if (length(idx) > 0) idx[1] else length(all_dates)
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-Start_Date_parsed, -End_Date_parsed)
  
  # Ensure indices are valid
  out$`Start Index` <- pmax(1L, pmin(out$`Start Index`, ncol(detection_matrix)))
  out$`End Index`   <- pmax(out$`Start Index`, pmin(out$`End Index`, ncol(detection_matrix)))
  
  camera_counts <- camera_days <- numeric(nrow(out))
  
  for (i in seq_len(nrow(out))) {
    idx_start <- out$`Start Index`[i]
    idx_end   <- out$`End Index`[i]
    
    if (idx_start > idx_end || idx_start < 1 || idx_end > ncol(detection_matrix)) {
      warning("Invalid indices for camera ", i, ": Start=", idx_start, ", End=", idx_end)
      camera_counts[i] <- 0
      camera_days[i]   <- 0
      next
    }
    
    camera_counts[i] <- sum(detection_matrix[i, idx_start:idx_end])
    camera_days[i]   <- length(idx_start:idx_end) -
      sum(camera_time_matrix[i, idx_start:idx_end])
  }
  
  list(
    out                = out,
    detection_matrix   = detection_matrix,
    camera_time_matrix = camera_time_matrix,
    camera_counts      = camera_counts,
    camera_days        = camera_days
  )
}

# -------------------------------------------------------------------
# 3. REM model (nimble) — mirrors Rmd structure
# -------------------------------------------------------------------
run_REM <- function(y,
                    r_km,
                    camera_days,
                    iter     = 6000,
                    burnin   = 1000,
                    thin     = 5,
                    n_chains = 3,
                    D_max    = 200,
                    log_v_mean = 1.339,
                    log_v_sd   = 0.2955,
                    sd_eps_max = 10,
                    theta_deg  = 55) {
  
  quiet_require("nimble")
  
  J <- length(y)
  if (length(r_km) != J || length(camera_days) != J) {
    stop("run_REM: y, r_km, and camera_days must have same length.", call. = FALSE)
  }
  
  code <- nimble::nimbleCode({
    
    # Priors
    D      ~ dunif(0, D_max)
    log_v  ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dunif(0, sd_eps_max)
    
    v     <- exp(log_v)
    theta <- theta_rad   # passed as constant
    
    for (j in 1:J) {
      
      eps[j] ~ dnorm(0, sd = sd_eps)
      
      log(lambda[j]) <- log(2 + theta) - log(pi_const) + log(v) +
        log(r[j]) + log(D) + log(camera_days[j]) + eps[j]
      
      y[j]     ~ dpois(lambda[j])
      y_sim[j] ~ dpois(lambda[j])
      
      pearson_obs[j] <- (y[j]     - lambda[j])^2 / lambda[j]
      pearson_sim[j] <- (y_sim[j] - lambda[j])^2 / lambda[j]
    }
    
    sum_obs <- sum(pearson_obs[1:J])
    sum_sim <- sum(pearson_sim[1:J])
    bp      <- step(sum_sim - sum_obs)
    
    D_mi2 <- D * 2.59
  })
  
  Const <- list(
    J          = J,
    r          = as.numeric(r_km),
    camera_days = as.numeric(camera_days),
    pi_const   = pi,
    theta_rad  = theta_deg * pi / 180,
    D_max      = D_max,
    log_v_mean = log_v_mean,
    log_v_sd   = log_v_sd,
    sd_eps_max = sd_eps_max
  )
  
  Data <- list(y = as.numeric(y))
  
  monitors <- c("D", "v", "sd_eps",
                "bp", "sum_obs", "sum_sim", "D_mi2")
  
  # Parallel chain execution (like Rmd)
  quiet_require("parallel")
  
  run_MCMC_chain <- function(index, code, constants, data, monitors,
                             niter, nburnin, thin) {
    library(nimble)
    nimble::nimbleMCMC(
      code      = code,
      constants = constants,
      data      = data,
      monitors  = monitors,
      niter     = niter,
      nburnin   = nburnin,
      thin      = thin,
      WAIC      = TRUE
    )
  }
  
  if (n_chains > 1) {
    cl_size <- min(n_chains, parallel::detectCores() - 1)
    cl_size <- max(1, cl_size)
    cl <- parallel::makeCluster(cl_size)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    fits <- parallel::parLapply(
      cl        = cl,
      X         = seq_len(n_chains),
      fun       = run_MCMC_chain,
      code      = code,
      constants = Const,
      data      = Data,
      monitors  = monitors,
      niter     = iter,
      nburnin   = burnin,
      thin      = thin
    )
  } else {
    fits <- list(run_MCMC_chain(1, code, Const, Data, monitors, iter, burnin, thin))
  }
  
  samples_list <- lapply(fits, `[[`, "samples")
  samples_all  <- do.call(rbind, samples_list)
  
  waic_vec <- vapply(
    fits,
    function(f) if (!is.null(f$WAIC)) f$WAIC$WAIC else NA_real_,
    numeric(1)
  )
  waic_mean <- mean(waic_vec, na.rm = TRUE)
  
  list(
    method       = "REM",
    samples_list = samples_list,
    samples_all  = samples_all,
    waic         = waic_mean
  )
}

# -------------------------------------------------------------------
# 4. TTE model (nimble) — mirrors Rmd structure
# -------------------------------------------------------------------
run_TTE <- function(y,
                    r_km,
                    camera_days,
                    iter     = 6000,
                    burnin   = 1000,
                    thin     = 5,
                    n_chains = 3,
                    D_max    = 200,
                    log_v_mean = 1.339,
                    log_v_sd   = 0.2955,
                    sd_eps_max = 10,
                    theta_deg  = 55) {
  
  quiet_require("nimble")
  
  J <- length(y)
  if (length(r_km) != J || length(camera_days) != J) {
    stop("run_TTE: y, r_km, and camera_days must have same length.", call. = FALSE)
  }
  
  code <- nimble::nimbleCode({
    
    D      ~ dunif(0, D_max)
    log_v  ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dunif(0, sd_eps_max)
    
    v     <- exp(log_v)
    theta <- theta_deg  # degrees
    
    for (j in 1:J) {
      
      a[j] <- pi_const * r[j]^2 * theta / 360
      mvw[j] <- 0.59 * r[j]
      time_unit[j]  <- mvw[j] / v
      tte_units[j]  <- camera_days[j] / time_unit[j]
      
      eps[j] ~ dnorm(0, sd = sd_eps)
      
      log(lambda[j]) <- log(D) + log(tte_units[j]) + log(a[j]) + eps[j]
      
      y[j]     ~ dpois(lambda[j])
      y_sim[j] ~ dpois(lambda[j])
      
      pearson_obs[j] <- (y[j]     - lambda[j])^2 / lambda[j]
      pearson_sim[j] <- (y_sim[j] - lambda[j])^2 / lambda[j]
    }
    
    sum_obs <- sum(pearson_obs[1:J])
    sum_sim <- sum(pearson_sim[1:J])
    
    bp <- step(sum_sim - sum_obs)
    
    D_mi2 <- D * 2.59
  })
  
  Const <- list(
    J           = J,
    pi_const    = pi,
    r           = as.numeric(r_km),
    camera_days = as.numeric(camera_days),
    D_max       = D_max,
    log_v_mean  = log_v_mean,
    log_v_sd    = log_v_sd,
    sd_eps_max  = sd_eps_max,
    theta_deg   = theta_deg
  )
  
  Data <- list(y = as.numeric(y))
  
  monitors <- c("D", "v", "sd_eps",
                "bp", "sum_obs", "sum_sim", "D_mi2")
  
  # Parallel chain execution (like Rmd)
  quiet_require("parallel")
  
  run_MCMC_chain <- function(index, code, constants, data, monitors,
                             niter, nburnin, thin) {
    library(nimble)
    nimble::nimbleMCMC(
      code      = code,
      constants = constants,
      data      = data,
      monitors  = monitors,
      niter     = niter,
      nburnin   = nburnin,
      thin      = thin,
      WAIC      = TRUE
    )
  }
  
  if (n_chains > 1) {
    cl_size <- min(n_chains, parallel::detectCores() - 1)
    cl_size <- max(1, cl_size)
    cl <- parallel::makeCluster(cl_size)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    fits <- parallel::parLapply(
      cl        = cl,
      X         = seq_len(n_chains),
      fun       = run_MCMC_chain,
      code      = code,
      constants = Const,
      data      = Data,
      monitors  = monitors,
      niter     = iter,
      nburnin   = burnin,
      thin      = thin
    )
  } else {
    fits <- list(run_MCMC_chain(1, code, Const, Data, monitors, iter, burnin, thin))
  }
  
  samples_list <- lapply(fits, `[[`, "samples")
  samples_all  <- do.call(rbind, samples_list)
  
  waic_vec <- vapply(
    fits,
    function(f) if (!is.null(f$WAIC)) f$WAIC$WAIC else NA_real_,
    numeric(1)
  )
  waic_mean <- mean(waic_vec, na.rm = TRUE)
  
  list(
    method       = "TTE",
    samples_list = samples_list,
    samples_all  = samples_all,
    waic         = waic_mean
  )
}

# -------------------------------------------------------------------
# USCR: state space, buffer, area (Camera_Unmarked_Analyses_updated.Rmd)
# -------------------------------------------------------------------
uscr_state_space_and_area <- function(out,
                                        buffer_chi_p = 0.99,
                                        buffer_scale = 0.40) {
  buffer <- buffer_scale * sqrt(stats::qchisq(buffer_chi_p, df = 2))
  buffer_sq <- buffer ^ 2

  xlim <- c(
    min(out$utm_e, na.rm = TRUE) - buffer,
    max(out$utm_e, na.rm = TRUE) + buffer
  )
  ylim <- c(
    min(out$utm_n, na.rm = TRUE) - buffer,
    max(out$utm_n, na.rm = TRUE) + buffer
  )

  if (all(c("Longitude", "Latitude") %in% names(out))) {
    quiet_require("sf")

    coords_df <- unique(out[, c("Longitude", "Latitude"), drop = FALSE])
    mean_lon <- mean(coords_df$Longitude, na.rm = TRUE)
    utm_zone <- floor((mean_lon + 180) / 6) + 1
    epsg_code <- 26900 + utm_zone

    cam_buff <- sf::st_as_sf(coords_df, coords = c("Longitude", "Latitude"), crs = 4269) |>
      sf::st_transform(crs = epsg_code) |>
      sf::st_buffer(buffer * 1000) |>
      sf::st_union()

    area_mi2 <- as.numeric(sf::st_area(cam_buff)) / (2.59 * 1e6)
  } else {
    # Simulated UTM-only grids: rectangular area (mi^2)
    area_mi2 <- (xlim[2] - xlim[1]) * (ylim[2] - ylim[1]) / 2.59
  }

  list(
    buffer = buffer,
    buffer_sq = buffer_sq,
    xlim = xlim,
    ylim = ylim,
    area_mi2 = area_mi2
  )
}

# -------------------------------------------------------------------
# 5. USCR model (nimble + data augmentation + in_ss) — Rmd parity
# -------------------------------------------------------------------
run_USCR <- function(out,
                     camera_counts,
                     camera_days,
                     iter           = 11000,
                     burnin         = 1000,
                     thin           = 10,
                     n_chains       = 3,
                     M              = 1000,
                     log_sigma_mean = -1.4442,
                     log_sigma_sd   = 0.1451,
                     log_lam0_sd    = 1,
                     sd_eps_shape   = 1,
                     buffer_chi_p   = 0.99,
                     buffer_scale   = 0.40,
                     adaptive        = TRUE,
                     compute_WAIC    = TRUE,
                     diagnostic_mode = FALSE,
                     rhat_target     = 1.1,
                     psi_threshold   = 0.9,
                     psi_prob_cutoff = 0.01,
                     max_adapt_rounds = 2,
                     tuning_n_chains = NULL) {

  # Delegate to the rewritten USCR drop-in while keeping this app's API.
  if (!exists("run_USCR_dropin")) {
    stop("run_USCR_dropin not found; ensure rewritten_uscr_dropin.R was sourced.", call. = FALSE)
  }

  if (is.null(tuning_n_chains)) tuning_n_chains <- min(2L, as.integer(n_chains))
  tuning_n_chains <- as.integer(tuning_n_chains)

  # Speed-aware tuning start: start with a fraction of user-specified iter.
  iter_tune <- max(as.integer(burnin) + 500L, as.integer(floor(as.numeric(iter) / 4)))
  iter_tune <- min(as.integer(iter), iter_tune)

  thin_tune <- 1L

  return(
    run_USCR_dropin(
      out = out,
      camera_counts = camera_counts,
      camera_days = camera_days,
      M = M,
      iter_tune = iter_tune,
      burnin = burnin,
      thin_tune = thin_tune,
      n_chains_tune = tuning_n_chains,
      n_chains_final = as.integer(n_chains),
      final_iter = NULL,
      final_thin = NULL,
      parallel_chains = TRUE,
      compute_WAIC = compute_WAIC,
      adapt = adaptive,
      max_adapt_rounds = max_adapt_rounds,
      rhat_target = rhat_target,
      psi_threshold = psi_threshold,
      psi_prob_cutoff = psi_prob_cutoff,
      monitor_latent = diagnostic_mode,
      log_sigma_mean = log_sigma_mean,
      log_sigma_sd = log_sigma_sd,
      log_lam0_sd = log_lam0_sd,
      sd_eps_shape = sd_eps_shape,
      buffer_chi_p = buffer_chi_p,
      buffer_scale = buffer_scale,
      seed = NULL,
      verbose = FALSE
    )
  )

  if (isTRUE(adaptive)) {
    current_M <- as.integer(M)
    current_iter <- as.integer(iter)
    current_burnin <- as.integer(burnin)
    current_thin <- thin

    if (is.null(tuning_n_chains)) {
      tuning_n_chains <- min(2L, as.integer(n_chains))
    } else {
      tuning_n_chains <- as.integer(tuning_n_chains)
    }

    # Only shorten the tuning/final run when we can actually compute Rhat
    # (i.e., at least 2 chains). Otherwise, keep the user-specified iter
    # to avoid under-sampling.
    if (tuning_n_chains > 1L) {
      current_iter <- min(
        current_iter,
        max(current_burnin + 500L, as.integer(floor(current_iter / 4)))
      )
    }

    adapt_round <- 0L
    rhat_too_high <- TRUE
    M_too_small <- TRUE

    while (adapt_round < as.integer(max_adapt_rounds) &&
           (rhat_too_high || M_too_small)) {
      adapt_round <- adapt_round + 1L

      # Short tuning run: disable WAIC during tuning for speed.
      tune_iter <- current_iter
      tune_burnin <- min(current_burnin, as.integer(floor(tune_iter / 3)))
      tune_burnin <- max(0L, tune_burnin)

      fit_tune <- tryCatch(
        run_USCR(
          out            = out,
          camera_counts  = camera_counts,
          camera_days    = camera_days,
          iter           = tune_iter,
          burnin         = tune_burnin,
          thin           = current_thin,
          n_chains       = tuning_n_chains,
          M              = current_M,
          log_sigma_mean = log_sigma_mean,
          log_sigma_sd   = log_sigma_sd,
          log_lam0_sd    = log_lam0_sd,
          sd_eps_shape   = sd_eps_shape,
          buffer_chi_p   = buffer_chi_p,
          buffer_scale   = buffer_scale,
          adaptive       = FALSE,
          compute_WAIC   = FALSE,
          diagnostic_mode = diagnostic_mode
        ),
        error = function(e) e
      )

      if (inherits(fit_tune, "error") || is.null(fit_tune$samples_all)) {
        # If tuning fails, fall back to a final run with the current settings.
        break
      }

      psi_post <- fit_tune$samples_all[, "psi"]
      M_too_small <- mean(psi_post > psi_threshold, na.rm = TRUE) > psi_prob_cutoff

      rhat_too_high <- FALSE
      if (!is.na(tuning_n_chains) && tuning_n_chains > 1L) {
        chain_samples <- lapply(fit_tune$samples_list, `[[`, "samples")
        quiet_require("MCMCvis")
        diag <- MCMCvis::MCMCsummary(
          chain_samples,
          params = "all",
          ISB = FALSE,
          exact = FALSE,
          n.eff = FALSE
        )
        rhat_val <- suppressWarnings(max(diag$Rhat, na.rm = TRUE))
        rhat_too_high <- !is.na(rhat_val) && (rhat_val > rhat_target)
      }

      # Update parameters for next loop
      if (rhat_too_high) current_iter <- current_iter * 2L
      if (M_too_small) current_M <- current_M * 2L

      if (!rhat_too_high && !M_too_small) break
    }

    final_burnin <- min(current_burnin, as.integer(floor(current_iter / 3)))
    final_burnin <- max(0L, final_burnin)

    return(
      run_USCR(
        out            = out,
        camera_counts  = camera_counts,
        camera_days    = camera_days,
        iter           = current_iter,
        burnin         = final_burnin,
        thin           = thin,
        n_chains       = n_chains,
        M              = current_M,
        log_sigma_mean = log_sigma_mean,
        log_sigma_sd   = log_sigma_sd,
        log_lam0_sd    = log_lam0_sd,
        sd_eps_shape   = sd_eps_shape,
        buffer_chi_p   = buffer_chi_p,
        buffer_scale   = buffer_scale,
        adaptive       = FALSE,
        compute_WAIC   = compute_WAIC,
        diagnostic_mode = diagnostic_mode
      )
    )
  }
  
  quiet_require("nimble")
  
  J <- nrow(out)
  if (length(camera_counts) != J || length(camera_days) != J) {
    stop("run_USCR: camera_counts and camera_days must have length nrow(out).",
         call. = FALSE)
  }
  
  ss <- uscr_state_space_and_area(out, buffer_chi_p = buffer_chi_p,
                                  buffer_scale = buffer_scale)
  
  # J == 1: nimble cannot index min(dist2[i,1:J]) or if(J>1) the same as J>1;
  # use a dedicated single-camera block (Rmd still allows one camera).
  code <- if (J == 1L) {
    nimble::nimbleCode({
      
      log_sigma ~ dnorm(log_sigma_mean, sd = log_sigma_sd)
      log_lam_0 ~ dnorm(0, sd = log_lam0_sd)
      psi       ~ dunif(0, 1)
      sd_eps    ~ dgamma(sd_eps_shape, sd_eps_shape)
      
      sigma <- exp(log_sigma)
      lam_0 <- exp(log_lam_0)
      
      for (i in 1:M) {
        
        z[i] ~ dbern(psi)
        
        hrc[i, 1] ~ dunif(xlim[1], xlim[2])
        hrc[i, 2] ~ dunif(ylim[1], ylim[2])
        
        dist2[i, 1] <- (hrc[i, 1] - cam[1, 1])^2 +
          (hrc[i, 2] - cam[1, 2])^2
        
        lambda[i, 1] <- z[i] * lam_0 *
          exp(-dist2[i, 1] / (2 * sigma^2))
        
        in_ss[i] ~ dconstraint(dist2[i, 1] < buffer_sq)
      }
      
      eps[1] ~ dnorm(0, sd = sd_eps)
      
      Lambda[1] <- sum(lambda[1:M, 1])
      log(mu[1]) <- log(Lambda[1]) + log(day1) + eps[1]
      
      y[1]     ~ dpois(mu[1])
      y_sim[1] ~ dpois(mu[1])
      
      pearson_obs[1] <- (y[1]     - mu[1])^2 / mu[1]
      pearson_sim[1] <- (y_sim[1] - mu[1])^2 / mu[1]
      
      sum_obs <- pearson_obs[1]
      sum_sim <- pearson_sim[1]
      
      bp <- step(sum_sim - sum_obs)
      
      N     <- sum(z[1:M])
      D_mi2 <- N / area_m2
    })
  } else {
    nimble::nimbleCode({
      
      log_sigma ~ dnorm(log_sigma_mean, sd = log_sigma_sd)
      log_lam_0 ~ dnorm(0, sd = log_lam0_sd)
      psi       ~ dunif(0, 1)
      sd_eps    ~ dgamma(sd_eps_shape, sd_eps_shape)
      
      sigma <- exp(log_sigma)
      lam_0 <- exp(log_lam_0)
      
      for (i in 1:M) {
        
        z[i] ~ dbern(psi)
        
        hrc[i, 1] ~ dunif(xlim[1], xlim[2])
        hrc[i, 2] ~ dunif(ylim[1], ylim[2])
        
        for (j in 1:J) {
          dist2[i, j] <- (hrc[i, 1] - cam[j, 1])^2 +
            (hrc[i, 2] - cam[j, 2])^2
          
          lambda[i, j] <- z[i] * lam_0 *
            exp(-dist2[i, j] / (2 * sigma^2))
        }
        
        min_dist2[i] <- min(dist2[i, 1:J])
        in_ss[i] ~ dconstraint(min_dist2[i] < buffer_sq)
      }
      
      for (j in 1:J) {
        
        eps[j] ~ dnorm(0, sd = sd_eps)
        
        Lambda[j] <- sum(lambda[1:M, j])
        log(mu[j]) <- log(Lambda[j]) + log(days_per_cam[j]) + eps[j]
        
        y[j]     ~ dpois(mu[j])
        y_sim[j] ~ dpois(mu[j])
        
        pearson_obs[j] <- (y[j]     - mu[j])^2 / mu[j]
        pearson_sim[j] <- (y_sim[j] - mu[j])^2 / mu[j]
      }
      
      sum_obs <- sum(pearson_obs[1:J])
      sum_sim <- sum(pearson_sim[1:J])
      
      bp <- step(sum_sim - sum_obs)
      
      N     <- sum(z[1:M])
      D_mi2 <- N / area_m2
    })
  }
  
  Const <- list(
    M              = as.integer(M),
    J              = J,
    xlim           = ss$xlim,
    ylim           = ss$ylim,
    buffer_sq      = ss$buffer_sq,
    area_m2        = ss$area_m2,
    cam            = as.matrix(cbind(out$utm_e, out$utm_n)),
    log_sigma_mean = log_sigma_mean,
    log_sigma_sd   = log_sigma_sd,
    log_lam0_sd    = log_lam0_sd,
    sd_eps_shape   = sd_eps_shape
  )
  if (J == 1L) {
    Const$day1 <- as.numeric(camera_days[1])
  } else {
    Const$days_per_cam <- as.numeric(camera_days)
  }
  
  Data <- list(
    y     = as.numeric(camera_counts),
    in_ss = rep(1L, M)
  )
  
  monitors <- c("log_sigma", "log_lam_0", "psi", "sd_eps",
                "sigma", "lam_0", "N", "D_mi2",
                "sum_obs", "sum_sim", "bp")
  if (isTRUE(diagnostic_mode)) {
    # Heavy diagnostics (can be slow / memory intensive)
    monitors <- c(monitors, "z", "hrc", "eps")
  }
  
  quiet_require("parallel")
  
  inits_fun <- function() {
    Mloc <- as.integer(Const$M)
    Jloc <- Const$J
    cam  <- Const$cam
    hrc  <- matrix(nrow = Mloc, ncol = 2)
    idx  <- sample.int(Jloc, Mloc, replace = TRUE)
    hrc[, 1] <- cam[idx, 1] + stats::runif(Mloc, -0.01, 0.01)
    hrc[, 2] <- cam[idx, 2] + stats::runif(Mloc, -0.01, 0.01)
    list(hrc = hrc)
  }
  
  run_MCMC_chain <- function(index, code, constants, data, monitors,
                             niter, nburnin, thin, inits_fun) {
    library(nimble)
    nimble::nimbleMCMC(
      code      = code,
      constants = constants,
      data      = data,
      inits     = inits_fun(),
      monitors  = monitors,
      niter     = niter,
      nburnin   = nburnin,
      thin      = thin,
      WAIC      = compute_WAIC
    )
  }
  
  if (n_chains > 1) {
    cl_size <- min(n_chains, parallel::detectCores() - 1)
    cl_size <- max(1, cl_size)
    cl <- parallel::makeCluster(cl_size)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    samples_list <- parallel::parLapply(
      cl         = cl,
      X          = seq_len(n_chains),
      fun        = run_MCMC_chain,
      code       = code,
      constants  = Const,
      data       = Data,
      monitors   = monitors,
      niter      = iter,
      nburnin    = burnin,
      thin       = thin,
      inits_fun  = inits_fun
    )
  } else {
    samples_list <- list(
      run_MCMC_chain(1, code, Const, Data, monitors, iter, burnin, thin, inits_fun)
    )
  }
  
  samples_all <- do.call(rbind, lapply(samples_list, `[[`, "samples"))
  
  waic_vec <- vapply(
    samples_list,
    function(f) if (!is.null(f$WAIC)) f$WAIC$WAIC else NA_real_,
    numeric(1)
  )
  waic_mean <- mean(waic_vec, na.rm = TRUE)
  
  list(
    method       = "USCR",
    samples_list = samples_list,
    samples_all  = samples_all,
    waic         = waic_mean
  )
}

# -------------------------------------------------------------------
# 6. WAIC-based model averaging — mirrors Rmd logic
# -------------------------------------------------------------------
waic_model_average <- function(rem_fit,
                               tte_fit,
                               uscr_fit,
                               prob_threshold = 20) {
  quiet_require("dplyr")
  quiet_require("tibble")
  
  # Expect samples_all with column D_mi2
  rem_D  <- rem_fit$samples_all[, "D_mi2"]
  tte_D  <- tte_fit$samples_all[, "D_mi2"]
  uscr_D <- uscr_fit$samples_all[, "D_mi2"]
  
  waic_tbl <- tibble::tibble(
    model = c("REM", "TTE", "USCR"),
    waic  = c(rem_fit$waic, tte_fit$waic, uscr_fit$waic)
  ) |>
    dplyr::mutate(
      deltaWAIC = waic - min(waic),
      rel_lik   = exp(-0.5 * deltaWAIC),
      w         = rel_lik / sum(rel_lik)
    ) |>
    dplyr::arrange(model)
  
  density_estimates <- tibble::tibble(
    REM  = rem_D,
    TTE  = tte_D,
    USCR = uscr_D
  ) |>
    dplyr::mutate(
      unweighted_mean = apply(cbind(sort(REM), sort(TTE), sort(USCR)), 1, mean),
      weighted_mean   = apply(
        waic_tbl$w * t(cbind(sort(REM), sort(TTE), sort(USCR))), 2, sum
      )
    )
  
  # Summaries per column (REM, TTE, USCR, unweighted, weighted)
  stats_mat <- apply(
    density_estimates,
    2,
    function(x) c(
      mean   = mean(x),
      lower  = stats::quantile(x, 0.025),
      upper  = stats::quantile(x, 0.975),
      probgt = mean(x > prob_threshold)
    )
  )
  
  table_out <- tibble::tibble(
    Method          = c(waic_tbl$model, "unweighted mean", "weighted mean"),
    `Mean density`  = stats_mat["mean", ],
    `Lower 2.5%`    = stats_mat["lower", ],
    `Upper 97.5%`   = stats_mat["upper", ],
    `Prob > 20 DPSM`= stats_mat["probgt", ],
    `WAIC weight`   = c(waic_tbl$w, NA, NA)
  )
  
  list(
    waic_table    = waic_tbl,
    draws_table   = density_estimates,
    summary_table = table_out
  )
}

