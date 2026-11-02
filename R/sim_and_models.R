## R/sim_and_models.R
## Simulation helpers + NPS model inputs + REM / USCR / TTE models
## Models follow Camera_Unmarked_Analyses_updated.Rmd, with a speed-aware adaptive USCR fit for app use.

# -------------------------------------------------------------------
# Package checks (loaded lazily where needed)
# -------------------------------------------------------------------

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }
}

# -------------------------------------------------------------------
# 1. Simulation helpers (for teaching/sim tab)
# -------------------------------------------------------------------
home_range_km2_to_sigma_m <- function(home_range_km2) {
  radius_km <- sqrt(home_range_km2 / pi)
  sigma_km <- radius_km / 2.45
  sigma_km * 1000
}

tte_detection_radius_multiplier <- function(theta_deg) {
  theta_deg <- as.numeric(theta_deg)
  3.324e-01 + 5.580e-03 * theta_deg - 1.454e-05 * theta_deg^2
}

simulate_camera_counts <- function(n_side   = 5,
                                   spacing_m = 300,
                                   days      = 21,
                                   D_per_km2 = 25,
                                   lambda0   = 0.20,
                                   home_range_km2 = 0.89,
                                   sigma_m   = NULL,
                                   buffer_m  = NULL,
                                   seed      = 1) {
  quiet_require("secr")
  
  set.seed(seed)
  if (is.null(sigma_m)) {
    sigma_m <- home_range_km2_to_sigma_m(home_range_km2)
  }
  if (is.null(buffer_m)) {
    buffer_m <- 4 * sigma_m
  }
  
  traps_obj <- secr::make.grid(
    nx       = n_side,
    ny       = n_side,
    spacing  = spacing_m,
    detector = "count"
  )
  
  mask_obj <- secr::make.mask(
    traps   = traps_obj,
    buffer  = buffer_m,
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
      home_range_km2 = home_range_km2,
      sigma_m   = sigma_m,
      buffer_m  = buffer_m,
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
      `Site Name`       = paste0("C", dplyr::row_number()),
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

build_teaching_sim_grid <- function(n_side, spacing_m, detection_radius_m, days) {
  coords <- expand.grid(
    x = seq(0, by = spacing_m, length.out = n_side),
    y = seq(0, by = spacing_m, length.out = n_side)
  )
  coords <- coords[order(coords$y, coords$x), , drop = FALSE]
  coords$x <- coords$x - mean(range(coords$x))
  coords$y <- coords$y - mean(range(coords$y))
  
  coords |>
    dplyr::mutate(
      `Site Name` = paste0("C", dplyr::row_number()),
      utm_e = x / 1000,
      utm_n = y / 1000,
      `Detection Distance` = detection_radius_m,
      `Start Index` = 1L,
      `End Index` = as.integer(days)
    )
}

simulate_teaching_counts <- function(model = c("REM", "TTE"),
                                     n_side = 5,
                                     spacing_m = 300,
                                     days = 21,
                                     D_per_km2 = 25,
                                     detection_radius_m = 12,
                                     theta_deg = 55,
                                     v_km_day = 4,
                                     sd_eps = 0.2,
                                     seed = 1) {
  model <- match.arg(model)
  quiet_require("dplyr")
  
  if (n_side < 1 || days < 1 || spacing_m <= 0 || detection_radius_m <= 0) {
    stop("Teaching simulator inputs must be positive.", call. = FALSE)
  }
  if (D_per_km2 < 0 || v_km_day <= 0 || sd_eps < 0) {
    stop("Density, movement speed, and heterogeneity inputs are out of range.", call. = FALSE)
  }
  
  set.seed(seed)
  
  out <- build_teaching_sim_grid(
    n_side = n_side,
    spacing_m = spacing_m,
    detection_radius_m = detection_radius_m,
    days = days
  )
  
  J <- nrow(out)
  camera_days <- rep(as.numeric(days), J)
  r_km <- rep(detection_radius_m / 1000, J)
  eps <- stats::rnorm(J, mean = 0, sd = sd_eps)
  
  if (identical(model, "REM")) {
    theta_rad <- theta_deg * pi / 180
    lambda <- ((2 + theta_rad) / pi) * v_km_day * r_km * D_per_km2 * camera_days * exp(eps)
  } else {
    area_km2 <- pi * r_km^2 * theta_deg / 360
    time_unit <- tte_detection_radius_multiplier(theta_deg) * r_km / v_km_day
    tte_units <- camera_days / time_unit
    lambda <- D_per_km2 * tte_units * area_km2 * exp(eps)
  }
  
  camera_counts <- stats::rpois(J, lambda = lambda)
  
  list(
    model = model,
    out = out,
    camera_counts = camera_counts,
    camera_days = camera_days,
    truth = list(
      model = model,
      D_per_km2 = D_per_km2,
      D_per_mi2 = D_per_km2 * 2.59,
      detection_radius_m = detection_radius_m,
      theta_deg = theta_deg,
      v_km_day = v_km_day,
      sd_eps = sd_eps,
      seed = seed,
      n_cams = J,
      n_side = n_side,
      spacing_m = spacing_m,
      days = days,
      mean_lambda = mean(lambda),
      total_expected = sum(lambda),
      total_observed = sum(camera_counts)
    )
  )
}

# Back-compatible alias used by older app code
build_sim_data_for_nimble <- function(ch, detection_radius_m) {
  sim_model_inputs(
    sim = list(ch = ch),
    detection_radius_m = detection_radius_m
  )
}

# -------------------------------------------------------------------
# 2. NPS inputs: format deployments + images as in the Rmd
# -------------------------------------------------------------------
build_nps_model_inputs <- function(
    deployments,
    images,
    species_to_analyze,
    max_days = NULL
    ) {
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
  
  deps <- format_deployments(deployments, max_days = max_days)

  # Timestamp is already parsed to POSIXct in check_images() before
  # images_checked() is ever set, so no re-parsing is needed here.

  # Latitude/Longitude presence and validity are already enforced as a hard
  # stop in check_deployments() at upload time, so all that's needed here is
  # numeric coercion (deployment$Longitude/Latitude may still be character
  # columns coming out of CSV import even though every value parses cleanly).
  deps <- deps |>
    dplyr::mutate(
      Longitude = as.numeric(Longitude),
      Latitude  = as.numeric(Latitude)
    )

  # UTM conversion and centering. Longitude/Latitude are treated as WGS84
  # (EPSG:4326, the standard GPS datum) and projected into the matching
  # WGS84 UTM zone, so this works for deployments anywhere in the world
  # rather than assuming North America. The UTM zone comes from the mean
  # longitude and the hemisphere (for choosing the Northern- vs.
  # Southern-Hemisphere EPSG series) comes from the mean latitude; as
  # before, this assumes all cameras in a single dataset fall within one
  # UTM zone, which holds for a spatially clustered deployment/array.
  mean_lon  <- mean(deps$Longitude, na.rm = TRUE)
  mean_lat  <- mean(deps$Latitude, na.rm = TRUE)
  utm_zone  <- floor((mean_lon + 180) / 6) + 1
  utm_zone  <- min(60L, max(1L, as.integer(utm_zone)))
  epsg_code <- if (mean_lat >= 0) 32600 + utm_zone else 32700 + utm_zone

  utm_coords <- deps |>
    sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) |>
    sf::st_transform(crs = epsg_code) |>
    sf::st_coordinates()
  
  deps <- deps |>
    dplyr::mutate(
      utm_e = (utm_coords[, "X"] -
                 mean(range(utm_coords[, "X"], na.rm = TRUE))) / 1000,
      utm_n = (utm_coords[, "Y"] -
                 mean(range(utm_coords[, "Y"], na.rm = TRUE))) / 1000
    )
  
  # Grouping by Site Name and Cluster ID together (rather than Cluster ID
  # alone) means a Cluster ID that happens to appear at more than one site
  # (e.g. an event-numbering scheme that isn't actually unique per camera)
  # produces separate groups per site instead of silently merging their
  # detections under one site. check_images() also flags this case directly
  # at upload time.
  seqs <- images |>
    dplyr::filter(Species == species_to_analyze) |>
    dplyr::group_by(`Site Name`, `Cluster ID`) |>
    dplyr::summarise(
      detection_date = lubridate::date(min(Timestamp)),
      species_count  = max(as.numeric(`Sighting Count`), na.rm = TRUE),
      start          = min(Timestamp),
      end            = max(Timestamp),
      .groups        = "drop"
    ) |>
    dplyr::mutate(
      cluster_length_days = as.numeric(difftime(end, start, units = "days"))
    )
  
  # Parse Start/End Date into real Date columns, carried on deps from here
  # on (through the join into `out` below) so they only need to be parsed
  # once. A plain as.Date() with a fixed format is sufficient here: format_
  # deployments() always re-normalizes End Date to mm/dd/yyyy, and
  # check_deployments() now hard-stops any Start Date/End Date that isn't
  # already in that format, so nothing reaches this function needing a
  # more permissive multi-format parser.
  deps <- deps |>
    dplyr::mutate(
      Start_Date_parsed = as.Date(`Start Date`, format = "%m/%d/%Y"),
      End_Date_parsed   = as.Date(`End Date`,   format = "%m/%d/%Y")
    )


  min_start <- min(deps$Start_Date_parsed, na.rm = TRUE)
  max_end   <- max(deps$End_Date_parsed, na.rm = TRUE)

  # Zero-count rows across all sites × all days in deployment window
  zero_counts <- expand.grid(
    `Site Name` = unique(deps$`Site Name`),
    detection_date = seq(
      from = min_start,
      to   = max_end,
      by   = "day"
    )
  ) |>
    dplyr::mutate(
      species_count       = 0,
      cluster_length_days = 0
    )

  counts_time <- seqs |>
    dplyr::bind_rows(zero_counts) |>
    dplyr::group_by(`Site Name`, detection_date) |>
    dplyr::summarise(
      detections       = sum(species_count),
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
    dplyr::left_join(deps, by = "Site Name")

  detection_matrix <- counts |>
    tibble::column_to_rownames("Site Name") |>
    as.matrix()

  camera_time_matrix <- counts_time |>
    dplyr::select(-detections) |>
    tidyr::pivot_wider(
      names_from  = detection_date,
      values_from = camera_time_days,
      values_fill = 0,
      names_sort  = TRUE
    ) |>
    tibble::column_to_rownames("Site Name") |>
    as.matrix()
  
  # Create Start Index and End Index columns
  # These map each camera's deployment dates to column indices in the detection matrix.
  # Derived directly from detection_matrix's own column names (rather than an
  # independently-sorted copy of counts_time$detection_date) so all_dates[k]
  # is guaranteed to be the date of column k by construction, not by relying
  # on pivot_wider's string-based names_sort happening to agree with a
  # separate chronological Date sort.
  all_dates <- as.Date(colnames(detection_matrix))

  # Start_Date_parsed/End_Date_parsed already arrived on `out` via the
  # left_join(deps, ...) above, so no re-parsing is needed here.
  out <- out |>
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
    out           = out,
    camera_counts = camera_counts,
    camera_days   = camera_days
  )
}

# -------------------------------------------------------------------
# 3. REM model (nimble) — mirrors Rmd structure
# -------------------------------------------------------------------
run_basic_nimble_chain <- function(code,
                                   constants,
                                   data,
                                   monitors,
                                   niter,
                                   nburnin,
                                   thin,
                                   compute_WAIC = TRUE) {
  quiet_require("nimble")

  out <- nimble::nimbleMCMC(
    code = code,
    constants = constants,
    data = data,
    monitors = monitors,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    WAIC = compute_WAIC,
    summary = FALSE,
    samplesAsCodaMCMC = FALSE
  )

  normalize_nimble_output(out)
}

run_basic_nimble_chains <- function(code,
                                    constants,
                                    data,
                                    monitors,
                                    niter,
                                    nburnin,
                                    thin,
                                    n_chains = 1,
                                    parallel_chains = TRUE,
                                    compute_WAIC = TRUE,
                                    seed = NULL) {
  quiet_require("parallel")

  if (n_chains <= 1L || !parallel_chains) {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    return(list(
      run_basic_nimble_chain(
        code = code,
        constants = constants,
        data = data,
        monitors = monitors,
        niter = niter,
        nburnin = nburnin,
        thin = thin,
        compute_WAIC = compute_WAIC
      )
    ))
  }

  cl_size <- min(n_chains, max(1L, parallel::detectCores() - 1L))
  cl <- parallel::makeCluster(cl_size)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  fits <- parallel::parLapply(
    cl = cl,
    X = seq_len(n_chains),
    fun = function(chain_id, code, constants, data, monitors,
                   niter, nburnin, thin, compute_WAIC) {
      library(nimble)

      out <- nimble::nimbleMCMC(
        code = code,
        constants = constants,
        data = data,
        monitors = monitors,
        niter = niter,
        nburnin = nburnin,
        thin = thin,
        WAIC = compute_WAIC,
        summary = FALSE,
        samplesAsCodaMCMC = FALSE
      )

      if (is.matrix(out) || is.data.frame(out)) {
        return(list(samples = as.matrix(out), WAIC = NULL))
      }
      out
    },
    code = code,
    constants = constants,
    data = data,
    monitors = monitors,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    compute_WAIC = compute_WAIC
  )

  fits
}

run_encounter_rate_model <- function(method,
                                     code,
                                     constants,
                                     data,
                                     monitors,
                                     iter,
                                     burnin,
                                     thin,
                                     n_chains,
                                     adaptive = TRUE,
                                     compute_WAIC = TRUE,
                                     rhat_target = 1.1,
                                     max_adapt_rounds = NULL,
                                     parallel_chains = TRUE,
                                     status_callback = NULL,
                                     seed = NULL,
                                     verbose = FALSE) {
  report_status <- function(round, iter, burnin, prev_elapsed_min = NULL, prev_iter = NULL) {
    if (is.function(status_callback)) {
      status_callback(
        round = round,
        iter = iter,
        burnin = burnin,
        prev_elapsed_min = prev_elapsed_min,
        prev_iter = prev_iter
      )
    }
  }

  iter <- as.integer(iter)
  burnin <- as.integer(burnin)
  thin <- as.integer(max(1, thin))
  n_chains <- as.integer(max(1, n_chains))

  current_iter <- iter
  current_burnin <- burnin
  current_thin <- thin
  round_i <- 0L
  current_rhat <- Inf
  adapt_log <- list()
  last_fit <- NULL
  prev_elapsed_min <- NULL
  prev_iter <- NULL

  repeat {
    round_i <- round_i + 1L

    if (!isTRUE(adaptive) && round_i > 1L) {
      break
    }

    if (!is.null(max_adapt_rounds) && round_i > as.integer(max_adapt_rounds)) {
      warning(
        method,
        " adaptive tuning reached max_adapt_rounds = ",
        max_adapt_rounds,
        " before convergence. Proceeding with the latest fit."
      )
      break
    }

    if (isTRUE(verbose)) {
      message(
        method, " MCMC round ", round_i,
        ": iter=", current_iter,
        ", burnin=", current_burnin,
        ", thin=", current_thin,
        ", chains=", n_chains
      )
    }
    report_status(
      round = round_i,
      iter = current_iter,
      burnin = current_burnin,
      prev_elapsed_min = prev_elapsed_min,
      prev_iter = prev_iter
    )
    round_started <- Sys.time()

    fit <- run_basic_nimble_chains(
      code = code,
      constants = constants,
      data = data,
      monitors = monitors,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      parallel_chains = parallel_chains,
      compute_WAIC = compute_WAIC,
      seed = if (is.null(seed)) NULL else seed + round_i
    )

    samples_list <- lapply(fit, `[[`, "samples")
    current_rhat <- safe_rhat_max(samples_list)
    converged <- is.na(current_rhat) || current_rhat <= rhat_target

    adapt_log[[round_i]] <- data.frame(
      round = round_i,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      rhat_max = current_rhat,
      stringsAsFactors = FALSE
    )

    last_fit <- fit
    prev_elapsed_min <- as.numeric(difftime(Sys.time(), round_started, units = "mins"))
    prev_iter <- current_iter

    if (converged || !isTRUE(adaptive)) {
      break
    }

    current_iter <- current_iter * 2L
    current_burnin <- current_burnin * 2L
    current_thin <- current_thin * 2L
  }

  samples_list <- lapply(last_fit, `[[`, "samples")
  samples_all <- do.call(rbind, samples_list)

  list(
    method = method,
    samples_list = samples_list,
    samples_all = samples_all,
    waic = extract_waic_mean(last_fit),
    fit_objects = last_fit,
    round_history = if (length(adapt_log) > 0L) do.call(rbind, adapt_log) else NULL,
    final_rhat_max = current_rhat,
    settings = list(
      iter = current_iter,
      burnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      adaptive = adaptive,
      max_adapt_rounds = max_adapt_rounds,
      rhat_target = rhat_target
    )
  )
}

run_REM <- function(y,
                    r_km,
                    camera_days,
                    iter     = 6000,
                    burnin   = 1000,
                    thin     = 5,
                    n_chains = 2,
                    D_max    = 200,
                    log_v_mean = 1.130,
                    log_v_sd   = 0.3372,
                    sd_eps_shape = 1,
                    sd_eps_rate  = 1,
                    theta_deg  = 55,
                    adaptive = TRUE,
                    compute_WAIC = TRUE,
                    rhat_target = 1.1,
                    max_adapt_rounds = NULL,
                    parallel_chains = TRUE,
                    status_callback = NULL,
                    seed = NULL,
                    verbose = FALSE) {
  
  quiet_require("nimble")
  
  J <- length(y)
  if (length(r_km) != J || length(camera_days) != J) {
    stop("run_REM: y, r_km, and camera_days must have same length.", call. = FALSE)
  }
  if (J < 2L) {
    stop("run_REM: at least 2 cameras are required to estimate density.", call. = FALSE)
  }
  if (length(theta_deg) == 1L) {
    theta_deg <- rep(as.numeric(theta_deg), J)
  } else if (length(theta_deg) != J) {
    stop("run_REM: theta_deg must have length 1 or the same length as y.", call. = FALSE)
  }
  
  code <- nimble::nimbleCode({
    
    # Priors
    D      ~ dunif(0, D_max)
    log_v  ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dgamma(sd_eps_shape, sd_eps_rate)
    
    v     <- exp(log_v)
    
    for (j in 1:J) {
      
      eps[j] ~ dnorm(0, sd = sd_eps)
      
      log(lambda[j]) <- log(2 + theta[j]) - log(pi_const) + log(v) +
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
    theta      = as.numeric(theta_deg) * pi / 180,
    D_max      = D_max,
    log_v_mean = log_v_mean,
    log_v_sd   = log_v_sd,
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate  = sd_eps_rate
  )
  
  Data <- list(y = as.numeric(y))
  
  monitors <- c("D", "v", "sd_eps",
                "bp", "sum_obs", "sum_sim", "D_mi2")

  run_encounter_rate_model(
    method = "REM",
    code = code,
    constants = Const,
    data = Data,
    monitors = monitors,
    iter = iter,
    burnin = burnin,
    thin = thin,
    n_chains = n_chains,
    adaptive = adaptive,
    compute_WAIC = compute_WAIC,
    rhat_target = rhat_target,
    max_adapt_rounds = max_adapt_rounds,
    parallel_chains = parallel_chains,
    status_callback = status_callback,
    seed = seed,
    verbose = verbose
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
                    n_chains = 2,
                    D_max    = 200,
                    log_v_mean = 1.130,
                    log_v_sd   = 0.3372,
                    sd_eps_shape = 1,
                    sd_eps_rate  = 1,
                    theta_deg  = 55,
                    adaptive = TRUE,
                    compute_WAIC = TRUE,
                    rhat_target = 1.1,
                    max_adapt_rounds = NULL,
                    parallel_chains = TRUE,
                    status_callback = NULL,
                    seed = NULL,
                    verbose = FALSE) {
  
  quiet_require("nimble")
  
  J <- length(y)
  if (length(r_km) != J || length(camera_days) != J) {
    stop("run_TTE: y, r_km, and camera_days must have same length.", call. = FALSE)
  }
  if (J < 2L) {
    stop("run_TTE: at least 2 cameras are required to estimate density.", call. = FALSE)
  }
  if (length(theta_deg) == 1L) {
    theta_deg <- rep(as.numeric(theta_deg), J)
  } else if (length(theta_deg) != J) {
    stop("run_TTE: theta_deg must have length 1 or the same length as y.", call. = FALSE)
  }
  
  code <- nimble::nimbleCode({
    
    D      ~ dunif(0, D_max)
    log_v  ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dgamma(sd_eps_shape, sd_eps_rate)
    
    v     <- exp(log_v)
    
    for (j in 1:J) {
      
      a[j] <- pi_const * r[j]^2 * theta[j] / 360
      mvw[j] <- angle_mult[j] * r[j]
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
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate  = sd_eps_rate,
    theta       = as.numeric(theta_deg),
    angle_mult  = tte_detection_radius_multiplier(theta_deg)
  )
  
  Data <- list(y = as.numeric(y))
  
  monitors <- c("D", "v", "sd_eps",
                "bp", "sum_obs", "sum_sim", "D_mi2")

  run_encounter_rate_model(
    method = "TTE",
    code = code,
    constants = Const,
    data = Data,
    monitors = monitors,
    iter = iter,
    burnin = burnin,
    thin = thin,
    n_chains = n_chains,
    adaptive = adaptive,
    compute_WAIC = compute_WAIC,
    rhat_target = rhat_target,
    max_adapt_rounds = max_adapt_rounds,
    parallel_chains = parallel_chains,
    status_callback = status_callback,
    seed = seed,
    verbose = verbose
  )
}

# -------------------------------------------------------------------
# USCR: buffer sizing from sigma (movement/space-use scale)
# -------------------------------------------------------------------
# 99% circular home-range radius (km) implied by a given sigma (km).
uscr_home_range_radius_km <- function(sigma_km, p = 0.99) {
  sigma_km * sqrt(stats::qchisq(p, df = 2))
}

# Suggested state-space buffer (m), given a sigma value (km) and a safety
# multiplier on the home-range radius. Used both to size the buffer up front
# (from the log(sigma) prior) and to check it after fitting (from the
# posterior). A bigger multiplier is safer but increases the state space
# (and run time).
uscr_buffer_from_sigma_km <- function(sigma_km, p = 0.99, multiplier = 1.05) {
  multiplier * uscr_home_range_radius_km(sigma_km, p = p) * 1000
}

# Default buffer (m) computed from the log(sigma) prior before any data are
# fit: uses the p-quantile of the prior as a quick, conservative guess at
# sigma, then sizes the buffer from that guess via uscr_buffer_from_sigma_km().
uscr_default_buffer_m <- function(log_sigma_mean, log_sigma_sd, p = 0.99, multiplier = 1.05) {
  sigma_q <- exp(log_sigma_mean + log_sigma_sd * stats::qnorm(p))
  uscr_buffer_from_sigma_km(sigma_q, p = p, multiplier = multiplier)
}

# -------------------------------------------------------------------
# USCR: state space, buffer, area
# Drop-in replacement for the USCR section in sim_and_models.R
# -------------------------------------------------------------------

uscr_state_space_and_area <- function(out,
                                      buffer_m = NULL,
                                      log_sigma_mean = -1.5269,
                                      log_sigma_sd = 0.1535) {
  if (!is.null(buffer_m) && is.finite(buffer_m) && buffer_m > 0) {
    buffer <- as.numeric(buffer_m) / 1000
  } else {
    buffer <- uscr_default_buffer_m(log_sigma_mean, log_sigma_sd) / 1000
  }
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

    site_col <- if ("Site Name" %in% names(out)) "Site Name" else NULL
    coords_df <- unique(out[, c(site_col, "Longitude", "Latitude"), drop = FALSE])
    # Latitude/Longitude validity is enforced at upload time (see
    # check_deployments()); simulated grids never reach this branch, since
    # they have no Longitude/Latitude columns at all. All that's needed here
    # is numeric coercion.
    coords_df$Longitude <- as.numeric(coords_df$Longitude)
    coords_df$Latitude  <- as.numeric(coords_df$Latitude)
    # Same WGS84/UTM approach as build_nps_model_inputs(): global rather
    # than North-America-only, with the hemisphere taken from mean latitude.
    mean_lon <- mean(coords_df$Longitude, na.rm = TRUE)
    mean_lat <- mean(coords_df$Latitude, na.rm = TRUE)
    utm_zone <- floor((mean_lon + 180) / 6) + 1
    utm_zone <- min(60L, max(1L, as.integer(utm_zone)))
    epsg_code <- if (mean_lat >= 0) 32600 + utm_zone else 32700 + utm_zone

    cam_buff <- sf::st_as_sf(coords_df, coords = c("Longitude", "Latitude"), crs = 4326) |>
      sf::st_transform(crs = epsg_code) |>
      sf::st_buffer(buffer * 1000) |>
      sf::st_union()

    area_mi2 <- as.numeric(sf::st_area(cam_buff)) / (2.59 * 1e6)
  } else {
    if (any(!is.finite(out$utm_e) | !is.finite(out$utm_n))) {
      bad <- which(!is.finite(out$utm_e) | !is.finite(out$utm_n))
      bad_sites <- if ("Site Name" %in% names(out)) as.character(out$`Site Name`[bad]) else as.character(bad)
      stop(
        "Missing or invalid projected coordinates for site(s): ",
        paste(unique(bad_sites), collapse = ", "),
        call. = FALSE
      )
    }
    # Simulated UTM-only grids: preserve current script behavior
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
# USCR model code builder
# -------------------------------------------------------------------

# USCR requires at least 2 cameras: sigma (the spatial decay of detection
# with distance) is estimated from spatial contrast in detections across
# cameras at different locations, so a single camera can't identify it.
# This is the only supported model form; run_USCR() enforces J >= 2 before
# calling this.
build_uscr_code <- function() {
  quiet_require("nimble")

  nimble::nimbleCode({
    log_sigma ~ dnorm(log_sigma_mean, sd = log_sigma_sd)
    log_lam_0 ~ dnorm(log_lam0_mean, sd = log_lam0_sd)
    psi ~ dunif(0, 1)
    sd_eps ~ dgamma(sd_eps_shape, sd_eps_rate)

    sigma <- exp(log_sigma)
    lam_0 <- exp(log_lam_0)

    for (i in 1:M) {
      z[i] ~ dbern(psi)

      hrc[i, 1] ~ dunif(xlim[1], xlim[2])
      hrc[i, 2] ~ dunif(ylim[1], ylim[2])

      for (j in 1:J) {
        dist2[i, j] <- (hrc[i, 1] - cam[j, 1]) ^ 2 +
          (hrc[i, 2] - cam[j, 2]) ^ 2

        lambda[i, j] <- z[i] * lam_0 *
          exp(-dist2[i, j] / (2 * sigma ^ 2))
      }

      min_dist2[i] <- min(dist2[i, 1:J])
      in_ss[i] ~ dconstraint(min_dist2[i] < buffer_sq)
    }

    for (j in 1:J) {
      eps[j] ~ dnorm(0, sd = sd_eps)
      Lambda[j] <- sum(lambda[1:M, j])
      log(mu[j]) <- log(Lambda[j]) + log(days_per_cam[j]) + eps[j]

      y[j] ~ dpois(mu[j])
      y_sim[j] ~ dpois(mu[j])

      pearson_obs[j] <- (y[j] - mu[j]) ^ 2 / mu[j]
      pearson_sim[j] <- (y_sim[j] - mu[j]) ^ 2 / mu[j]
    }

    sum_obs <- sum(pearson_obs[1:J])
    sum_sim <- sum(pearson_sim[1:J])
    bp <- step(sum_sim - sum_obs)

    N <- sum(z[1:M])
    D_mi2 <- N / area_mi2
  })
}

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

make_uscr_constants <- function(out,
                                camera_days,
                                M,
                                log_sigma_mean,
                                log_sigma_sd,
                                log_lam0_mean,
                                log_lam0_sd,
                                sd_eps_shape,
                                sd_eps_rate,
                                buffer_m) {
  ss <- uscr_state_space_and_area(
    out,
    buffer_m = buffer_m,
    log_sigma_mean = log_sigma_mean,
    log_sigma_sd = log_sigma_sd
  )

  const <- list(
    M = as.integer(M),
    J = nrow(out),
    xlim = ss$xlim,
    ylim = ss$ylim,
    buffer_sq = ss$buffer_sq,
    area_mi2 = ss$area_mi2,
    cam = as.matrix(cbind(out$utm_e, out$utm_n)),
    log_sigma_mean = log_sigma_mean,
    log_sigma_sd = log_sigma_sd,
    log_lam0_mean = log_lam0_mean,
    log_lam0_sd = log_lam0_sd,
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate = sd_eps_rate
  )

  const$days_per_cam <- as.numeric(camera_days)

  const
}

make_uscr_inits <- function(constants) {
  Mloc <- as.integer(constants$M)
  Jloc <- as.integer(constants$J)
  cam <- constants$cam

  idx <- sample.int(Jloc, Mloc, replace = TRUE)
  hrc <- cbind(
    cam[idx, 1] + stats::runif(Mloc, -0.01, 0.01),
    cam[idx, 2] + stats::runif(Mloc, -0.01, 0.01)
  )

  list(hrc = hrc)
}

normalize_nimble_output <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(list(samples = as.matrix(x), WAIC = NULL))
  }
  x
}

run_uscr_chain <- function(code,
                           constants,
                           data,
                           monitors,
                           niter,
                           nburnin,
                           thin,
                           compute_WAIC = FALSE) {
  quiet_require("nimble")

  out <- nimble::nimbleMCMC(
    code = code,
    constants = constants,
    data = data,
    inits = make_uscr_inits(constants),
    monitors = monitors,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    WAIC = compute_WAIC,
    summary = FALSE,
    samplesAsCodaMCMC = FALSE
  )

  normalize_nimble_output(out)
}

run_uscr_chains <- function(code,
                            constants,
                            data,
                            monitors,
                            niter,
                            nburnin,
                            thin,
                            n_chains = 1,
                            parallel_chains = TRUE,
                            compute_WAIC = FALSE,
                            seed = NULL) {
  quiet_require("parallel")

  if (n_chains <= 1L || !parallel_chains) {
    if (!is.null(seed)) {
      set.seed(seed)
    }

    return(list(
      run_uscr_chain(
        code = code,
        constants = constants,
        data = data,
        monitors = monitors,
        niter = niter,
        nburnin = nburnin,
        thin = thin,
        compute_WAIC = compute_WAIC
      )
    ))
  }

  cl_size <- min(n_chains, max(1L, parallel::detectCores() - 1L))
  cl <- parallel::makeCluster(cl_size)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  fits <- parallel::parLapply(
    cl = cl,
    X = seq_len(n_chains),
    fun = function(chain_id, code, constants, data, monitors,
                   niter, nburnin, thin, compute_WAIC) {
      library(nimble)

      Mloc <- as.integer(constants$M)
      Jloc <- as.integer(constants$J)
      cam <- constants$cam
      idx <- sample.int(Jloc, Mloc, replace = TRUE)
      hrc <- cbind(
        cam[idx, 1] + runif(Mloc, -0.01, 0.01),
        cam[idx, 2] + runif(Mloc, -0.01, 0.01)
      )

      out <- nimble::nimbleMCMC(
        code = code,
        constants = constants,
        data = data,
        inits = list(hrc = hrc),
        monitors = monitors,
        niter = niter,
        nburnin = nburnin,
        thin = thin,
        WAIC = compute_WAIC,
        summary = FALSE,
        samplesAsCodaMCMC = FALSE
      )

      if (is.matrix(out) || is.data.frame(out)) {
        return(list(samples = as.matrix(out), WAIC = NULL))
      }
      out
    },
    code = code,
    constants = constants,
    data = data,
    monitors = monitors,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    compute_WAIC = compute_WAIC
  )

  fits
}

safe_rhat_max <- function(samples_list) {
  if (length(samples_list) < 2L) {
    return(NA_real_)
  }

  quiet_require("MCMCvis")
  diag <- tryCatch(
    MCMCvis::MCMCsummary(samples_list),
    error = function(e) NULL
  )

  if (is.null(diag) || !"Rhat" %in% names(diag)) {
    return(NA_real_)
  }

  rhat <- suppressWarnings(max(diag$Rhat, na.rm = TRUE))
  if (is.infinite(rhat)) {
    NA_real_
  } else {
    rhat
  }
}

extract_waic_mean <- function(fits) {
  waic_vec <- vapply(
    fits,
    function(x) {
      if (!is.null(x$WAIC) && !is.null(x$WAIC$WAIC)) {
        x$WAIC$WAIC
      } else {
        NA_real_
      }
    },
    numeric(1)
  )

  if (all(is.na(waic_vec))) {
    NA_real_
  } else {
    mean(waic_vec, na.rm = TRUE)
  }
}

# -------------------------------------------------------------------
# Power-analysis helpers: compare posterior draws to known simulated truth
# -------------------------------------------------------------------

# Pull posterior draws for a given parameter name, normalizing density to
# D_km2 regardless of whether the fit monitors D (REM/TTE) or D_mi2 (USCR).
extract_param_draws <- function(samples_all, param_name) {
  cn <- colnames(samples_all)
  if (identical(param_name, "D_km2")) {
    if ("D" %in% cn) return(as.numeric(samples_all[, "D"]))
    if ("D_mi2" %in% cn) return(as.numeric(samples_all[, "D_mi2"]) / 2.59)
    return(NULL)
  }
  if (param_name %in% cn) return(as.numeric(samples_all[, param_name]))
  NULL
}

# Build a truth-vs-estimate table for the power-analysis blocks.
# truth: named list, e.g. list(D_km2 = 25, sigma = 0.217, lam_0 = 0.2)
summarize_vs_truth <- function(fit, truth) {
  if (is.null(fit) || is.null(fit$samples_all)) return(NULL)
  
  rows <- lapply(names(truth), function(param_name) {
    true_val <- truth[[param_name]]
    if (is.null(true_val) || !is.finite(true_val)) return(NULL)
    
    draws <- extract_param_draws(fit$samples_all, param_name)
    if (is.null(draws) || !length(draws)) return(NULL)
    
    q <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      Parameter    = param_name,
      Truth        = true_val,
      Mean         = mean(draws, na.rm = TRUE),
      `Lower 2.5%` = q[[1]],
      `Upper 97.5%`= q[[2]],
      Covered      = true_val >= q[[1]] && true_val <= q[[2]],
      check.names  = FALSE
    )
  })
  
  do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
}

# -------------------------------------------------------------------
# Main USCR wrapper
# Direct app implementation of the collaborator Rmd while-loop logic.
# -------------------------------------------------------------------

run_USCR <- function(out,
                     camera_counts,
                     camera_days,
                     iter = 6000,
                     burnin = 1000,
                     thin = 5,
                     n_chains = 2,
                     M = 100,
                     log_sigma_mean = -1.5269,
                     log_sigma_sd = 0.1535,
                     log_lam0_mean = 0,
                     log_lam0_sd = 1,
                     sd_eps_shape = 1,
                     sd_eps_rate = 1,
                     buffer_m = NULL,
                     adaptive = TRUE,
                     compute_WAIC = TRUE,
                     diagnostic_mode = FALSE,
                     rhat_target = 1.1,
                     psi_threshold = 0.9,
                     psi_prob_cutoff = 0.01,
                     buffer_p = 0.99,
                     buffer_multiplier = 1.05,
                     buffer_safety_factor = 1.01,
                     max_adapt_rounds = NULL,
                     iter_cap = NULL,
                     M_cap = NULL,
                     parallel_chains = TRUE,
                     status_callback = NULL,
                     seed = NULL,
                     verbose = FALSE) {

  quiet_require("nimble")

  report_status <- function(round, iter, burnin, M = NULL, buffer_m = NULL, prev_elapsed_min = NULL, prev_iter = NULL) {
    if (is.function(status_callback)) {
      status_callback(
        round = round,
        iter = iter,
        burnin = burnin,
        M = M,
        buffer_m = buffer_m,
        prev_elapsed_min = prev_elapsed_min,
        prev_iter = prev_iter
      )
    }
  }

  build_uscr_result <- function(fit_objects,
                                samples_list,
                                samples_all,
                                round_history,
                                final_rhat_max,
                                settings_M,
                                settings_iter,
                                settings_burnin,
                                settings_thin,
                                settings_buffer_m,
                                area_mi2 = NA_real_) {
    list(
      method = "USCR",
      samples_list = samples_list,
      samples_all = samples_all,
      waic = extract_waic_mean(fit_objects),
      fit_objects = fit_objects,
      round_history = round_history,
      final_rhat_max = final_rhat_max,
      settings = list(
        M = settings_M,
        iter = settings_iter,
        burnin = settings_burnin,
        thin = settings_thin,
        iter_cap = iter_cap,
        M_cap = M_cap,
        n_chains = n_chains,
        compute_WAIC = compute_WAIC,
        diagnostic_mode = diagnostic_mode,
        adaptive = adaptive,
        max_adapt_rounds = max_adapt_rounds,
        rhat_target = rhat_target,
        psi_threshold = psi_threshold,
        psi_prob_cutoff = psi_prob_cutoff,
        buffer_m = settings_buffer_m,
        buffer_p = buffer_p,
        buffer_multiplier = buffer_multiplier,
        area_mi2 = area_mi2,
        log_sigma_mean = log_sigma_mean,
        log_sigma_sd = log_sigma_sd,
        log_lam0_mean = log_lam0_mean,
        log_lam0_sd = log_lam0_sd,
        sd_eps_shape = sd_eps_shape,
        sd_eps_rate = sd_eps_rate
      )
    )
  }

  J <- nrow(out)
  if (length(camera_counts) != J || length(camera_days) != J) {
    stop(
      "run_USCR: camera_counts and camera_days must have length nrow(out).",
      call. = FALSE
    )
  }
  if (J < 2L) {
    stop(
      "run_USCR: at least 2 cameras are required. USCR estimates the spatial ",
      "scale parameter sigma from detection contrast across cameras at ",
      "different locations, which a single camera cannot provide.",
      call. = FALSE
    )
  }

  iter <- as.integer(iter)
  burnin <- as.integer(burnin)
  thin <- as.integer(max(1, thin))
  n_chains <- as.integer(max(1, n_chains))
  M <- as.integer(max(1, M))
  if (!is.null(iter_cap)) {
    iter_cap <- as.integer(max(1, iter_cap))
    iter <- min(iter, iter_cap)
  }
  if (!is.null(M_cap)) {
    M_cap <- as.integer(max(1, M_cap))
    M <- min(M, M_cap)
  }

  base_monitors <- c(
    "log_sigma", "log_lam_0", "psi", "sd_eps",
    "sigma", "lam_0", "N", "D_mi2",
    "sum_obs", "sum_sim", "bp"
  )
  latent_monitors <- c("z", "hrc", "eps")
  final_monitors <- if (isTRUE(diagnostic_mode)) {
    c(base_monitors, latent_monitors)
  } else {
    base_monitors
  }

  code <- build_uscr_code()

  current_M <- M
  current_iter <- iter
  current_burnin <- burnin
  current_thin <- thin
  current_buffer_m <- if (!is.null(buffer_m) && is.finite(buffer_m) && buffer_m > 0) {
    as.numeric(buffer_m)
  } else {
    uscr_default_buffer_m(log_sigma_mean, log_sigma_sd, p = buffer_p, multiplier = buffer_multiplier)
  }
  round_log <- list()
  current_rhat <- Inf
  round_i <- 0L
  prev_elapsed_min <- NULL
  prev_iter <- NULL
  fit <- NULL
  samples_list <- NULL
  samples_all <- NULL
  const <- NULL

  repeat {
    round_i <- round_i + 1L

    if (!isTRUE(adaptive) && round_i > 1L) {
      break
    }

    if (!is.null(max_adapt_rounds) && round_i > as.integer(max_adapt_rounds)) {
      warning(
        "USCR reached max_adapt_rounds = ",
        max_adapt_rounds,
        " before both the Rhat and M checks cleared. Proceeding with the latest fit."
      )
      break
    }

    const <- make_uscr_constants(
      out = out,
      camera_days = camera_days,
      M = current_M,
      log_sigma_mean = log_sigma_mean,
      log_sigma_sd = log_sigma_sd,
      log_lam0_mean = log_lam0_mean,
      log_lam0_sd = log_lam0_sd,
      sd_eps_shape = sd_eps_shape,
      sd_eps_rate = sd_eps_rate,
      buffer_m = current_buffer_m
    )

    data_list <- list(
      y = as.numeric(camera_counts),
      in_ss = rep(1L, const$M)
    )

    if (isTRUE(verbose)) {
      message(
        "USCR round ", round_i,
        ": M=", const$M,
        ", buffer=", round(current_buffer_m), "m",
        ", iter=", current_iter,
        ", burnin=", current_burnin,
        ", thin=", current_thin,
        ", chains=", n_chains
      )
    }
    report_status(
      round = round_i,
      iter = current_iter,
      burnin = current_burnin,
      M = const$M,
      buffer_m = current_buffer_m,
      prev_elapsed_min = prev_elapsed_min,
      prev_iter = prev_iter
    )
    round_started <- Sys.time()

    fit <- run_uscr_chains(
      code = code,
      constants = const,
      data = data_list,
      monitors = final_monitors,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      parallel_chains = parallel_chains,
      compute_WAIC = compute_WAIC,
      seed = if (is.null(seed)) NULL else seed + round_i
    )

    samples_list <- lapply(fit, `[[`, "samples")
    samples_all <- do.call(rbind, samples_list)
    psi_post <- samples_all[, "psi"]
    sigma_post_max <- max(samples_all[, "sigma"], na.rm = TRUE)

    current_rhat <- safe_rhat_max(samples_list)
    M_too_small <- mean(psi_post > psi_threshold, na.rm = TRUE) > psi_prob_cutoff
    converged <- is.na(current_rhat) || current_rhat <= rhat_target

    required_buffer_m <- uscr_buffer_from_sigma_km(sigma_post_max, p = buffer_p, multiplier = buffer_multiplier)
    buffer_too_small <- is.finite(required_buffer_m) && current_buffer_m < required_buffer_m

    round_log[[round_i]] <- data.frame(
      round = round_i,
      M = const$M,
      buffer_m = current_buffer_m,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      rhat_max = current_rhat,
      M_too_small = M_too_small,
      buffer_too_small = buffer_too_small,
      stringsAsFactors = FALSE
    )

    prev_elapsed_min <- as.numeric(difftime(Sys.time(), round_started, units = "mins"))
    prev_iter <- current_iter

    if (converged && !M_too_small && !buffer_too_small) {
      break
    }

    if (!isTRUE(adaptive)) {
      break
    }

    if (!converged) {
      current_iter <- current_iter * 2L
      current_burnin <- current_burnin * 2L
      current_thin <- current_thin * 2L
      if (!is.null(iter_cap)) {
        current_iter <- min(current_iter, iter_cap)
      }
    }

    if (M_too_small) {
      current_M <- current_M * 2L
      if (!is.null(M_cap)) {
        current_M <- min(current_M, M_cap)
      }
    }

    if (buffer_too_small) {
      current_buffer_m <- required_buffer_m * buffer_safety_factor
    }
  }

  build_uscr_result(
    fit_objects = fit,
    samples_list = samples_list,
    samples_all = samples_all,
    round_history = if (length(round_log) > 0L) do.call(rbind, round_log) else NULL,
    final_rhat_max = current_rhat,
    settings_M = const$M,
    settings_iter = current_iter,
    settings_burnin = current_burnin,
    settings_thin = current_thin,
    settings_buffer_m = current_buffer_m,
    area_mi2 = const$area_mi2
  )
}

# -------------------------------------------------------------------
# Suggested fast call for app / simulation contexts
# -------------------------------------------------------------------
# uscr_fit <- run_USCR(
#   out = out,
#   camera_counts = camera_counts,
#   camera_days = camera_days,
#   M = 100,
#   iter = 6000,
#   burnin = 1000,
#   thin = 5,
#   n_chains = 1,
#   adaptive = TRUE,
#   compute_WAIC = FALSE,
#   diagnostic_mode = FALSE,
#   seed = 123,
#   verbose = TRUE
# )
