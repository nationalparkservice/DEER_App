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

# Simulate USCR camera counts on a square camera grid. This mirrors the NIMBLE
# model in build_uscr_code(), so the simulated data and the fitted model share
# the same state space and observation process:
#   - Cameras: the same centered square grid the REM/TTE teaching simulator
#     uses (build_teaching_sim_grid()).
#   - Activity centers: N_rect ~ Poisson(D * rectangle area) points drawn with
#     runif() over the camera rectangle +/- buffer, keeping only those within
#     `buffer` of at least one camera (the model's in_ss constraint). The kept
#     points are a Poisson process with density D over exactly the model's
#     state space (the union of circles of radius `buffer` around cameras).
#   - Counts: y_j ~ Poisson(days_j * sum_i lambda0 * exp(-d_ij^2 / (2 sigma^2))
#     * exp(eps_j)), with camera-level noise eps_j ~ Normal(0, sd_eps).
# All distances are in km internally; buffer_m and sigma_m are in metres.
simulate_camera_counts <- function(n_side   = 5,
                                   spacing_m = 300,
                                   days      = 21,
                                   D_per_km2 = 25,
                                   lambda0   = 0.20,
                                   home_range_km2 = 0.89,
                                   sigma_m   = NULL,
                                   buffer_m  = NULL,
                                   sd_eps    = 0,
                                   seed      = 1) {
  if (n_side < 2 || days < 1 || spacing_m <= 0) {
    stop("USCR simulator needs at least a 2 x 2 grid, positive spacing, and at least 1 day.",
         call. = FALSE)
  }
  if (D_per_km2 < 0 || lambda0 <= 0 || sd_eps < 0) {
    stop("USCR simulator density, detection rate, and heterogeneity inputs are out of range.",
         call. = FALSE)
  }

  set.seed(seed)

  if (is.null(sigma_m)) {
    sigma_m <- home_range_km2_to_sigma_m(home_range_km2)
  }
  sigma_km <- sigma_m / 1000
  if (is.null(buffer_m) || !is.finite(buffer_m) || buffer_m <= 0) {
    buffer_m <- uscr_buffer_from_sigma_km(sigma_km)
  }
  buffer_km <- buffer_m / 1000

  # Cameras (detection radius is filled in later by sim_model_inputs()).
  out <- build_teaching_sim_grid(
    n_side = n_side,
    spacing_m = spacing_m,
    detection_radius_m = NA_real_,
    days = days
  )
  cam <- cbind(out$utm_e, out$utm_n)
  J <- nrow(cam)

  # Activity centers: uniform over the rectangle, then keep those within the
  # buffer of at least one camera.
  xlim <- range(cam[, 1]) + c(-1, 1) * buffer_km
  ylim <- range(cam[, 2]) + c(-1, 1) * buffer_km
  N_rect <- stats::rpois(1, D_per_km2 * diff(xlim) * diff(ylim))
  ac <- cbind(
    stats::runif(N_rect, xlim[1], xlim[2]),
    stats::runif(N_rect, ylim[1], ylim[2])
  )
  dist2 <- outer(ac[, 1], cam[, 1], "-")^2 + outer(ac[, 2], cam[, 2], "-")^2
  dist2 <- matrix(dist2, nrow = N_rect, ncol = J)
  keep <- if (N_rect > 0) apply(dist2, 1, min) < buffer_km^2 else logical(0)
  ac <- ac[keep, , drop = FALSE]
  dist2 <- dist2[keep, , drop = FALSE]

  # Counts.
  camera_days <- rep(as.numeric(days), J)
  Lambda <- colSums(lambda0 * exp(-dist2 / (2 * sigma_km^2)))  # per day, per camera
  eps <- stats::rnorm(J, mean = 0, sd = sd_eps)
  mu <- camera_days * Lambda * exp(eps)
  camera_counts <- stats::rpois(J, mu)

  area_km2 <- uscr_buffered_area_km2(cam, buffer_km)

  list(
    out = out,
    camera_counts = camera_counts,
    camera_days = camera_days,
    activity_centers = data.frame(x = ac[, 1] * 1000, y = ac[, 2] * 1000),  # metres
    state_space_outline = uscr_buffered_outline(cam, buffer_km, scale = 1000),  # metres
    truth = list(
      D_per_km2 = D_per_km2,
      N = nrow(ac),
      area_km2 = area_km2,
      realized_D_per_km2 = nrow(ac) / area_km2,
      lambda0   = lambda0,
      home_range_km2 = home_range_km2,
      sigma_m   = sigma_m,
      buffer_m  = buffer_m,
      sd_eps    = sd_eps,
      days      = days,
      spacing_m = spacing_m,
      n_side    = n_side,
      n_cams    = J
    )
  )
}

# For simulations: create "NPS-like" inputs for models
sim_model_inputs <- function(sim, detection_radius_m) {
  out <- sim$out
  out$`Detection Distance` <- detection_radius_m

  list(
    out           = out,
    camera_counts = sim$camera_counts,
    camera_days   = sim$camera_days
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
                                     max_adapt_rounds = 6,
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
  hit_max_adapt_rounds <- FALSE
  if (!is.null(max_adapt_rounds)) {
    max_adapt_rounds <- as.integer(max(1, max_adapt_rounds))
  }
  last_fit <- NULL
  prev_elapsed_min <- NULL
  prev_iter <- NULL

  repeat {
    round_i <- round_i + 1L

    if (!isTRUE(adaptive) && round_i > 1L) {
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

    # Stop here (before doubling) so the reported settings match the last fit.
    if (!is.null(max_adapt_rounds) && round_i >= max_adapt_rounds) {
      hit_max_adapt_rounds <- TRUE
      warning(
        method,
        " stopped after max_adapt_rounds = ",
        max_adapt_rounds,
        " without converging (Rhat). Proceeding with the latest fit."
      )
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
    hit_max_adapt_rounds = hit_max_adapt_rounds,
    unmet_checks = if (hit_max_adapt_rounds) "Rhat" else character(0),
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
                    max_adapt_rounds = 6,
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
                "bp", "sum_obs", "sum_sim")

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
                    max_adapt_rounds = 6,
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
                "bp", "sum_obs", "sum_sim")

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

# Union of circles of radius buffer_km around each camera (cam: J x 2 matrix
# of planar coordinates in km). This is the USCR state space implied by the
# in_ss constraint in build_uscr_code(). Planar sf geometry, no CRS needed.
uscr_buffered_union <- function(cam, buffer_km) {
  quiet_require("sf")
  pts <- sf::st_as_sf(
    data.frame(x = cam[, 1], y = cam[, 2]),
    coords = c("x", "y")
  )
  sf::st_union(sf::st_buffer(pts, dist = buffer_km))
}

uscr_buffered_area_km2 <- function(cam, buffer_km) {
  as.numeric(sf::st_area(uscr_buffered_union(cam, buffer_km)))
}

# Outline of the buffered state space as a data frame (x, y, group) for
# plotting; `scale` converts km to the plot's units (1000 for metres).
uscr_buffered_outline <- function(cam, buffer_km, scale = 1) {
  xy <- sf::st_coordinates(uscr_buffered_union(cam, buffer_km))
  grp_cols <- intersect(c("L1", "L2", "L3"), colnames(xy))
  data.frame(
    x = xy[, "X"] * scale,
    y = xy[, "Y"] * scale,
    group = interaction(as.data.frame(xy[, grp_cols, drop = FALSE]), drop = TRUE)
  )
}

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

    area_km2 <- as.numeric(sf::st_area(cam_buff)) / 1e6  # m^2 -> km^2
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
    # Simulated grids (planar km coordinates, no Longitude/Latitude): same
    # union-of-circles state space as field data, computed in km directly.
    area_km2 <- uscr_buffered_area_km2(cbind(out$utm_e, out$utm_n), buffer)
  }

  list(
    buffer = buffer,
    buffer_sq = buffer_sq,
    xlim = xlim,
    ylim = ylim,
    area_km2 = area_km2
  )
}

# -------------------------------------------------------------------
# USCR: density prior via psi ~ Beta(a, b)
# -------------------------------------------------------------------
# With data augmentation, N = sum(z[1:M]) and z[i] | psi ~ Bernoulli(psi), so
# psi ~ Beta(a, b) gives N ~ Beta-Binomial(M, a, b), and density is
# D = N / area_km2. a = b = 1 (psi ~ Uniform(0, 1)) is the default: a uniform
# prior on N over 0..M.
#
# An informative prior is specified as a 95% interval for density
# (animals/km^2). Beta(a, b) has two parameters, so both bounds can usually be
# matched exactly; the prior's mean and median follow from the interval and
# are reported back. Because M and the
# state-space area can change between adaptive rounds, run_USCR() re-solves
# the prior every round so the prior on density stays the same.
#
# Adapted from dbbinom(), bb_summary() and design_psi_prior() in the
# Poisson_Lognormal_USCR_Simulation project, with these changes: the target is
# the 95% interval only (the original pinned the mean and matched the interval
# approximately); quantiles are interpolated between integer N values, so the
# search is smooth (exact integer quantiles are piecewise constant); and misses
# are measured as fractions of the interval's width.

# Beta-Binomial(M, a, b) probability mass function
dbbinom <- function(k, M, a, b) {
  exp(lchoose(M, k) + lbeta(k + a, M - k + b) - lbeta(a, b))
}

# Mean, median, SD and (interpolated) 95% interval of Beta-Binomial(M, a, b)
bb_summary <- function(M, a, b) {
  k <- 0:M
  p <- dbbinom(k, M, a, b)
  p <- p / sum(p)
  cdf <- cumsum(p)
  q <- function(prob) {
    stats::approx(cdf, k, xout = prob, ties = min, rule = 2)$y
  }
  mean_N <- sum(k * p)
  list(
    mean   = mean_N,
    median = q(0.5),
    sd     = sqrt(sum((k - mean_N)^2 * p)),
    ci     = c(lo = q(0.025), hi = q(0.975))
  )
}

# Check user-supplied density prior bounds. Returns NULL if valid, otherwise a
# message suitable for showing to app users.
uscr_check_density_prior <- function(lo, hi) {
  vals <- suppressWarnings(as.numeric(c(lo, hi)))
  if (length(vals) != 2L || any(!is.finite(vals))) {
    return("Enter a lower 2.5% bound and an upper 97.5% bound for density (animals/km²).")
  }
  if (vals[1] < 0) {
    return("Density prior bounds can't be negative.")
  }
  if (!(vals[1] < vals[2])) {
    return("The lower 2.5% bound must be less than the upper 97.5% bound.")
  }
  NULL
}

# Starting M needed for an informative density prior: about 2x the prior's
# upper (97.5%) abundance, rounded up to a multiple of 50.
uscr_prior_start_M <- function(D_hi, area_km2) {
  as.integer(ceiling(2 * D_hi * area_km2 / 50) * 50)
}

# Uniform prior (psi ~ Beta(1, 1)) in the same format as the informative one.
uscr_uniform_density_prior <- function(M, area_km2) {
  list(uniform = TRUE, a = 1, b = 1, M = as.integer(M), area_km2 = area_km2)
}

# Solve for psi ~ Beta(a, b) so that D = N / area_km2 has the requested 95%
# prior interval, for a given M. With two Beta parameters and two targets, both
# bounds can usually be matched exactly; the prior's mean and median then
# follow from the interval (for strongly skewed intervals they may not be
# where a user expects, so they are reported back).
design_uscr_density_prior <- function(D_lo, D_hi, area_km2, M,
                                      miss_tolerance = 0.10) {
  msg <- uscr_check_density_prior(D_lo, D_hi)
  if (!is.null(msg)) stop(msg, call. = FALSE)
  M <- as.integer(M)

  N_lo <- D_lo * area_km2
  N_hi <- D_hi * area_km2

  # N can never exceed M, so the upper bound must sit below M.
  if (N_hi >= M) {
    stop(sprintf(paste0(
      "The density prior's upper bound (%.3g animals/km² over a %.2f km² ",
      "state space, about %.0f animals) needs M above %d. Increase M or lower ",
      "the upper bound."), D_hi, area_km2, N_hi, M), call. = FALSE)
  }

  # Misses are measured as fractions of the interval's width, matching the
  # off-target check below.
  width <- N_hi - N_lo
  miss_fn <- function(s) {
    ((s$ci[["lo"]] - N_lo) / width)^2 + ((s$ci[["hi"]] - N_hi) / width)^2
  }

  # Parameterize by the Beta mean (logit scale) and concentration a + b (log
  # scale). Start from the interval midpoint with the best concentration for
  # that mean, then search both together.
  p_start <- min(max(((N_lo + N_hi) / 2) / M, 1e-6), 1 - 1e-6)
  k_start <- exp(stats::optimize(function(lk) {
    miss_fn(bb_summary(M, p_start * exp(lk), (1 - p_start) * exp(lk)))
  }, interval = log(c(0.01, 1e6)))$minimum)
  obj <- function(par) {
    p <- stats::plogis(par[1])
    k <- exp(par[2])
    miss_fn(bb_summary(M, p * k, (1 - p) * k))
  }
  opt <- stats::optim(c(stats::qlogis(p_start), log(k_start)), obj,
                      method = "Nelder-Mead",
                      control = list(maxit = 2000, reltol = 1e-12))
  p <- stats::plogis(opt$par[1])
  k <- exp(opt$par[2])
  a <- p * k
  b <- (1 - p) * k

  s <- bb_summary(M, a, b)
  achieved <- c(
    lo     = s$ci[["lo"]] / area_km2,
    hi     = s$ci[["hi"]] / area_km2,
    mean   = s$mean / area_km2,
    median = s$median / area_km2
  )
  # How far each bound moved, as a fraction of the requested interval's width
  # (so asking for 24-26 and getting 22-28 counts as a large miss).
  miss <- c(
    lo = abs(achieved[["lo"]] - D_lo) / (D_hi - D_lo),
    hi = abs(achieved[["hi"]] - D_hi) / (D_hi - D_lo)
  )

  list(
    uniform = FALSE,
    a = a,
    b = b,
    M = M,
    area_km2 = area_km2,
    target = c(lo = D_lo, hi = D_hi),
    achieved = achieved,
    miss = miss,
    off_target = any(miss > miss_tolerance),
    miss_tolerance = miss_tolerance
  )
}

# One-line description of a density prior, for run-status panels.
describe_uscr_density_prior <- function(design) {
  if (is.null(design)) return(NULL)
  if (isTRUE(design$uniform)) {
    return(sprintf("Uniform (psi ~ Beta(1, 1); N uniform on 0 to M = %d)", design$M))
  }
  sprintf(paste0(
    "Informative. Target 95%% interval %.3g to %.3g animals/km²; achieved ",
    "%.3g to %.3g (mean %.3g, median %.3g) with psi ~ Beta(%.3g, %.3g), M = %d, ",
    "state-space area %.2f km²"),
    design$target[["lo"]], design$target[["hi"]],
    design$achieved[["lo"]], design$achieved[["hi"]],
    design$achieved[["mean"]], design$achieved[["median"]],
    design$a, design$b, design$M, design$area_km2)
}

# Plot the prior on density implied by a design from design_uscr_density_prior()
# (exact Beta-Binomial probabilities, rescaled from N to animals/km^2), with the
# target 95% bounds (dashed) and the prior's mean (solid).
plot_D_prior_predictive <- function(design, caption = NULL) {
  quiet_require("ggplot2")
  k <- 0:design$M
  p <- dbbinom(k, design$M, design$a, design$b)
  p <- p / sum(p)
  cdf <- cumsum(p)
  # Trim the far tails so the plot shows the bulk of the prior.
  keep <- cdf >= 0.0005 & c(0, cdf[-length(cdf)]) <= 0.9995
  df <- data.frame(D = k[keep] / design$area_km2, prob = p[keep])

  fig <- ggplot2::ggplot(df, ggplot2::aes(x = D, y = prob)) +
    ggplot2::geom_col(width = 1.02 / design$area_km2, fill = "#609048", alpha = 0.75,
                      position = "identity")  # slight overlap avoids hairline gaps +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Density (animals/km²)",
      y = "Prior probability",
      title = "Prior on density implied by your inputs",
      subtitle = "Dashed lines: your 95% interval. Solid line: the prior's mean.",
      caption = caption
    )

  if (!is.null(design$target)) {
    fig <- fig +
      ggplot2::geom_vline(xintercept = design$achieved[["mean"]], color = "#303018") +
      ggplot2::geom_vline(xintercept = design$target[c("lo", "hi")],
                          color = "firebrick", linetype = "dashed")
  }
  fig
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
    psi ~ dbeta(psi_a, psi_b)  # prior on N (and density); Beta(1, 1) = uniform
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
    D <- N / area_km2  # animals per km^2
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
                                buffer_m,
                                psi_a = 1,
                                psi_b = 1) {
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
    area_km2 = ss$area_km2,
    cam = as.matrix(cbind(out$utm_e, out$utm_n)),
    log_sigma_mean = log_sigma_mean,
    log_sigma_sd = log_sigma_sd,
    log_lam0_mean = log_lam0_mean,
    log_lam0_sd = log_lam0_sd,
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate = sd_eps_rate,
    psi_a = psi_a,
    psi_b = psi_b
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
# D_km2. All three models monitor density as "D", in animals per km^2.
extract_param_draws <- function(samples_all, param_name) {
  cn <- colnames(samples_all)
  if (identical(param_name, "D_km2")) {
    if ("D" %in% cn) return(as.numeric(samples_all[, "D"]))
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
      Median       = stats::median(draws, na.rm = TRUE),
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
                     density_prior = NULL,
                     M_check_fraction = 0.9,
                     M_check_prob = 0.01,
                     buffer_p = 0.99,
                     buffer_multiplier = 1.05,
                     buffer_safety_factor = 1.01,
                     max_adapt_rounds = 6,
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
                                area_km2 = NA_real_) {
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
        n_chains = n_chains,
        compute_WAIC = compute_WAIC,
        diagnostic_mode = diagnostic_mode,
        adaptive = adaptive,
        max_adapt_rounds = max_adapt_rounds,
        rhat_target = rhat_target,
        M_check_fraction = M_check_fraction,
        M_check_prob = M_check_prob,
        buffer_m = settings_buffer_m,
        buffer_p = buffer_p,
        buffer_multiplier = buffer_multiplier,
        area_km2 = area_km2,
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

  # density_prior: NULL or list(uniform = TRUE) for the default uniform prior,
  # or list(uniform = FALSE, lo = , hi = ): a 95% interval in animals/km^2.
  informative_prior <- !is.null(density_prior) && identical(density_prior$uniform, FALSE)
  if (informative_prior) {
    prior_msg <- uscr_check_density_prior(density_prior$lo, density_prior$hi)
    if (!is.null(prior_msg)) stop(prior_msg, call. = FALSE)
  }
  prior_design <- NULL

  base_monitors <- c(
    "log_sigma", "log_lam_0", "psi", "sd_eps",
    "sigma", "lam_0", "N", "D",
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
  hit_max_adapt_rounds <- FALSE
  unmet_checks <- character(0)
  if (!is.null(max_adapt_rounds)) {
    max_adapt_rounds <- as.integer(max(1, max_adapt_rounds))
  }
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

    # Re-solve the density prior for this round's state-space area and M, so
    # the prior on density stays fixed as the buffer and M change.
    round_area_km2 <- uscr_state_space_and_area(
      out,
      buffer_m = current_buffer_m,
      log_sigma_mean = log_sigma_mean,
      log_sigma_sd = log_sigma_sd
    )$area_km2
    if (informative_prior) {
      current_M <- max(current_M, uscr_prior_start_M(density_prior$hi, round_area_km2))
      prior_design <- design_uscr_density_prior(
        D_lo = density_prior$lo,
        D_hi = density_prior$hi,
        area_km2 = round_area_km2,
        M = current_M
      )
    } else {
      prior_design <- uscr_uniform_density_prior(current_M, round_area_km2)
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
      buffer_m = current_buffer_m,
      psi_a = prior_design$a,
      psi_b = prior_design$b
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
    N_post <- samples_all[, "N"]
    sigma_post_max <- max(samples_all[, "sigma"], na.rm = TRUE)

    current_rhat <- safe_rhat_max(samples_list)
    # M is too small if the posterior for N presses against it.
    M_too_small <- mean(N_post > M_check_fraction * const$M, na.rm = TRUE) > M_check_prob
    converged <- is.na(current_rhat) || current_rhat <= rhat_target

    required_buffer_m <- uscr_buffer_from_sigma_km(sigma_post_max, p = buffer_p, multiplier = buffer_multiplier)
    buffer_too_small <- is.finite(required_buffer_m) && current_buffer_m < required_buffer_m

    round_log[[round_i]] <- data.frame(
      round = round_i,
      M = const$M,
      buffer_m = current_buffer_m,
      psi_a = prior_design$a,
      psi_b = prior_design$b,
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

    # Stop here (before growing iter/M/buffer) so the reported settings match
    # the last fit.
    if (!is.null(max_adapt_rounds) && round_i >= max_adapt_rounds) {
      hit_max_adapt_rounds <- TRUE
      unmet_checks <- c(
        if (!converged) "Rhat",
        if (M_too_small) "M",
        if (buffer_too_small) "state-space buffer"
      )
      warning(
        "USCR stopped after max_adapt_rounds = ",
        max_adapt_rounds,
        " without passing these checks: ",
        paste(unmet_checks, collapse = ", "),
        ". Proceeding with the latest fit."
      )
      break
    }

    if (!converged) {
      current_iter <- current_iter * 2L
      current_burnin <- current_burnin * 2L
      current_thin <- current_thin * 2L
    }

    if (M_too_small) {
      current_M <- current_M * 2L
    }

    if (buffer_too_small) {
      current_buffer_m <- required_buffer_m * buffer_safety_factor
    }
  }

  result <- build_uscr_result(
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
    area_km2 = const$area_km2
  )
  result$hit_max_adapt_rounds <- hit_max_adapt_rounds
  result$unmet_checks <- unmet_checks
  result$density_prior <- prior_design  # the prior used in the final round
  result
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
