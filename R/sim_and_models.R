## R/sim_and_models.R
## Simulation helpers, model inputs for uploaded data, and the USCR, REM,
## and TTE models (which follow Camera_Unmarked_Analyses_updated.Rmd, with
## adaptive reruns for app use).

# -------------------------------------------------------------------
# Package checks (loaded lazily where needed)
# -------------------------------------------------------------------

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.", call. = FALSE)
  }
}

# -------------------------------------------------------------------
# 1. Simulation helpers (Simulate data tab)
# -------------------------------------------------------------------

# Average length of a straight-line path across a camera viewshed (a circular
# sector with full angle theta_deg degrees), as a multiple of the viewshed's
# radius. The quadratic was fitted to simulations of random straight-line paths
# across circular sectors with central angles from 5 to 90 degrees.
tte_path_length_multiplier <- function(theta_deg) {
  theta_deg <- as.numeric(theta_deg)
  3.324e-01 + 5.580e-03 * theta_deg - 1.454e-05 * theta_deg^2
}

# Square grid of n_side x n_side cameras spacing_m apart, centered on the
# origin: x and y in meters (for plots and tables), utm_e and utm_n in km
# (for the models).
build_sim_camera_grid <- function(n_side, spacing_m) {
  coords <- expand.grid(
    x = seq(0, by = spacing_m, length.out = n_side),
    y = seq(0, by = spacing_m, length.out = n_side)
  )
  coords <- coords[order(coords$y, coords$x), , drop = FALSE]

  coords |>
    dplyr::mutate(
      x = x - mean(range(x)),
      y = y - mean(range(y)),
      `Site Name` = paste0("C", dplyr::row_number()),
      utm_e = x / 1000,
      utm_n = y / 1000
    )
}

# Simulate USCR camera counts on a square camera grid. This mirrors the NIMBLE
# model in build_uscr_code(), so the simulated data and the fitted model share
# the same state space and observation process:
#   - Cameras: the grid from build_sim_camera_grid().
#   - Activity centers: N_rect ~ Poisson(D * rectangle area) points drawn with
#     runif() over the camera rectangle +/- buffer, keeping only those within
#     the buffer of at least one camera (the model's in_ss constraint). The
#     kept points are a Poisson process with density D over exactly the
#     model's state space (the union of circles of radius buffer around the
#     cameras).
#   - Counts: y_j ~ Poisson(days * sum_i lambda0 * exp(-d_ij^2 / (2 sigma^2))
#     * exp(use_eps * eps_j)), with camera-level noise eps_j ~ Normal(0, sd_eps)
#     and use_eps = 1 when camera_noise is on, 0 when off. eps is always
#     drawn, so the same seed gives the same activity centers and counts
#     (apart from the noise) with noise on or off, mirroring the model.
# sigma_km is in km and buffer_m in meters; the activity centers and the
# state-space outline are returned in meters, for plotting.
simulate_camera_counts <- function(n_side,
                                   spacing_m,
                                   days,
                                   D_per_km2,
                                   lambda0,
                                   sigma_km,
                                   buffer_m,
                                   sd_eps,
                                   camera_noise,
                                   seed) {
  set.seed(seed)

  out <- build_sim_camera_grid(n_side, spacing_m)
  cam <- cbind(out$utm_e, out$utm_n)
  J <- nrow(cam)
  buffer_km <- buffer_m / 1000

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
  keep <- if (N_rect > 0) apply(dist2, 1, min) < buffer_km^2 else logical(0)
  ac <- ac[keep, , drop = FALSE]
  dist2 <- dist2[keep, , drop = FALSE]

  # Counts.
  camera_days <- rep(as.numeric(days), J)
  Lambda <- colSums(lambda0 * exp(-dist2 / (2 * sigma_km^2)))  # per camera-day
  eps <- stats::rnorm(J, mean = 0, sd = sd_eps)
  use_eps <- as.numeric(isTRUE(camera_noise))
  mu <- camera_days * Lambda * exp(use_eps * eps)

  list(
    out = out,
    camera_counts = stats::rpois(J, mu),
    camera_days = camera_days,
    activity_centers = data.frame(x = ac[, 1] * 1000, y = ac[, 2] * 1000),
    state_space_outline = uscr_buffered_outline(cam, buffer_km, scale = 1000),
    N = nrow(ac),
    area_km2 = uscr_buffered_area_km2(cam, buffer_km)
  )
}

# Simulate REM or TTE camera counts on a square camera grid, from the
# expected count per camera that the model uses (see run_REM() and
# run_TTE()):
#   - REM: D * days * v * r * (2 + theta) / pi, with theta in radians.
#   - TTE: D * days * U * A, where U = v / (f(theta) * r) is the number of
#     viewshed-crossing times per day, f() is tte_path_length_multiplier(),
#     and A = pi * r^2 * theta / 360 is the viewshed area, with theta in
#     degrees.
# Each camera's expected count is multiplied by exp(use_eps * eps_j), with
# camera-level noise eps_j ~ Normal(0, sd_eps) and use_eps = 1 when
# camera_noise is on, 0 when off. eps is always drawn, so the same seed gives
# the same counts (apart from the noise) with noise on or off, mirroring the
# model.
simulate_encounter_counts <- function(model,
                                      n_side,
                                      spacing_m,
                                      days,
                                      D_per_km2,
                                      detection_radius_m,
                                      theta_deg,
                                      v_km_day,
                                      sd_eps,
                                      camera_noise,
                                      seed) {
  set.seed(seed)

  out <- build_sim_camera_grid(n_side, spacing_m) |>
    dplyr::mutate(`Detection Distance` = detection_radius_m)
  J <- nrow(out)
  camera_days <- rep(as.numeric(days), J)
  r_km <- detection_radius_m / 1000
  eps <- stats::rnorm(J, mean = 0, sd = sd_eps)
  use_eps <- as.numeric(isTRUE(camera_noise))

  per_day <- if (identical(model, "REM")) {
    D_per_km2 * v_km_day * r_km * (2 + theta_deg * pi / 180) / pi
  } else {
    crossings_per_day <- v_km_day / (tte_path_length_multiplier(theta_deg) * r_km)
    viewshed_km2 <- pi * r_km^2 * theta_deg / 360
    D_per_km2 * crossings_per_day * viewshed_km2
  }
  mu <- camera_days * per_day * exp(use_eps * eps)

  list(
    out = out,
    camera_counts = stats::rpois(J, mu),
    camera_days = camera_days
  )
}

# -------------------------------------------------------------------
# 2. Uploaded data: per-camera counts and effort for the models
# -------------------------------------------------------------------
build_nps_model_inputs <- function(
    deployments,
    images,
    species_to_analyze,
    max_days = NULL
    ) {
  quiet_require("dplyr")
  quiet_require("sf")

  # Each camera's analysis window (start/end date-times, optional trim) comes
  # from analysis_windows() in R/data_checks.R; check_deployments() has
  # already blocked any window with no time in it.
  windows <- analysis_windows(deployments, max_days = max_days)
  deps <- deployments |>
    dplyr::mutate(
      window_start = windows$window_start,
      window_end   = windows$window_end,
      window_days  = windows$days
    )

  # check_images() has already parsed every Timestamp to POSIXct.
  # Latitude/Longitude are validated in check_deployments(), but may still be
  # character columns from the CSV import.
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

  # Detection events (see detection_events() in R/data_checks.R): an event
  # counts if its first image is inside the camera's analysis window, and the
  # part of its duration inside the window is removed from the camera's effort.
  events <- detection_events(images, deployments, max_days = max_days) |>
    dplyr::filter(Species == species_to_analyze, in_window) |>
    dplyr::group_by(`Site Name`) |>
    dplyr::summarize(
      events_count = sum(count),
      event_days   = sum(days_in_window),
      .groups      = "drop"
    )

  out <- deps |>
    dplyr::left_join(events, by = "Site Name") |>
    dplyr::mutate(
      events_count = dplyr::coalesce(events_count, 0),
      event_days   = dplyr::coalesce(event_days, 0)
    )

  list(
    out           = dplyr::select(out, -events_count, -event_days),
    camera_counts = out$events_count,
    camera_days   = out$window_days - out$event_days
  )
}

# -------------------------------------------------------------------
# 3. MCMC chain runner (all models)
# -------------------------------------------------------------------

# nimbleMCMC() returns just the samples matrix when WAIC is off; return the
# same list(samples, WAIC) shape either way.
normalize_nimble_output <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(list(samples = as.matrix(x), WAIC = NULL))
  }
  x
}

# Runs the chains in parallel worker processes (one per chain, up to the
# number of cores minus one). make_inits(constants), if given, returns each
# chain's initial values (USCR uses it to start activity centers near
# cameras); otherwise nimble chooses initial values itself.
run_nimble_chains <- function(code,
                              constants,
                              data,
                              monitors,
                              niter,
                              nburnin,
                              thin,
                              n_chains,
                              compute_WAIC,
                              make_inits = NULL) {
  quiet_require("parallel")

  cl_size <- min(n_chains, max(1L, parallel::detectCores() - 1L))
  cl <- parallel::makeCluster(cl_size)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::parLapply(
    cl = cl,
    X = seq_len(n_chains),
    fun = function(chain_id, code, constants, data, monitors,
                   niter, nburnin, thin, compute_WAIC,
                   make_inits, normalize) {
      library(nimble)

      mcmc_args <- list(
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
      # nimbleMCMC() doesn't accept inits = NULL, so add inits only when
      # there are some. quote = TRUE passes the model code as is, rather than
      # evaluating it.
      if (!is.null(make_inits)) {
        mcmc_args$inits <- make_inits(constants)
      }
      normalize(do.call(nimble::nimbleMCMC, mcmc_args, quote = TRUE))
    },
    code = code,
    constants = constants,
    data = data,
    monitors = monitors,
    niter = niter,
    nburnin = nburnin,
    thin = thin,
    compute_WAIC = compute_WAIC,
    make_inits = make_inits,
    normalize = normalize_nimble_output
  )
}

# -------------------------------------------------------------------
# 4. REM and TTE: adaptive fit
# -------------------------------------------------------------------

# Fits REM or TTE (model code and constants from run_REM() or run_TTE()), then
# reruns it while any check fails, up to max_adapt_rounds rounds: iterations,
# burn-in, and thinning double if Rhat is above rhat_target, and D_max doubles
# if the posterior for D presses against it. The checks' thresholds keep their
# defaults.
run_encounter_rate_model <- function(method,
                                     code,
                                     constants,
                                     data,
                                     monitors,
                                     iter,
                                     burnin,
                                     thin,
                                     n_chains,
                                     D_max,
                                     density_prior,
                                     max_adapt_rounds,
                                     compute_WAIC,
                                     status_callback,
                                     rhat_target = 1.1,
                                     D_check_fraction = 0.9,
                                     D_check_prob = 0.01) {
  iter <- as.integer(iter)
  burnin <- as.integer(burnin)
  thin <- as.integer(max(1, thin))
  n_chains <- as.integer(max(1, n_chains))
  max_adapt_rounds <- as.integer(max(1, max_adapt_rounds))

  # density_prior: list(uniform = TRUE) for D ~ Uniform(0, D_max), or
  # list(uniform = FALSE, lo = , hi = ): a 95% interval in animals/km^2 (the
  # app checks the interval before running). D_max is the starting upper
  # limit under the uniform prior; an informative prior starts it at twice
  # the interval's upper bound instead.
  informative_prior <- identical(density_prior$uniform, FALSE)
  current_D_max <- if (informative_prior) {
    encounter_prior_start_D_max(density_prior$hi)
  } else {
    as.numeric(D_max)
  }

  current_iter <- iter
  current_burnin <- burnin
  current_thin <- thin
  round_i <- 0L
  current_rhat <- Inf
  round_log <- list()
  hit_max_adapt_rounds <- FALSE
  unmet_checks <- character(0)
  prior_design <- NULL
  prev_elapsed_min <- NULL
  prev_iter <- NULL
  fit <- NULL
  samples_list <- NULL
  samples_all <- NULL

  repeat {
    round_i <- round_i + 1L

    status_callback(
      round = round_i,
      iter = current_iter,
      burnin = current_burnin,
      D_max = current_D_max,
      prev_elapsed_min = prev_elapsed_min,
      prev_iter = prev_iter
    )
    round_started <- Sys.time()

    # Re-solve the density prior for this round's D_max, so the requested
    # interval holds as D_max grows.
    prior_design <- if (informative_prior) {
      design_beta_density_prior(density_prior$lo, density_prior$hi, current_D_max)
    } else {
      uniform_beta_density_prior(current_D_max)
    }
    constants$D_max <- current_D_max
    constants$D_a <- prior_design$a
    constants$D_b <- prior_design$b

    fit <- run_nimble_chains(
      code = code,
      constants = constants,
      data = data,
      monitors = monitors,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      compute_WAIC = compute_WAIC
    )

    samples_list <- lapply(fit, `[[`, "samples")
    samples_all <- do.call(rbind, samples_list)
    current_rhat <- safe_rhat_max(samples_list)
    converged <- is.na(current_rhat) || current_rhat <= rhat_target

    # D_max is too small if the posterior for D presses against it.
    D_post <- samples_all[, "D"]
    D_too_small <- mean(D_post > D_check_fraction * current_D_max, na.rm = TRUE) > D_check_prob

    round_log[[round_i]] <- data.frame(
      round = round_i,
      D_max = current_D_max,
      D_a = prior_design$a,
      D_b = prior_design$b,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      rhat_max = current_rhat,
      D_too_small = D_too_small,
      stringsAsFactors = FALSE
    )

    prev_elapsed_min <- as.numeric(difftime(Sys.time(), round_started, units = "mins"))
    prev_iter <- current_iter

    if (converged && !D_too_small) {
      break
    }

    # Stop here (before doubling) so the reported settings match the last fit.
    if (round_i >= max_adapt_rounds) {
      hit_max_adapt_rounds <- TRUE
      unmet_checks <- c(if (!converged) "Rhat", if (D_too_small) "D_max")
      warning(
        method,
        " stopped after max_adapt_rounds = ",
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
    if (D_too_small) {
      current_D_max <- current_D_max * 2
    }
  }

  list(
    method = method,
    samples_list = samples_list,
    samples_all = samples_all,
    waic = extract_waic_mean(fit),
    fit_objects = fit,
    round_history = do.call(rbind, round_log),
    final_rhat_max = current_rhat,
    hit_max_adapt_rounds = hit_max_adapt_rounds,
    unmet_checks = unmet_checks,
    density_prior = prior_design,  # the prior used in the final round
    settings = list(
      iter = current_iter,
      burnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      max_adapt_rounds = max_adapt_rounds,
      rhat_target = rhat_target,
      D_max = current_D_max,
      D_check_fraction = D_check_fraction,
      D_check_prob = D_check_prob
    )
  )
}

# -------------------------------------------------------------------
# 5. REM model
# -------------------------------------------------------------------
# Random encounter model (Rowcliffe et al. 2008): the expected count at camera
# j is D * days_j * v * r_j * (2 + theta_j) / pi, with theta_j in radians.

# The app supplies every setting without a default here (its constants live
# in app.R) and checks the inputs before running. theta_deg is one detection
# angle for every camera (simulated data) or one per camera (uploaded data).
run_REM <- function(y,
                    r_km,
                    camera_days,
                    theta_deg,
                    iter,
                    burnin,
                    thin,
                    n_chains,
                    D_max,
                    log_v_mean,
                    log_v_sd,
                    sd_eps_shape,
                    sd_eps_rate,
                    density_prior,
                    camera_noise,
                    max_adapt_rounds,
                    compute_WAIC,
                    status_callback) {
  quiet_require("nimble")

  J <- length(y)
  theta_deg <- rep_len(as.numeric(theta_deg), J)

  code <- nimble::nimbleCode({
    p_D ~ dbeta(D_a, D_b)  # density prior: D = D_max * p_D; Beta(1, 1) = Uniform(0, D_max)
    D <- D_max * p_D
    log_v ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dgamma(sd_eps_shape, sd_eps_rate)

    v <- exp(log_v)

    for (j in 1:J) {
      eps[j] ~ dnorm(0, sd = sd_eps)

      # Camera-level noise (Poisson-lognormal) only when use_eps = 1. With
      # use_eps = 0, eps and sd_eps are sampled from their priors but don't
      # affect the data, and they are not monitored.
      log(lambda[j]) <- log(D) + log(camera_days[j]) + log(v) + log(r[j]) +
        log(2 + theta[j]) - log(pi_const) + use_eps * eps[j]

      y[j] ~ dpois(lambda[j])
      y_sim[j] ~ dpois(lambda[j])

      # Posterior predictive check (see ppc_summary()).
      pearson_obs[j] <- (y[j] - lambda[j])^2 / lambda[j]
      pearson_sim[j] <- (y_sim[j] - lambda[j])^2 / lambda[j]
    }

    sum_obs <- sum(pearson_obs[1:J])
    sum_sim <- sum(pearson_sim[1:J])
    bp <- step(sum_sim - sum_obs)
  })

  const <- list(
    J = J,
    r = as.numeric(r_km),
    camera_days = as.numeric(camera_days),
    pi_const = pi,
    theta = theta_deg * pi / 180,
    log_v_mean = log_v_mean,
    log_v_sd = log_v_sd,
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate = sd_eps_rate,
    use_eps = as.numeric(isTRUE(camera_noise))
  )

  data_list <- list(y = as.numeric(y))

  # sd_eps is only meaningful (and only monitored) with camera-level noise on.
  monitors <- c("D", "v", if (isTRUE(camera_noise)) "sd_eps",
                "bp", "sum_obs", "sum_sim")

  fit <- run_encounter_rate_model(
    method = "REM",
    code = code,
    constants = const,
    data = data_list,
    monitors = monitors,
    iter = iter,
    burnin = burnin,
    thin = thin,
    n_chains = n_chains,
    D_max = D_max,
    density_prior = density_prior,
    max_adapt_rounds = max_adapt_rounds,
    compute_WAIC = compute_WAIC,
    status_callback = status_callback
  )
  fit$settings$camera_noise <- isTRUE(camera_noise)
  fit
}

# -------------------------------------------------------------------
# 6. TTE model
# -------------------------------------------------------------------
# Time-to-event model (Moeller et al. 2018), fit as the equivalent Poisson
# count model: the expected count at camera j is D * days_j * U_j * A_j, where
# U_j = v / (f(theta_j) * r_j) is the number of viewshed-crossing times per
# day, f() is tte_path_length_multiplier(), A_j = pi * r_j^2 * theta_j / 360
# is the viewshed area, and theta_j is in degrees.

# The app supplies every setting without a default here (its constants live
# in app.R) and checks the inputs before running. theta_deg is one detection
# angle for every camera (simulated data) or one per camera (uploaded data).
run_TTE <- function(y,
                    r_km,
                    camera_days,
                    theta_deg,
                    iter,
                    burnin,
                    thin,
                    n_chains,
                    D_max,
                    log_v_mean,
                    log_v_sd,
                    sd_eps_shape,
                    sd_eps_rate,
                    density_prior,
                    camera_noise,
                    max_adapt_rounds,
                    compute_WAIC,
                    status_callback) {
  quiet_require("nimble")

  J <- length(y)
  r_km <- as.numeric(r_km)
  theta_deg <- rep_len(as.numeric(theta_deg), J)

  code <- nimble::nimbleCode({
    p_D ~ dbeta(D_a, D_b)  # density prior: D = D_max * p_D; Beta(1, 1) = Uniform(0, D_max)
    D <- D_max * p_D
    log_v ~ dnorm(log_v_mean, sd = log_v_sd)
    sd_eps ~ dgamma(sd_eps_shape, sd_eps_rate)

    v <- exp(log_v)

    for (j in 1:J) {
      U[j] <- v / path_length[j]
      eps[j] ~ dnorm(0, sd = sd_eps)

      # Camera-level noise (Poisson-lognormal) only when use_eps = 1. With
      # use_eps = 0, eps and sd_eps are sampled from their priors but don't
      # affect the data, and they are not monitored.
      log(lambda[j]) <- log(D) + log(camera_days[j]) + log(U[j]) + log(A[j]) +
        use_eps * eps[j]

      y[j] ~ dpois(lambda[j])
      y_sim[j] ~ dpois(lambda[j])

      # Posterior predictive check (see ppc_summary()).
      pearson_obs[j] <- (y[j] - lambda[j])^2 / lambda[j]
      pearson_sim[j] <- (y_sim[j] - lambda[j])^2 / lambda[j]
    }

    sum_obs <- sum(pearson_obs[1:J])
    sum_sim <- sum(pearson_sim[1:J])
    bp <- step(sum_sim - sum_obs)
  })

  const <- list(
    J = J,
    camera_days = as.numeric(camera_days),
    path_length = tte_path_length_multiplier(theta_deg) * r_km,
    A = pi * r_km^2 * theta_deg / 360,
    log_v_mean = log_v_mean,
    log_v_sd = log_v_sd,
    sd_eps_shape = sd_eps_shape,
    sd_eps_rate = sd_eps_rate,
    use_eps = as.numeric(isTRUE(camera_noise))
  )

  data_list <- list(y = as.numeric(y))

  # sd_eps is only meaningful (and only monitored) with camera-level noise on.
  monitors <- c("D", "v", if (isTRUE(camera_noise)) "sd_eps",
                "bp", "sum_obs", "sum_sim")

  fit <- run_encounter_rate_model(
    method = "TTE",
    code = code,
    constants = const,
    data = data_list,
    monitors = monitors,
    iter = iter,
    burnin = burnin,
    thin = thin,
    n_chains = n_chains,
    D_max = D_max,
    density_prior = density_prior,
    max_adapt_rounds = max_adapt_rounds,
    compute_WAIC = compute_WAIC,
    status_callback = status_callback
  )
  fit$settings$camera_noise <- isTRUE(camera_noise)
  fit
}

# -------------------------------------------------------------------
# USCR: buffer sizing from sigma (movement/space-use scale)
# -------------------------------------------------------------------
# Buffer rule: the state-space buffer must be at least uscr_buffer_multiplier
# times the radius of a circular home range containing uscr_buffer_p of an
# animal's activity. When the post-fit check finds the buffer too small,
# run_USCR() sets it uscr_buffer_safety_factor times above that minimum, so
# the next round passes if sigma barely changes.
uscr_buffer_p <- 0.99
uscr_buffer_multiplier <- 1.05
uscr_buffer_safety_factor <- 1.01

# Circular home-range radius (km) containing proportion p of an animal's
# activity, for a given sigma (km).
uscr_home_range_radius_km <- function(sigma_km, p = uscr_buffer_p) {
  sigma_km * sqrt(stats::qchisq(p, df = 2))
}

# The inverse: sigma (km) for a circular home range of area home_range_km2
# (km^2) containing proportion p of an animal's activity.
uscr_sigma_km_from_home_range <- function(home_range_km2, p) {
  sqrt(home_range_km2 / pi) / sqrt(stats::qchisq(p, df = 2))
}

# Suggested state-space buffer (m), given a sigma value (km) and a safety
# multiplier on the home-range radius. Used both to size the buffer up front
# (from the log(sigma) prior) and to check it after fitting (from the
# posterior). A bigger multiplier is safer but increases the state space
# (and run time).
uscr_buffer_from_sigma_km <- function(sigma_km, p = uscr_buffer_p,
                                      multiplier = uscr_buffer_multiplier) {
  multiplier * uscr_home_range_radius_km(sigma_km, p = p) * 1000
}

# Default buffer (m) computed from the log(sigma) prior before any data are
# fit: uses the p-quantile of the prior as a quick, conservative guess at
# sigma, then sizes the buffer from that guess via uscr_buffer_from_sigma_km().
uscr_default_buffer_m <- function(log_sigma_mean, log_sigma_sd, p = uscr_buffer_p,
                                  multiplier = uscr_buffer_multiplier) {
  sigma_q <- exp(log_sigma_mean + log_sigma_sd * stats::qnorm(p))
  uscr_buffer_from_sigma_km(sigma_q, p = p, multiplier = multiplier)
}

# -------------------------------------------------------------------
# USCR: state space, buffer, area
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
# plotting; `scale` converts km to the plot's units (1000 for meters).
uscr_buffered_outline <- function(cam, buffer_km, scale = 1) {
  xy <- sf::st_coordinates(uscr_buffered_union(cam, buffer_km))
  grp_cols <- intersect(c("L1", "L2", "L3"), colnames(xy))
  data.frame(
    x = xy[, "X"] * scale,
    y = xy[, "Y"] * scale,
    group = interaction(as.data.frame(xy[, grp_cols, drop = FALSE]), drop = TRUE)
  )
}

# USCR state space: the union of circles of radius buffer_m around the
# cameras, using their projected coordinates (utm_e, utm_n, in km), plus the
# bounding box for the activity-center priors and the area of the union.
uscr_state_space_and_area <- function(out, buffer_m) {
  cam <- cbind(out$utm_e, out$utm_n)
  buffer <- buffer_m / 1000
  list(
    buffer_sq = buffer^2,
    xlim = range(cam[, 1]) + c(-1, 1) * buffer,
    ylim = range(cam[, 2]) + c(-1, 1) * buffer,
    area_km2 = uscr_buffered_area_km2(cam, buffer)
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
check_density_prior <- function(lo, hi) {
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
  msg <- check_density_prior(D_lo, D_hi)
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

# -------------------------------------------------------------------
# REM / TTE: density prior via D = D_max * p_D, p_D ~ Beta(a, b)
# -------------------------------------------------------------------
# This is the continuous counterpart of the USCR prior above (where N / M has
# a Beta prior), so an informative interval gives all three models nearly the
# same prior on density. a = b = 1 is D ~ Uniform(0, D_max). D_max is an
# upper limit on density: it starts at 200 animals/km^2 under the uniform
# prior, or at 2x the interval's upper bound under an informative prior, and
# run_encounter_rate_model() doubles it (re-solving the prior) whenever the
# posterior for D presses against it.

encounter_prior_start_D_max <- function(D_hi) {
  2 * D_hi
}

uniform_beta_density_prior <- function(D_max) {
  list(uniform = TRUE, a = 1, b = 1, D_max = D_max)
}

# Solve for p_D ~ Beta(a, b) so that D = D_max * p_D has the requested 95%
# interval. Both bounds can usually be matched exactly; the mean and median
# follow from the interval and are reported back.
design_beta_density_prior <- function(D_lo, D_hi, D_max, miss_tolerance = 0.10) {
  msg <- check_density_prior(D_lo, D_hi)
  if (!is.null(msg)) stop(msg, call. = FALSE)
  if (D_hi >= D_max) {
    stop(sprintf("The density prior's upper bound (%.3g animals/km²) must be below D_max (%.3g).",
                 D_hi, D_max), call. = FALSE)
  }

  width <- D_hi - D_lo
  q_of <- function(p, k) D_max * stats::qbeta(c(0.025, 0.975), p * k, (1 - p) * k)
  miss_fn <- function(q) ((q[1] - D_lo) / width)^2 + ((q[2] - D_hi) / width)^2

  p_start <- min(max(((D_lo + D_hi) / 2) / D_max, 1e-6), 1 - 1e-6)
  k_start <- exp(stats::optimize(function(lk) miss_fn(q_of(p_start, exp(lk))),
                                 interval = log(c(0.01, 1e6)))$minimum)
  opt <- stats::optim(
    c(stats::qlogis(p_start), log(k_start)),
    function(par) miss_fn(q_of(stats::plogis(par[1]), exp(par[2]))),
    method = "Nelder-Mead",
    control = list(maxit = 2000, reltol = 1e-12)
  )
  p <- stats::plogis(opt$par[1])
  k <- exp(opt$par[2])
  a <- p * k
  b <- (1 - p) * k

  q <- q_of(p, k)
  achieved <- c(
    lo     = q[1],
    hi     = q[2],
    mean   = D_max * p,
    median = D_max * stats::qbeta(0.5, a, b)
  )
  miss <- c(lo = abs(q[1] - D_lo) / width, hi = abs(q[2] - D_hi) / width)

  list(
    uniform = FALSE,
    a = a,
    b = b,
    D_max = D_max,
    target = c(lo = D_lo, hi = D_hi),
    achieved = achieved,
    miss = miss,
    off_target = any(miss > miss_tolerance),
    miss_tolerance = miss_tolerance
  )
}

describe_beta_density_prior <- function(design) {
  if (is.null(design)) return(NULL)
  if (isTRUE(design$uniform)) {
    return(sprintf("Uniform (D ~ Uniform(0, D_max), D_max = %.3g animals/km²)", design$D_max))
  }
  sprintf(paste0(
    "Informative. Target 95%% interval %.3g to %.3g animals/km²; achieved ",
    "%.3g to %.3g (mean %.3g, median %.3g) with D = D_max * Beta(%.3g, %.3g), ",
    "D_max = %.3g"),
    design$target[["lo"]], design$target[["hi"]],
    design$achieved[["lo"]], design$achieved[["hi"]],
    design$achieved[["mean"]], design$achieved[["median"]],
    design$a, design$b, design$D_max)
}

# Plot the density priors for REM/TTE (curve) and, if available, USCR (bars:
# exact Beta-Binomial probabilities for N, rescaled to a density over
# animals/km^2 so both are on the same scale), with the target 95% interval
# (dashed lines).
plot_density_priors <- function(beta_design, uscr_design = NULL) {
  quiet_require("ggplot2")

  x_max <- beta_design$D_max * stats::qbeta(0.9995, beta_design$a, beta_design$b)
  if (!is.null(uscr_design)) {
    k <- 0:uscr_design$M
    p <- dbbinom(k, uscr_design$M, uscr_design$a, uscr_design$b)
    p <- p / sum(p)
    cdf <- cumsum(p)
    keep <- cdf >= 0.0005 & c(0, cdf[-length(cdf)]) <= 0.9995
    uscr_df <- data.frame(D = k[keep] / uscr_design$area_km2,
                          dens = p[keep] * uscr_design$area_km2,
                          model = "USCR")
    x_max <- max(x_max, max(uscr_df$D))
  }

  x <- seq(0, x_max, length.out = 400)
  beta_df <- data.frame(
    D = x,
    dens = stats::dbeta(x / beta_design$D_max, beta_design$a, beta_design$b) / beta_design$D_max,
    model = "REM & TTE"
  )
  beta_df <- beta_df[is.finite(beta_df$dens), ]

  colors <- c("REM & TTE" = "#8C4F2D", "USCR" = "#4F6FA3")
  fig <- ggplot2::ggplot()
  if (!is.null(uscr_design)) {
    fig <- fig +
      ggplot2::geom_col(
        data = uscr_df,
        ggplot2::aes(x = D, y = dens, fill = model),
        width = 1.02 / uscr_design$area_km2,  # slight overlap avoids hairline gaps
        position = "identity"
      ) +
      ggplot2::scale_fill_manual(values = c("USCR" = "#B9C6DF"), name = NULL)
  }
  fig <- fig +
    ggplot2::geom_line(data = beta_df, ggplot2::aes(x = D, y = dens, color = model), linewidth = 1) +
    ggplot2::scale_color_manual(values = colors, name = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(
      x = "Density (animals/km²)",
      y = "Prior density",
      title = "Prior on density implied by your interval",
      subtitle = "Dashed lines: your 95% interval."
    )
  if (!is.null(beta_design$target)) {
    fig <- fig + ggplot2::geom_vline(xintercept = beta_design$target[c("lo", "hi")],
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
      # Camera-level noise (Poisson-lognormal) only when use_eps = 1. With
      # use_eps = 0, eps and sd_eps are sampled from their priors but don't
      # affect the data, and they are not monitored.
      log(mu[j]) <- log(Lambda[j]) + log(days_per_cam[j]) + use_eps * eps[j]

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
                                psi_a,
                                psi_b,
                                camera_noise) {
  ss <- uscr_state_space_and_area(out, buffer_m = buffer_m)

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
    psi_b = psi_b,
    use_eps = as.numeric(isTRUE(camera_noise))
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
# Simulate data tab: compare posterior draws to the simulated truth
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

# Posterior predictive check. The Bayesian p-value is the posterior mean of
# bp = step(sum_sim - sum_obs), comparing a Pearson discrepancy for the
# observed counts with one for counts simulated from the fitted model. Values
# near 0 mean the observed counts vary more than the model expects (e.g.
# overdispersion); values near 1 mean they vary less. The check "fails" when
# the p-value falls outside pass_range (lower and upper limits).
ppc_summary <- function(fit, pass_range) {
  if (is.null(fit) || is.null(fit$samples_all) || !"bp" %in% colnames(fit$samples_all)) {
    return(NULL)
  }
  p <- mean(fit$samples_all[, "bp"], na.rm = TRUE)
  lower <- pass_range[1]
  upper <- pass_range[2]
  list(p = p, fail = is.finite(p) && (p < lower || p > upper), lower = lower, upper = upper)
}

# Truth-vs-estimate table for a model fit to simulated data.
# truth: named list, e.g. list(D_km2 = 25, sigma = 0.217, lam_0 = 1)
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
# -------------------------------------------------------------------
# Fits USCR, then reruns it while any check fails, up to max_adapt_rounds
# rounds: iterations, burn-in, and thinning double if Rhat is above
# rhat_target; M doubles if the posterior for N presses against it; and the
# state-space buffer grows if it is too small for the largest posterior sigma.

# The app supplies every setting without a default here (its constants live
# in app.R). The checks' thresholds and the buffer rule keep their defaults.
run_USCR <- function(out,
                     camera_counts,
                     camera_days,
                     iter,
                     burnin,
                     thin,
                     n_chains,
                     M,
                     log_sigma_mean,
                     log_sigma_sd,
                     log_lam0_mean,
                     log_lam0_sd,
                     sd_eps_shape,
                     sd_eps_rate,
                     buffer_m,
                     density_prior,
                     camera_noise,
                     max_adapt_rounds,
                     compute_WAIC,
                     status_callback,
                     rhat_target = 1.1,
                     M_check_fraction = 0.9,
                     M_check_prob = 0.01,
                     buffer_p = uscr_buffer_p,
                     buffer_multiplier = uscr_buffer_multiplier,
                     buffer_safety_factor = uscr_buffer_safety_factor) {

  quiet_require("nimble")

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
        sd_eps_rate = sd_eps_rate,
        camera_noise = isTRUE(camera_noise)
      )
    )
  }

  iter <- as.integer(iter)
  burnin <- as.integer(burnin)
  thin <- as.integer(max(1, thin))
  n_chains <- as.integer(max(1, n_chains))
  M <- as.integer(max(1, M))

  # density_prior: list(uniform = TRUE) for the default uniform prior, or
  # list(uniform = FALSE, lo = , hi = ): a 95% interval in animals/km^2 (the
  # app checks the interval before running).
  informative_prior <- identical(density_prior$uniform, FALSE)
  prior_design <- NULL

  # sd_eps is only meaningful (and only monitored) with camera-level noise on.
  monitors <- c(
    "log_sigma", "log_lam_0", "psi", if (isTRUE(camera_noise)) "sd_eps",
    "sigma", "lam_0", "N", "D",
    "sum_obs", "sum_sim", "bp"
  )

  code <- build_uscr_code()

  current_M <- M
  current_iter <- iter
  current_burnin <- burnin
  current_thin <- thin
  current_buffer_m <- as.numeric(buffer_m)
  round_log <- list()
  current_rhat <- Inf
  round_i <- 0L
  hit_max_adapt_rounds <- FALSE
  unmet_checks <- character(0)
  max_adapt_rounds <- as.integer(max(1, max_adapt_rounds))
  prev_elapsed_min <- NULL
  prev_iter <- NULL
  prev_M <- NULL
  fit <- NULL
  samples_list <- NULL
  samples_all <- NULL
  const <- NULL

  repeat {
    round_i <- round_i + 1L

    # Re-solve the density prior for this round's state-space area and M, so
    # the prior on density stays fixed as the buffer and M change.
    round_area_km2 <- uscr_state_space_and_area(out, buffer_m = current_buffer_m)$area_km2
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
      psi_b = prior_design$b,
      camera_noise = camera_noise
    )

    data_list <- list(
      y = as.numeric(camera_counts),
      in_ss = rep(1L, const$M)
    )

    status_callback(
      round = round_i,
      iter = current_iter,
      burnin = current_burnin,
      M = const$M,
      buffer_m = current_buffer_m,
      prev_elapsed_min = prev_elapsed_min,
      prev_iter = prev_iter,
      prev_M = prev_M
    )
    round_started <- Sys.time()

    fit <- run_nimble_chains(
      code = code,
      constants = const,
      data = data_list,
      monitors = monitors,
      niter = current_iter,
      nburnin = current_burnin,
      thin = current_thin,
      n_chains = n_chains,
      compute_WAIC = compute_WAIC,
      make_inits = make_uscr_inits
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
    prev_M <- const$M

    if (converged && !M_too_small && !buffer_too_small) {
      break
    }

    # Stop here (before growing iter/M/buffer) so the reported settings match
    # the last fit.
    if (round_i >= max_adapt_rounds) {
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
    round_history = do.call(rbind, round_log),
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
