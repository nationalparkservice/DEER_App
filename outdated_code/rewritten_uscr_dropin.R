# -------------------------------------------------------------------
# USCR: state space, buffer, area
# Drop-in replacement for the USCR section in sim_and_models.R
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

build_uscr_code <- function(J) {
  quiet_require("nimble")

  if (J == 1L) {
    nimble::nimbleCode({
      log_sigma ~ dnorm(log_sigma_mean, sd = log_sigma_sd)
      log_lam_0 ~ dnorm(0, sd = log_lam0_sd)
      psi ~ dunif(0, 1)
      sd_eps ~ dgamma(sd_eps_shape, sd_eps_shape)

      sigma <- exp(log_sigma)
      lam_0 <- exp(log_lam_0)

      for (i in 1:M) {
        z[i] ~ dbern(psi)

        hrc[i, 1] ~ dunif(xlim[1], xlim[2])
        hrc[i, 2] ~ dunif(ylim[1], ylim[2])

        dist2[i, 1] <- (hrc[i, 1] - cam[1, 1]) ^ 2 +
          (hrc[i, 2] - cam[1, 2]) ^ 2

        lambda[i, 1] <- z[i] * lam_0 *
          exp(-dist2[i, 1] / (2 * sigma ^ 2))

        in_ss[i] ~ dconstraint(dist2[i, 1] < buffer_sq)
      }

      eps[1] ~ dnorm(0, sd = sd_eps)
      Lambda[1] <- sum(lambda[1:M, 1])
      log(mu[1]) <- log(Lambda[1]) + log(day1) + eps[1]

      y[1] ~ dpois(mu[1])
      y_sim[1] ~ dpois(mu[1])

      pearson_obs[1] <- (y[1] - mu[1]) ^ 2 / mu[1]
      pearson_sim[1] <- (y_sim[1] - mu[1]) ^ 2 / mu[1]

      sum_obs <- pearson_obs[1]
      sum_sim <- pearson_sim[1]
      bp <- step(sum_sim - sum_obs)

      N <- sum(z[1:M])
      D_mi2 <- N / area_mi2
    })
  } else {
    nimble::nimbleCode({
      log_sigma ~ dnorm(log_sigma_mean, sd = log_sigma_sd)
      log_lam_0 ~ dnorm(0, sd = log_lam0_sd)
      psi ~ dunif(0, 1)
      sd_eps ~ dgamma(sd_eps_shape, sd_eps_shape)

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
}

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

make_uscr_constants <- function(out,
                                camera_days,
                                M,
                                log_sigma_mean,
                                log_sigma_sd,
                                log_lam0_sd,
                                sd_eps_shape,
                                buffer_chi_p,
                                buffer_scale) {
  ss <- uscr_state_space_and_area(
    out,
    buffer_chi_p = buffer_chi_p,
    buffer_scale = buffer_scale
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
    log_lam0_sd = log_lam0_sd,
    sd_eps_shape = sd_eps_shape
  )

  if (const$J == 1L) {
    const$day1 <- as.numeric(camera_days[1])
  } else {
    const$days_per_cam <- as.numeric(camera_days)
  }

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
  diag <- MCMCvis::MCMCsummary(samples_list)

  if (!"Rhat" %in% names(diag)) {
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
# Main USCR wrapper
# -------------------------------------------------------------------

run_USCR <- function(out,
                     camera_counts,
                     camera_days,
                     M = 300,
                     iter = NULL,
                     thin = NULL,
                     n_chains = NULL,
                     iter_tune = 2000,
                     burnin = 1000,
                     thin_tune = 1,
                     n_chains_tune = 1,
                     n_chains_final = 3,
                     final_iter = NULL,
                     final_thin = NULL,
                     parallel_chains = TRUE,
                     compute_WAIC = TRUE,
                     adapt = TRUE,
                     adaptive = NULL,
                     max_adapt_rounds = 4,
                     rhat_target = 1.1,
                     psi_threshold = 0.9,
                     psi_prob_cutoff = 0.01,
                     monitor_latent = FALSE,
                     diagnostic_mode = NULL,
                     tuning_n_chains = NULL,
                     log_sigma_mean = -1.4442,
                     log_sigma_sd = 0.1451,
                     log_lam0_sd = 1,
                     sd_eps_shape = 1,
                     buffer_chi_p = 0.99,
                     buffer_scale = 0.40,
                     seed = NULL,
                     verbose = TRUE,
                     ...) {

  # Compatibility shim: older and newer app wrappers use slightly
  # different argument names for the same USCR settings.
  if (!is.null(iter) && is.null(final_iter)) {
    final_iter <- as.integer(iter)
  }
  if (!is.null(thin) && is.null(final_thin)) {
    final_thin <- as.integer(thin)
  }
  if (!is.null(n_chains) && missing(n_chains_final)) {
    n_chains_final <- as.integer(n_chains)
  }
  if (!is.null(tuning_n_chains)) {
    n_chains_tune <- as.integer(tuning_n_chains)
  }
  if (!is.null(adaptive)) {
    adapt <- isTRUE(adaptive)
  }
  if (!is.null(diagnostic_mode)) {
    monitor_latent <- isTRUE(diagnostic_mode)
  }

  quiet_require("nimble")

  J <- nrow(out)
  if (length(camera_counts) != J || length(camera_days) != J) {
    stop(
      "run_USCR: camera_counts and camera_days must have length nrow(out).",
      call. = FALSE
    )
  }

  base_monitors <- c(
    "log_sigma", "log_lam_0", "psi", "sd_eps",
    "sigma", "lam_0", "N", "D_mi2",
    "sum_obs", "sum_sim", "bp"
  )
  latent_monitors <- c("z", "hrc", "eps")
  final_monitors <- if (isTRUE(monitor_latent)) {
    c(base_monitors, latent_monitors)
  } else {
    base_monitors
  }

  code <- build_uscr_code(J)

  current_M <- as.integer(M)
  current_iter <- as.integer(iter_tune)
  current_thin <- as.integer(thin_tune)
  adapt_log <- list()
  last_tune_fit <- NULL

  if (adapt) {
    for (round_i in seq_len(max_adapt_rounds)) {
      const <- make_uscr_constants(
        out = out,
        camera_days = camera_days,
        M = current_M,
        log_sigma_mean = log_sigma_mean,
        log_sigma_sd = log_sigma_sd,
        log_lam0_sd = log_lam0_sd,
        sd_eps_shape = sd_eps_shape,
        buffer_chi_p = buffer_chi_p,
        buffer_scale = buffer_scale
      )

      data_list <- list(
        y = as.numeric(camera_counts),
        in_ss = rep(1L, const$M)
      )

      if (verbose) {
        message(
          "USCR tuning round ", round_i,
          ": M=", const$M,
          ", iter=", current_iter,
          ", thin=", current_thin,
          ", chains=", n_chains_tune
        )
      }

      tune_fit <- run_uscr_chains(
        code = code,
        constants = const,
        data = data_list,
        monitors = base_monitors,
        niter = current_iter,
        nburnin = burnin,
        thin = current_thin,
        n_chains = n_chains_tune,
        parallel_chains = parallel_chains,
        compute_WAIC = FALSE,
        seed = if (is.null(seed)) NULL else seed + round_i
      )

      samples_list <- lapply(tune_fit, `[[`, "samples")
      samples_all <- do.call(rbind, samples_list)
      psi_post <- samples_all[, "psi"]

      rhat <- safe_rhat_max(samples_list)
      M_too_small <- mean(psi_post > psi_threshold, na.rm = TRUE) > psi_prob_cutoff
      converged <- is.na(rhat) || rhat <= rhat_target

      adapt_log[[round_i]] <- data.frame(
        round = round_i,
        M = const$M,
        niter = current_iter,
        nburnin = burnin,
        thin = current_thin,
        n_chains = n_chains_tune,
        rhat_max = rhat,
        M_too_small = M_too_small,
        stringsAsFactors = FALSE
      )

      last_tune_fit <- tune_fit

      if (converged && !M_too_small) {
        break
      }

      if (!converged) {
        current_iter <- current_iter * 2L
        current_thin <- max(1L, floor((current_iter - burnin) / 1000))
      } else if (M_too_small) {
        current_M <- current_M * 2L
      }
    }
  }

  if (is.null(final_iter)) {
    final_iter <- current_iter
  }
  if (is.null(final_thin)) {
    final_thin <- current_thin
  }

  final_const <- make_uscr_constants(
    out = out,
    camera_days = camera_days,
    M = current_M,
    log_sigma_mean = log_sigma_mean,
    log_sigma_sd = log_sigma_sd,
    log_lam0_sd = log_lam0_sd,
    sd_eps_shape = sd_eps_shape,
    buffer_chi_p = buffer_chi_p,
    buffer_scale = buffer_scale
  )

  final_data <- list(
    y = as.numeric(camera_counts),
    in_ss = rep(1L, final_const$M)
  )

  if (verbose) {
    message(
      "USCR final run: M=", final_const$M,
      ", iter=", final_iter,
      ", thin=", final_thin,
      ", chains=", n_chains_final,
      ", WAIC=", compute_WAIC,
      ", monitor_latent=", monitor_latent
    )
  }

  final_fit <- run_uscr_chains(
    code = code,
    constants = final_const,
    data = final_data,
    monitors = final_monitors,
    niter = final_iter,
    nburnin = burnin,
    thin = final_thin,
    n_chains = n_chains_final,
    parallel_chains = parallel_chains,
    compute_WAIC = compute_WAIC,
    seed = if (is.null(seed)) NULL else seed + 1000L
  )

  samples_list <- lapply(final_fit, `[[`, "samples")
  samples_all <- do.call(rbind, samples_list)

  list(
    method = "USCR",
    samples_list = samples_list,
    samples_all = samples_all,
    waic = extract_waic_mean(final_fit),
    fit_objects = final_fit,
    tuning_fit_objects = last_tune_fit,
    tuning_history = if (length(adapt_log) > 0L) do.call(rbind, adapt_log) else NULL,
    final_rhat_max = safe_rhat_max(samples_list),
    settings = list(
      M = final_const$M,
      iter_tune = iter_tune,
      final_iter = final_iter,
      burnin = burnin,
      thin_tune = thin_tune,
      final_thin = final_thin,
      n_chains_tune = n_chains_tune,
      n_chains_final = n_chains_final,
      compute_WAIC = compute_WAIC,
      monitor_latent = monitor_latent,
      adapt = adapt,
      max_adapt_rounds = max_adapt_rounds,
      rhat_target = rhat_target,
      psi_threshold = psi_threshold,
      psi_prob_cutoff = psi_prob_cutoff,
      buffer_chi_p = buffer_chi_p,
      buffer_scale = buffer_scale,
      area_mi2 = final_const$area_mi2
    )
  )
}

# -------------------------------------------------------------------
# Suggested fast call for app / simulation contexts
# -------------------------------------------------------------------
# uscr_fit <- run_USCR(
#   out = out,
#   camera_counts = camera_counts,
#   camera_days = camera_days,
#   M = 300,
#   iter_tune = 2000,
#   burnin = 1000,
#   thin_tune = 1,
#   n_chains_tune = 1,
#   n_chains_final = 3,
#   compute_WAIC = TRUE,
#   adapt = TRUE,
#   monitor_latent = FALSE,
#   seed = 123,
#   verbose = TRUE
# )
