# R/data_checks.R
# Cleaning, QC, trimming, and summary helpers for NPS camera data

# -------------------------------------------------------------------
# PRE-QC CLEANING HELPERS
# -------------------------------------------------------------------

parse_timestamp_robust <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = tz))
  }
  if (inherits(x, "Date")) {
    return(as.POSIXct(x, tz = tz))
  }
  
  x_chr <- trimws(as.character(x))
  x_chr[x_chr %in% c("", "NA", "NaN")] <- NA_character_
  
  parsed <- suppressWarnings(
    lubridate::parse_date_time(
      x_chr,
      orders = c(
        "Ymd HMS", "Ymd HM", "Ymd",
        "mdY HMS", "mdY HM", "mdY",
        "mdy HMS", "mdy HM", "mdy",
        "m/d/y HMS", "m/d/y HM", "m/d/y",
        "m/d/Y HMS", "m/d/Y HM", "m/d/Y"
      ),
      tz = tz,
      quiet = TRUE
    )
  )
  
  as.POSIXct(parsed, tz = tz)
}

# Strict mm/dd/yyyy parser for deployment dates. as.Date(x, format = "%m/%d/%Y")
# on its own is too forgiving: "6/15/24" silently becomes the year 0024 and
# trailing junk ("06/15/2024abc") is ignored, so both would pass validation and
# then break camera-day and window calculations later. This requires the whole
# value to be m/d/yyyy (1–2 digit month/day, 4-digit year) and returns NA otherwise.
parse_mdY_strict <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(format(x, "%Y-%m-%d")))
  x <- trimws(as.character(x))
  ok <- !is.na(x) & grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", x)
  out <- rep(as.Date(NA), length(x))
  out[ok] <- as.Date(x[ok], format = "%m/%d/%Y")
  out
}

clean_deployment_import <- function(deployments) {
  # Standardize column names: strip dots & weird whitespace, collapse spaces
  nm <- names(deployments)
  nm <- gsub("\\.+$", "", nm)
  nm <- gsub("[\\s\\p{Z}]+$", "", nm, perl = TRUE)
  nm <- gsub("^[\\s\\p{Z}]+", "", nm, perl = TRUE)
  nm <- gsub("[\\s\\p{Z}]+", " ", nm, perl = TRUE)
  nm <- trimws(nm)
  names(deployments) <- nm
  
  # Trim leading/trailing whitespace from all character columns
  deployments[] <- lapply(deployments, function(x) {
    if (is.character(x)) {
      gsub("^[\\s\\p{Z}]+|[\\s\\p{Z}]+$", "", x, perl = TRUE)
    } else x
  })
  
  # Camera Functioning: accept common variants → Yes / No (QC expects these)
  if ("Camera Functioning" %in% names(deployments)) {
    x <- deployments$`Camera Functioning`
    if (is.logical(x)) {
      deployments$`Camera Functioning` <- ifelse(x, "Yes", "No")
    } else {
      xc  <- trimws(as.character(x))
      low <- tolower(xc)
      yn <- xc
      yes_i <- low %in% c("yes", "y", "true", "t", "1") | xc %in% c("TRUE", "1")
      no_i  <- low %in% c("no", "n", "false", "f", "0") | xc %in% c("FALSE", "0")
      yn[yes_i] <- "Yes"
      yn[no_i]  <- "No"
      deployments$`Camera Functioning` <- yn
    }
  }
  
  deployments
}

clean_images_import <- function(images) {
  # Check for NULL or empty input
  if (is.null(images)) {
    stop("clean_images_import: images is NULL")
  }
  if (!is.data.frame(images)) {
    stop("clean_images_import: images must be a data frame")
  }
  if (nrow(images) == 0) {
    warning("clean_images_import: images data frame is empty")
    return(images)
  }
  
  # Standardize column names
  nm <- names(images)
  if (is.null(nm) || length(nm) == 0) {
    stop("clean_images_import: images has no column names")
  }
  nm <- gsub("\\.+$", "", nm)
  nm <- trimws(nm)
  names(images) <- nm
  
  # Strip explicit " UTC" suffix if present in Timestamp, then parse. The
  # column is only converted when EVERY non-blank value parses. Previously it
  # was converted if any value parsed, which turned unparseable timestamps into
  # NA here — check_images() then reported them as "Timestamp missing" instead
  # of "Bad Timestamp format", and the original text was lost. Leaving the
  # column as text in that case lets check_images() report them correctly.
  if ("Timestamp" %in% names(images) && !inherits(images$Timestamp, "POSIXt")) {
    raw_ts <- sub(" UTC$", "", trimws(as.character(images$Timestamp)))
    parsed_ts <- parse_timestamp_robust(raw_ts)
    has_input <- !is.na(raw_ts) & !raw_ts %in% c("", "NA", "NaN")
    if (any(has_input) && all(!is.na(parsed_ts[has_input]))) {
      images$Timestamp <- parsed_ts
    } else {
      images$Timestamp <- raw_ts
    }
  }
  
  # Split multi-species rows like "Deer|Raccoon" into separate rows, pairing
  # each Species entry with the Sighting Count entry in the same position
  # (e.g. Species = "Deer|Squirrel", Sighting Count = "1|4" becomes one row
  # for 1 deer and one row for 4 squirrels; all other columns are repeated
  # unchanged on both new rows).
  if (all(c("Species", "Sighting Count") %in% names(images))) {
    images <- images %>%
      dplyr::mutate(`Sighting Count` = as.character(`Sighting Count`))
    
    # tidyr::separate_rows() hard-errors with an unhelpful "can't recycle"
    # message if a row's Species and Sighting Count split into different
    # numbers of pipe-delimited pieces, so check for that explicitly first
    # and only for rows that actually use the pipe-delimited feature (a
    # single, non-piped Species with a blank Sighting Count is a separate,
    # already-handled "missing field" issue, not a mismatch).
    has_pipe <- grepl("|", images$Species, fixed = TRUE) |
      grepl("|", images$`Sighting Count`, fixed = TRUE)
    has_pipe[is.na(has_pipe)] <- FALSE
    if (any(has_pipe)) {
      species_n <- lengths(strsplit(images$Species[has_pipe], "|", fixed = TRUE))
      count_n   <- lengths(strsplit(images$`Sighting Count`[has_pipe], "|", fixed = TRUE))
      mismatch_rows <- which(has_pipe)[species_n != count_n]
      if (length(mismatch_rows) > 0) {
        stop(
          "clean_images_import: Species and Sighting Count have a different number ",
          "of pipe-delimited values in row(s) ", paste(mismatch_rows, collapse = ", "),
          ". Each pipe-delimited Species entry needs a matching Sighting Count entry ",
          "(for example Species = 'Deer|Squirrel' needs Sighting Count = '1|4')."
        )
      }
    }
    
    images <- images %>%
      tidyr::separate_rows(Species, `Sighting Count`, sep = "\\|")
    
    # Sighting Count must be numeric. Checked here, right after splitting and
    # before the as.numeric() coercion below, while the original text is
    # still available to report — as.numeric() would otherwise silently turn
    # a value like "unknown" into NA, losing the offending value and letting
    # build_nps_model_inputs() silently produce -Inf/NA detection counts
    # downstream instead of failing clearly here.
    bad_count <- which(
      !is.na(images$`Sighting Count`) & images$`Sighting Count` != "" &
        is.na(suppressWarnings(as.numeric(images$`Sighting Count`)))
    )
    if (length(bad_count) > 0) {
      stop(
        "clean_images_import: Sighting Count has non-numeric value(s) in row(s) ",
        paste(bad_count, collapse = ", "), ": ",
        paste(images$`Sighting Count`[bad_count], collapse = ", "), ". ",
        "Sighting Count must be a number (or pipe-delimited numbers for multi-species rows, ",
        "for example '1|4')."
      )
    }
    
    images <- images %>%
      dplyr::mutate(`Sighting Count` = as.numeric(`Sighting Count`))
  }
  
  images
}

# -------------------------------------------------------------------
# QC: check_deployments (as before)
# -------------------------------------------------------------------

check_deployments <- function(deployment, images = NULL, hemisphere = "Western",
                              lat_hemisphere = "Northern") {
  issues <- list()
  hemisphere <- if (identical(hemisphere, "Eastern")) "Eastern" else "Western"
  lat_hemisphere <- if (identical(lat_hemisphere, "Southern")) "Southern" else "Northern"
  
  # Rows with a missing/unparseable Start Date or End Date, or a missing/
  # invalid Latitude or Longitude, are collected here and treated as a hard
  # stop after the per-row loop (see below) rather than a warning, because
  # every model run depends on these four fields: Start/End Date drive the
  # deployment-window trim and the detection-matrix date range in
  # build_nps_model_inputs(), and Latitude/Longitude drive the UTM
  # projection used by every model's spatial/distance calculations.
  blocking_issues <- character()
  
  # ---- Column presence ----
  required_cols <- c(
    "Site Name",
    "Start Date", "Start Time", "End Date", "End Time",
    "Latitude", "Longitude", "Camera Functioning",
    "Camera Malfunction Date", "Detection Distance"
  )
  
  missing_cols <- setdiff(required_cols, names(deployment))
  if (length(missing_cols) > 0) {
    stop(
      "Deployment file is missing required column(s): ", paste(missing_cols, collapse = ", "), ". ",
      "If these columns are truly missing from your dataset (for example, no Start Time or End Time was ",
      "ever recorded), add the column names to your CSV file and leave the cells blank. If you believe ",
      "these columns exist but are labeled slightly differently (e.g. 'Start_Date' instead of 'Start Date'), ",
      "rename them to match exactly. After fixing the file, re-upload it.",
      call. = FALSE
    )
  }
  
  # ---- Site Name uniqueness (hard stop) ----
  # A duplicate Site Name would silently collide downstream (camera counts,
  # coordinate lookups, the trap array used by USCR, etc.), so this is a hard
  # error rather than a warning like the checks below.
  site_names <- trimws(deployment$`Site Name`)
  dup_sites <- unique(site_names[duplicated(site_names) & !is.na(site_names) & site_names != ""])
  if (length(dup_sites) > 0) {
    stop(
      "Deployment file has duplicate Site Name value(s): ", paste(dup_sites, collapse = ", "), ". ",
      "Each row must have a unique Site Name. If two cameras were deployed at the same physical site at ",
      "different times, give them distinct Site Names (for example, add a suffix). If this is a duplicate ",
      "row from a data-entry mistake, remove the extra row(s). After fixing the file, re-upload it.",
      call. = FALSE
    )
  }
  
  # ---- Per-row checks ----
  deployment_cols <- setdiff(required_cols, "Site Name")
  
  for (i in seq_len(nrow(deployment))) {
    row_values <- deployment[i, deployment_cols]
    
    # Skip row if all deployment columns except Site Name and Notes are blank
    if (all(is.na(row_values) | row_values == "")) next
    
    # --- Dates ---
    date_cols <- c("Start Date", "End Date")
    for (col in date_cols) {
      val <- deployment[[col]][i]
      if (is.na(val) || val == "") {
        msg <- paste("❌", col, "missing in row", i)
        issues <- c(issues, msg)
        blocking_issues <- c(blocking_issues, msg)
      } else if (is.na(parse_mdY_strict(val))) {
        msg <- paste("❌ Bad date in", col, "row", i, ":", val, " — should be mm/dd/yyyy (4-digit year)")
        issues <- c(issues, msg)
        blocking_issues <- c(blocking_issues, msg)
      }
    }
    
    # End Date must not be before Start Date. build_nps_model_inputs() maps
    # each camera's Start/End Date to column indices in a shared date matrix,
    # so a swapped/backward date range would otherwise silently produce a
    # collapsed or wrong deployment window for that camera instead of an
    # error.
    start_parsed <- parse_mdY_strict(deployment$`Start Date`[i])
    end_parsed   <- parse_mdY_strict(deployment$`End Date`[i])
    if (!is.na(start_parsed) && !is.na(end_parsed) && end_parsed < start_parsed) {
      msg <- paste("❌ End Date is before Start Date in row", i, ":",
                   deployment$`Start Date`[i], "→", deployment$`End Date`[i])
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    }
    
    # Camera Malfunction Date is required whenever Camera Functioning = No.
    # This only depends on columns already in the deployment file, so it
    # doesn't need images to be uploaded first.
    cam_func <- deployment$`Camera Functioning`[i]
    if (!is.na(cam_func) && tolower(cam_func) == "no") {
      val <- deployment$`Camera Malfunction Date`[i]
      if (is.na(val) || val == "") {
        issues <- c(issues, paste("❌ Camera Malfunction Date missing in row", i,
                                  " — required because Camera Functioning = No"))
      } else if (is.na(parse_mdY_strict(val))) {
        issues <- c(issues, paste("❌ Bad date in Camera Malfunction Date row", i, ":", val,
                                  " — should be mm/dd/yyyy (4-digit year). End Date was NOT replaced",
                                  "for this row."))
      } else {
        # The malfunction date replaces End Date below, so it should fall
        # within the deployment; otherwise the camera's window gets stretched
        # or inverted.
        mal_parsed <- parse_mdY_strict(val)
        if ((!is.na(start_parsed) && mal_parsed < start_parsed) ||
            (!is.na(end_parsed) && mal_parsed > end_parsed)) {
          issues <- c(issues, paste("⚠️ Camera Malfunction Date in row", i, ":", val,
                                    "is outside the Start Date–End Date range",
                                    paste0("(", deployment$`Start Date`[i], " – ",
                                           deployment$`End Date`[i], ")")))
        }
      }
    }
    
    # --- Times ---
    time_cols <- c("Start Time", "End Time")
    for (col in time_cols) {
      val <- deployment[[col]][i]
      if (is.na(val) || val == "") {
        issues <- c(issues, paste("❌", col, "missing in row", i))
      } else if (!grepl("^(?:[01]?[0-9]|2[0-3]):[0-5][0-9](:[0-5][0-9])?$", val)) {
        issues <- c(issues, paste("❌ Bad time in", col, "row", i, ":", val, " — should be HH:MM or HH:MM:SS 24h"))
      }
    }
    
    # ---- Numeric checks with suppression ----
    dd_val <- deployment$`Detection Distance`[i]
    if (is.na(dd_val) || dd_val == "") {
      issues <- c(issues, paste("❌ Detection Distance missing in row", i))
    } else if (suppressWarnings(is.na(as.numeric(dd_val)))) {
      issues <- c(issues, paste("❌ Non-numeric value in Detection Distance row", i, ":", dd_val))
    }
    
    if ("Camera Detection Angle" %in% names(deployment)) {
      theta_val <- deployment$`Camera Detection Angle`[i]
      theta_num <- suppressWarnings(as.numeric(theta_val))
      if (!is.na(theta_val) && theta_val != "") {
        if (is.na(theta_num)) {
          issues <- c(issues, paste("❌ Non-numeric value in Camera Detection Angle row", i, ":", theta_val))
        } else if (theta_num <= 0 || theta_num > 360) {
          issues <- c(issues, paste("❌ Camera Detection Angle out of range in row", i, ":", theta_val, "— must be between 0 and 360 degrees"))
        }
      }
    }
    
    if ("Camera Height" %in% names(deployment)) {
      ch_val <- deployment$`Camera Height`[i]
      if (!is.na(ch_val) && ch_val != "" && suppressWarnings(is.na(as.numeric(ch_val)))) {
        issues <- c(issues, paste("❌ Non-numeric value in Camera Height row", i, ":", ch_val))
      }
    }
    
    # --- Latitude check (must exist and be non-zero; auto-fix sign to match
    # the selected hemisphere in Model settings) ---
    lat_val <- suppressWarnings(as.numeric(deployment$Latitude[i]))
    if (is.na(lat_val) || lat_val == 0) {
      msg <- paste("❌ Latitude missing or invalid (0) in row", i)
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    } else {
      fixed_lat <- if (identical(lat_hemisphere, "Southern")) -abs(lat_val) else abs(lat_val)
      if (fixed_lat != lat_val) {
        deployment$Latitude[i] <- fixed_lat
        message("🛠 Fixed Latitude in row ", i, ": ", lat_val, " → ", fixed_lat)
      }
    }
    
    # --- Longitude check (must exist and be non-zero; auto-fix sign to match
    # the selected hemisphere in Model settings) ---
    long_val <- suppressWarnings(as.numeric(deployment$Longitude[i]))
    if (is.na(long_val) || long_val == 0) {
      msg <- paste("❌ Longitude missing or invalid (0) in row", i)
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    } else {
      fixed_val <- if (identical(hemisphere, "Eastern")) abs(long_val) else -abs(long_val)
      if (fixed_val != long_val) {
        deployment$Longitude[i] <- fixed_val
        message("🛠 Fixed Longitude in row ", i, ": ", long_val, " → ", fixed_val)
      }
    }
    
    # --- Camera Orientation ---
    if ("Camera Orientation" %in% names(deployment)) {
      val <- deployment$`Camera Orientation`[i]
      valid_cardinals <- c("N","NE","E","SE","S","SW","W","NW")
      if (!is.na(val) && val != "") {
        val_upper <- toupper(val)
        val_numeric <- suppressWarnings(as.numeric(val))
        if (!(val_upper %in% valid_cardinals | (!is.na(val_numeric) & val_numeric >= 0 & val_numeric <= 359))) {
          issues <- c(issues, paste("❌ Invalid Camera Orientation in row", i, ":", val,
                                    "— must be N/NE/.../NW or 0–359 degrees"))
        }
      }
    }
    
    # --- Camera Functioning ---
    val <- deployment$`Camera Functioning`[i]
    if (is.na(val) || val == "") {
      issues <- c(issues, paste("❌ Camera Functioning missing in row", i))
    } else {
      val_lower <- tolower(val)
      if (!(val_lower %in% c("yes","no"))) {
        issues <- c(issues, paste("❌ Invalid Camera Functioning value in row", i, ":", val, " — must be Yes or No"))
      }
    }
    
  } # end row loop
  
  # ---- Duplicate coordinates check (warning) ----
  # Two different Site Names sharing the exact same Latitude AND Longitude
  # usually means a copy-paste error created two records for what's actually
  # one physical camera location. It's fine for Latitude values or Longitude
  # values to repeat individually (e.g. cameras on the same north-south
  # transect) — only a duplicated (Latitude, Longitude) pair is flagged. Runs
  # after the per-row loop so it checks the sign-corrected coordinates.
  lat_num_all <- suppressWarnings(as.numeric(deployment$Latitude))
  lon_num_all <- suppressWarnings(as.numeric(deployment$Longitude))
  has_coords  <- !is.na(lat_num_all) & !is.na(lon_num_all)
  coord_key   <- paste(lat_num_all, lon_num_all)
  dup_keys    <- unique(coord_key[has_coords][duplicated(coord_key[has_coords])])
  for (k in dup_keys) {
    dup_rows <- which(has_coords & coord_key == k)
    issues <- c(issues, paste0(
      "⚠️ Duplicate coordinates in ", length(dup_rows), " row(s) (Site Name: ",
      paste(deployment$`Site Name`[dup_rows], collapse = ", "),
      ") — same Latitude/Longitude; verify these are meant to be different camera locations."
    ))
  }
  
  # ---- Hard stop: missing/invalid Start Date, End Date, Latitude, or
  # Longitude ----
  # Unlike the row-level issues above (which are surfaced as warnings so
  # users can see everything at once and decide whether to proceed), these
  # four fields are load-bearing for every model run, so a missing or
  # invalid value blocks upload entirely rather than silently propagating
  # into the deployment-window trim or the UTM projection.
  if (length(blocking_issues) > 0) {
    stop(
      "Deployment file has missing or invalid Start Date, End Date, Latitude, and/or Longitude ",
      "value(s) in ", length(blocking_issues), " place(s):\n",
      paste0("  - ", blocking_issues, collapse = "\n"),
      "\nStart Date and End Date must be present, in mm/dd/yyyy format, and End Date must be on or after ",
      "Start Date for every row; they set each camera's deployment window. Latitude and Longitude must be ",
      "present and non-zero for every row; ",
      "they're used to project camera locations for every model. After fixing the file, re-upload it.",
      call. = FALSE
    )
  }
  
  
  #Replace End Date with Camera Malfunction Date for Site Names where Camera Functioning is "No"
  # The original End Date (when the camera was picked up) is kept in
  # `Original End Date` first. check_images() needs it: a malfunctioning
  # camera usually keeps taking (unusable) pictures until pickup, so those
  # images are expected and are NOT a clock problem. Only created if absent,
  # because check_deployments() runs a second time on already-checked data
  # when images are uploaded, and by then End Date has been overwritten.
  if (!"Original End Date" %in% names(deployment)) {
    deployment$`Original End Date` <- deployment$`End Date`
  }
  # Only valid mm/dd/yyyy malfunction dates are used. Previously any non-blank
  # value was copied in, so a badly formatted malfunction date (which is only a
  # warning above) would overwrite a good End Date with an unparseable one
  # *after* the End Date hard-stop had already passed.
  malfunction_idx <- which(
    tolower(deployment$`Camera Functioning`) == "no" &
      !is.na(parse_mdY_strict(deployment$`Camera Malfunction Date`))
  )
  
  deployment$`End Date`[malfunction_idx] <-
    deployment$`Camera Malfunction Date`[malfunction_idx]
  
  # ---- Output ----
  if (length(issues) == 0) {
    message("✅ Deployments file is formatted correctly!")
  } else {
    warning("⚠️ Issues found in deployments file:\n", paste(issues, collapse = "\n"))
  }
  return(deployment)
}

# -------------------------------------------------------------------
# QC: check_images (as before)
# -------------------------------------------------------------------

check_images <- function(images, deployments, survey_year = NULL, hemisphere = "Western",
                         lat_hemisphere = "Northern", grace_days = 1) {
  # survey_year is retained for backward compatibility with older callers.
  # The current app relies on parsed image timestamps instead.
  issues <- c()
  fixes  <- c()
  hemisphere <- if (identical(hemisphere, "Eastern")) "Eastern" else "Western"
  lat_hemisphere <- if (identical(lat_hemisphere, "Southern")) "Southern" else "Northern"
  
  # ---- Required columns ----
  required_cols <- c("Site Name", "Timestamp", "Species", "Sighting Count", "Cluster ID")
  missing_cols <- setdiff(required_cols, names(images))
  if (length(missing_cols) > 0) {
    stop(paste0(
      "❌ Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      "\nPlease fix the images file and re-run the function."
    ))
  }
  
  # ---- Trim whitespace from character columns ----
  char_cols <- names(images)[sapply(images, is.character)]
  images[char_cols] <- lapply(images[char_cols], trimws)
  
  
  # ---- Latitude / Longitude checks (optional in image files) ----
  if (all(c("Latitude", "Longitude") %in% names(images))) {
    lat_num <- suppressWarnings(as.numeric(images$Latitude))
    lon_num <- suppressWarnings(as.numeric(images$Longitude))
    
    bad_lat <- which(is.na(lat_num) | lat_num == 0)
    if (length(bad_lat) > 0) {
      issues <- c(issues, paste0("❌ ", length(bad_lat), " image(s) have missing or zero Latitude"))
    }
    
    bad_lon <- which(is.na(lon_num) | lon_num == 0)
    if (length(bad_lon) > 0) {
      issues <- c(issues, paste0("❌ ", length(bad_lon), " image(s) have missing or zero Longitude"))
    }
    
    # Auto-fix latitude sign to match the selected hemisphere in Model settings
    fixed_lat <- if (identical(lat_hemisphere, "Southern")) -abs(lat_num) else abs(lat_num)
    fix_lat <- which(!is.na(lat_num) & lat_num != 0 & fixed_lat != lat_num)
    if (length(fix_lat) > 0) {
      images$Latitude[fix_lat] <- fixed_lat[fix_lat]
      message("🛠 Fixed Latitude for ", length(fix_lat), " image(s) ")
    }
    
    # Auto-fix longitude sign to match the selected hemisphere in Model settings
    fixed_lon <- if (identical(hemisphere, "Eastern")) abs(lon_num) else -abs(lon_num)
    fix_lon <- which(!is.na(lon_num) & lon_num != 0 & fixed_lon != lon_num)
    if (length(fix_lon) > 0) {
      images$Longitude[fix_lon] <- fixed_lon[fix_lon]
      message("🛠 Fixed Longitude for ", length(fix_lon), " image(s) ")
    }
  }
  
  # ---- Required fields check ----
  for (field in c("Species", "Sighting Count", "Cluster ID")) {
    missing_rows <- which(is.na(images[[field]]) | images[[field]] == "")
    if (length(missing_rows) > 0) {
      issues <- c(issues, paste0("❌ Missing ", field, " in ", length(missing_rows), " image(s)"))
    }
  }
  
  # Cluster ID may be any format (character or numeric) — it just needs to be
  # present, which the "Required fields check" above already covers. Multiple
  # rows legitimately share the same Cluster ID (e.g. burst photos, or a
  # multi-species image split into one row per species), so no uniqueness
  # check is applied here. It should, however, always identify detections at
  # a single camera: build_nps_model_inputs() groups by Site Name and
  # Cluster ID together, so a Cluster ID spanning more than one Site Name no
  # longer gets silently merged across cameras, but it usually still
  # signals a real problem with how Cluster IDs were assigned upstream
  # (e.g. an event-numbering scheme that isn't actually unique per camera),
  # so it's flagged here for the user to check.
  cluster_site_counts <- images %>%
    dplyr::filter(!is.na(`Cluster ID`) & `Cluster ID` != "") %>%
    dplyr::group_by(`Cluster ID`) %>%
    dplyr::summarise(n_sites = dplyr::n_distinct(`Site Name`), .groups = "drop")
  bad_clusters <- cluster_site_counts$`Cluster ID`[cluster_site_counts$n_sites > 1]
  if (length(bad_clusters) > 0) {
    issues <- c(issues, paste0(
      "❌ ", length(bad_clusters), " Cluster ID value(s) appear at more than one Site Name: ",
      paste(utils::head(bad_clusters, 10), collapse = ", "),
      if (length(bad_clusters) > 10) ", ..." else ""
    ))
  }
  
  # ---- Cross-check Site Names with deployments ----
  # This is a hard stop rather than a warning: build_nps_model_inputs() left-
  # joins the wide detection matrix onto deployment metadata by Site Name, so
  # an image Site Name with no matching deployment row would silently pick up
  # NA Start Date/End Date/Latitude/Longitude/Detection Distance and surface
  # as a confusing NA/NaN error deep inside the model instead of a clear
  # message here.
  # Blank Site Names are reported separately: they'd otherwise fall into the
  # "not found" list below and show up as a confusing "NA".
  missing_site_idx <- which(is.na(images$`Site Name`) | images$`Site Name` == "")
  bad_site_match <- setdiff(which(!images$`Site Name` %in% deployments$`Site Name`),
                            missing_site_idx)
  if (length(missing_site_idx) > 0 || length(bad_site_match) > 0) {
    bad_site_names <- unique(images$`Site Name`[bad_site_match])
    stop(
      if (length(missing_site_idx) > 0) paste0(
        "Image file has ", length(missing_site_idx), " row(s) with a blank Site Name ",
        "(first few rows: ", paste(utils::head(missing_site_idx, 10), collapse = ", "), "). "
      ) else "",
      if (length(bad_site_match) > 0) paste0(
        "Image file has ", length(bad_site_match), " row(s) with Site Name value(s) not found ",
        "in the deployment file: ", paste(utils::head(bad_site_names, 10), collapse = ", "),
        if (length(bad_site_names) > 10) ", ..." else "", ". "
      ) else "",
      "Every Site Name in the images file must exactly match a Site Name in the deployment file. ",
      "Fix the blank or mismatched Site Name(s) (check for typos or a missing deployment row), then re-upload.",
      call. = FALSE
    )
  }
  
  # ---- Timestamp parsing ----
  ts_parsed <- parse_timestamp_robust(images$Timestamp)
  ts_missing_input <- !is.na(images$Timestamp) & trimws(as.character(images$Timestamp)) != ""
  missing_ts_idx <- which(!ts_missing_input)
  if (length(missing_ts_idx) > 0) {
    issues <- c(issues, paste0("❌ Timestamp missing in ", length(missing_ts_idx), " image(s)"))
  }
  bad_ts <- which(ts_missing_input & is.na(ts_parsed))
  if (length(bad_ts) > 0) {
    site_counts <- table(images$`Site Name`[bad_ts])
    for (sn in names(site_counts)) {
      issues <- c(issues, paste0("❌ Bad Timestamp format at site ", sn,
                                 " — ", site_counts[[sn]], " occurrence(s); use a recognizable date-time like yyyy-mm-dd HH:MM:SS or mm/dd/yyyy HH:MM."))
    }
  }
  images$Timestamp <- ts_parsed
  
  # ---- Image dates vs deployment window ----
  # Deployment Start/End Times are approximate (they're whatever was written
  # down when someone walked up to the camera), so an image a few minutes or
  # hours past End Time is normal and not worth flagging. What this check is
  # for is timestamps that are genuinely wrong — usually a camera whose
  # date/clock was never set or was set to the wrong day/month/year.
  #
  # So the comparison is done on calendar DATES only and ignores Start/End
  # Time entirely (a blank End Time used to make the whole window NA, which
  # silently skipped the check for that camera). An image is flagged only if
  # its date is more than `grace_days` days before Start Date or after the
  # date the camera was picked up. With the default grace_days = 1, any image
  # on the pickup date, or the following day, passes.
  #
  # Malfunctioning cameras: check_deployments() replaces End Date with the
  # Camera Malfunction Date, but the camera often keeps firing until it's
  # picked up. Those images are expected, not a clock error, so the window
  # here uses `Original End Date` (the pickup date) when it exists. Images
  # between the malfunction date and pickup are reported as an informational
  # message instead, since the model already leaves them out (the camera's
  # End Index stops at the malfunction date in build_nps_model_inputs()).
  grace_days <- suppressWarnings(as.integer(grace_days[1]))
  if (is.na(grace_days) || grace_days < 0) grace_days <- 1L
  
  pickup_col <- if ("Original End Date" %in% names(deployments)) "Original End Date" else "End Date"
  
  deployment_windows <- data.frame(
    `Site Name` = deployments$`Site Name`,
    dep_start   = parse_mdY_strict(deployments$`Start Date`),
    dep_end     = parse_mdY_strict(deployments[[pickup_col]]),
    analysis_end = parse_mdY_strict(deployments$`End Date`),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Calendar date as shown on the timestamp (no time zone shifting).
  img_dates <- images %>%
    dplyr::filter(!is.na(Timestamp)) %>%
    dplyr::transmute(
      `Site Name`,
      img_date = as.Date(format(Timestamp, "%Y-%m-%d"))
    )
  
  window_check <- img_dates %>%
    dplyr::inner_join(deployment_windows, by = "Site Name") %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarise(
      dep_start    = dplyr::first(dep_start),
      dep_end      = dplyr::first(dep_end),
      analysis_end = dplyr::first(analysis_end),
      first_image  = min(img_date),
      last_image   = max(img_date),
      n_images     = dplyr::n(),
      n_before     = sum(img_date < dep_start - grace_days, na.rm = TRUE),
      n_after      = sum(img_date > dep_end   + grace_days, na.rm = TRUE),
      # Images after the malfunction date but not past pickup (+ grace), i.e.
      # expected images that the model will drop. Zero for working cameras,
      # where analysis_end == dep_end.
      n_post_malfunction = sum(img_date > analysis_end &
                                 img_date <= dep_end + grace_days, na.rm = TRUE),
      .groups = "drop"
    )
  
  fmt_d <- function(d) format(d, "%m/%d/%Y")
  
  for (i in seq_len(nrow(window_check))) {
    w <- window_check[i, ]
    
    if (w$n_before > 0) {
      issues <- c(issues, paste0(
        "⚠️ Site ", w$`Site Name`, ": ", w$n_before, " of ", w$n_images,
        " image(s) are dated more than ", grace_days, " day(s) before the deployment Start Date ",
        "(Start Date = ", fmt_d(w$dep_start), "; earliest image date = ", fmt_d(w$first_image), "). ",
        "Check that the camera's date/time was set correctly."
      ))
    }
    
    if (w$n_after > 0) {
      issues <- c(issues, paste0(
        "⚠️ Site ", w$`Site Name`, ": ", w$n_after, " of ", w$n_images,
        " image(s) are dated more than ", grace_days, " day(s) after the camera was picked up ",
        "(End Date = ", fmt_d(w$dep_end), "; latest image date = ", fmt_d(w$last_image), "). ",
        "Check that the camera's date/time was set correctly."
      ))
    }
    
    if (w$n_post_malfunction > 0) {
      message(
        "ℹ️ Site ", w$`Site Name`, ": ", w$n_post_malfunction,
        " image(s) are dated after the Camera Malfunction Date (", fmt_d(w$analysis_end),
        ") and will be excluded from the analysis."
      )
    }
  }
  
  # ---- Final output ----
  if (length(issues) == 0) {
    message("✅ Images file is formatted correctly and all Site Names match deployments!")
  } else {
    warning("⚠️ Issues found in images file:\n", paste(issues, collapse = "\n"))
  }
  
  return(images)
}

# -------------------------------------------------------------------
# FORMAT DEPLOYMENTS (optional trimming)
# -------------------------------------------------------------------

format_deployments <- function(deployments, max_days = NULL) {
  # Ensure dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for format_deployments()", call. = FALSE)
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required for format_deployments()", call. = FALSE)
  }
  
  # Parse dates (mm/dd/yyyy, mm/dd/yy, or yyyy-mm-dd).
  # The old version had an as.Date(tryFormats = ...) fallback that ran
  # whenever ANY value was NA (e.g. one blank row) and then re-parsed the
  # WHOLE column with the first format that fit the first value — "%m/%d/%y",
  # which reads "06/15/2024" as 2020-06-15 — throwing away the correct
  # lubridate results. It could also error outright. lubridate alone already
  # handles 2- and 4-digit years, so the fallback is removed.
  parse_date_safe <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXt")) return(as.Date(format(x, "%Y-%m-%d")))
    x <- trimws(as.character(x))
    x[x %in% c("", "NA")] <- NA_character_
    as.Date(lubridate::parse_date_time(x, orders = c("mdy", "ymd"), quiet = TRUE))
  }
  
  deps <- deployments %>%
    dplyr::mutate(
      Start_Date = parse_date_safe(`Start Date`),
      End_Date   = parse_date_safe(`End Date`)
    )
  
  deps <- deps %>%
    dplyr::mutate(
      operational_days = as.numeric(difftime(End_Date, Start_Date, units = "days"))
    )
  
  if (!is.null(max_days) && is.finite(max_days) && max_days > 0) {
    max_days <- as.integer(max_days[[1]])
    deps <- deps %>%
      dplyr::mutate(
        End_Date = dplyr::if_else(
          operational_days > max_days,
          Start_Date + max_days,
          End_Date
        )
      )
  }
  
  deps <- deps %>%
    dplyr::mutate(`End Date` = format(End_Date, "%m/%d/%Y")) %>%
    dplyr::select(-Start_Date, -End_Date, -operational_days)
  
  return(deps)
}

# -------------------------------------------------------------------
# SUMMARY HELPERS (deployment, images, deer)
# -------------------------------------------------------------------

summarize_deployments <- function(deployments) {
  deployments %>%
    dplyr::mutate(
      Start_Date = parse_mdY_strict(`Start Date`),
      End_Date   = parse_mdY_strict(`End Date`),
      operational_days = as.numeric(difftime(End_Date, Start_Date,
                                             units = "days"))
    ) %>%
    dplyr::summarise(
      n_cameras         = dplyr::n(),
      mean_days         = mean(operational_days, na.rm = TRUE),
      total_camera_days = sum(operational_days, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_images_by_species <- function(images) {
  images %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(
      total_images     = dplyr::n(),
      total_detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      cameras_detected = dplyr::n_distinct(`Site Name`),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_detections))
}

filter_species_rows <- function(images, species_name) {
  species_name <- trimws(as.character(species_name)[1])
  images %>%
    dplyr::filter(
      !is.na(Species),
      tolower(trimws(Species)) == tolower(species_name)
    )
}

species_summary_per_site <- function(images, species_name) {
  species_images <- filter_species_rows(images, species_name)
  
  species_images %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarise(
      images = dplyr::n(),
      detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )
}

first_finite_or_na <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals[1] else NA_real_
}

species_counts_per_camera <- function(images, species_name, deployments) {
  counts <- deployments %>%
    dplyr::left_join(
      filter_species_rows(images, species_name) %>%
        dplyr::group_by(`Site Name`) %>%
        dplyr::summarise(
          total_detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Site Name"
    ) %>%
    tidyr::replace_na(list(total_detections = 0))
  return(counts)
}

species_daily_detections <- function(images, species_name) {
  filter_species_rows(images, species_name) %>%
    dplyr::mutate(Date = as.Date(Timestamp)) %>%
    dplyr::group_by(`Site Name`, Date) %>%
    dplyr::summarise(
      detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )
}