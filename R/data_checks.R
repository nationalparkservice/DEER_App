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
  
  # Strip explicit " UTC" suffix if present in Timestamp
  if ("Timestamp" %in% names(images)) {
    images$Timestamp <- sub(" UTC$", "", images$Timestamp)
    parsed_ts <- parse_timestamp_robust(images$Timestamp)
    if (sum(!is.na(parsed_ts)) > 0) {
      images$Timestamp <- parsed_ts
    }
  }
  
  # Split multi-species rows like "Deer|Raccoon" into separate rows
  if (all(c("Species", "Sighting Count") %in% names(images))) {
    images <- images %>%
      dplyr::mutate(`Sighting Count` = as.character(`Sighting Count`)) %>%
      tidyr::separate_rows(Species, `Sighting Count`, sep = "\\|") %>%
      dplyr::mutate(`Sighting Count` = suppressWarnings(as.numeric(`Sighting Count`)))
  }
  
  images
}

# -------------------------------------------------------------------
# QC: check_deployments (as before)
# -------------------------------------------------------------------

check_deployments <- function(deployment, images = NULL) {
  issues <- list()
  
  # ---- Column presence ----
  required_cols <- c(
    "Site Name",
    "Start Date", "Start Time", "End Date", "End Time",
    "Latitude", "Longitude", "Camera Functioning",
    "Camera Malfunction Date", "Detection Distance"
  )
  
  missing_cols <- setdiff(required_cols, names(deployment))
  if (length(missing_cols) > 0) {
    cat(
      "❌ Missing required column(s):", paste(missing_cols, collapse = ", "), "\n\n",
      "💡 Guidance:\n",
      "• If these columns are truly missing from your dataset (for example, no Start Time or End Time was ever recorded),\n",
      "  please **add the column names** to your CSV file and leave the cells blank.\n",
      "• If you believe these columns exist but are labeled slightly differently (e.g., 'Start_Date' instead of 'Start Date'),\n",
      "  please rename them to match exactly the expected column names above.\n\n",
      "After fixing the file, re-run `check_deployments()`.\n"
    )
    stop("Deployment file is missing one or more required columns.")
  }
  
  # ---- Auto-fix Site Names missing leading zero ----
  deployment$`Site Name` <- vapply(deployment$`Site Name`, function(sn) {
    if (is.na(sn) || sn == "") return(sn)
    
    # Match patterns like HOFU_1 or FRSP_WILD_9 and add a leading zero.
    if (grepl("^[A-Z]{2,}(?:_[A-Z]{2,})*_[0-9]$", sn)) {
      corrected <- sub("_(\\d)$", "_0\\1", sn)
      message("🛠 Fixed Site Name: ", sn, " → ", corrected)
      return(corrected)
    }
    sn
  }, character(1))
  
  # --- Site Name validation ---
  valid_regex <- paste0(
    "^(",
    "[A-Z]{4}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{4}_[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{4}[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    ")$"
  )
  for (i in seq_len(nrow(deployment))) {
    sn <- trimws(deployment$`Site Name`[i])
    if (!is.na(sn) && sn != "" && !grepl(valid_regex, sn)) {
      issues <- c(issues, paste0(
        "Invalid Site_Name: ", sn,
        " → must follow an allowed format such as 'PARK_##', 'UNIT_##', 'PARK_UNIT_##', or 'PARKUNIT_##'."
      ))
    }
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
      if (!is.na(val) && val != "") {
        if (is.na(as.Date(val, format = "%m/%d/%Y"))) {
          issues <- c(issues, paste("❌ Bad date in", col, "row", i, ":", val, " — should be mm/dd/yyyy"))
        }
      }
    }
    
    # Camera Malfunction Date: only check once images are available,
    # and only for malfunctioning sites that actually have images.
    cam_func <- deployment$`Camera Functioning`[i]
    if (!is.na(cam_func) && tolower(cam_func) == "no") {
      val <- deployment$`Camera Malfunction Date`[i]
      cam_id <- deployment$`Site Name`[i]
      if (!is.null(images) && nrow(images) > 0 && cam_id %in% images$`Site Name`) {
        if (is.na(val) || val == "") {
          issues <- c(issues, paste("❌ Camera Malfunction Date missing in row", i,
                                    " — required because Camera Functioning = No"))
        } else if (is.na(as.Date(val, format = "%m/%d/%Y"))) {
          issues <- c(issues, paste("❌ Bad date in Camera Malfunction Date row", i, ":", val,
                                    " — should be mm/dd/yyyy"))
        }
      }
    }
    
    # --- Times ---
    time_cols <- c("Start Time", "End Time")
    for (col in time_cols) {
      val <- deployment[[col]][i]
      if (!is.na(val) && val != "") {
        if (!grepl("^(?:[01]?[0-9]|2[0-3]):[0-5][0-9](:[0-5][0-9])?$", val)) {
          issues <- c(issues, paste("❌ Bad time in", col, "row", i, ":", val, " — should be HH:MM or HH:MM:SS 24h"))
        }
      }
    }
    
    # ---- Numeric checks with suppression ----
    dd_val <- deployment$`Detection Distance`[i]
    if (!is.na(dd_val) && dd_val != "" && suppressWarnings(is.na(as.numeric(dd_val)))) {
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
    
    # --- Latitude check (must exist and be non-zero) ---
    lat_val <- suppressWarnings(as.numeric(deployment$Latitude[i]))
    if (is.na(lat_val) || lat_val == 0) {
      issues <- c(issues, paste("❌ Latitude missing or invalid (0) in row", i))
    }
    
    # --- Longitude check (must exist and be non-zero; auto-fix sign if positive) ---
    long_val <- suppressWarnings(as.numeric(deployment$Longitude[i]))
    if (is.na(long_val) || long_val == 0) {
      issues <- c(issues, paste("❌ Longitude missing or invalid (0) in row", i))
    } else if (long_val > 0) {
      # Auto-fix longitude sign
      fixed_val <- -abs(long_val)
      deployment$Longitude[i] <- fixed_val
      message("🛠 Fixed Longitude in row ", i, ": ", long_val, " → ", fixed_val)
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
    if (!is.na(val) && val != "") {
      val_lower <- tolower(val)
      if (!(val_lower %in% c("yes","no"))) {
        issues <- c(issues, paste("❌ Invalid Camera Functioning value in row", i, ":", val, " — must be Yes or No"))
      }
    }
    
  } # end row loop
  
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

check_images <- function(images, deployments, survey_year = NULL) {
  # survey_year is retained for backward compatibility with older callers.
  # The current app relies on parsed image timestamps instead.
  issues <- c()
  fixes  <- c()
  
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
  
  # ---- Fix Site Names missing leading zero ----
  fix_idx <- which(grepl("^[A-Z]+_[0-9]$", images$`Site Name`))
  if (length(fix_idx) > 0) {
    corrected_names <- sub("^(.*_)([0-9])$", "\\10\\2", images$`Site Name`[fix_idx])
    fixes <- paste(images$`Site Name`[fix_idx], "→", corrected_names)
    images$`Site Name`[fix_idx] <- corrected_names
    message("🛠 Fixed Site Names:\n", paste(unique(fixes), collapse = "\n"))
  }
  
  # ---- Site Name validation ----
  valid_regex <- paste0(
    "^(",
    "[A-Z]{4}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{4}_[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    "|",
    "[A-Z]{4}[A-Z]{2,}_(0[1-9]|[1-9][0-9]|[1-9][0-9]{2})",
    ")$"
  )
  invalid_idx <- which(!grepl(valid_regex, images$`Site Name`))
  if (length(invalid_idx) > 0) {
    site_counts <- table(images$`Site Name`[invalid_idx])
    for (sn in names(site_counts)) {
      issues <- c(issues, paste0(
        "❌ Invalid Site_Name: ", sn,
        " — found ", site_counts[[sn]], " occurrence(s); must follow an allowed site format such as 'PARK_##', 'UNIT_##', 'PARK_UNIT_##', or 'PARKUNIT_##'."
      ))
    }
  }
  
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
    
    # Auto-fix positive longitudes
    fix_lon <- which(!is.na(lon_num) & lon_num > 0 & lon_num <= 180)
    if (length(fix_lon) > 0) {
      images$Longitude[fix_lon] <- -abs(lon_num[fix_lon])
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
  
  # ---- Cluster ID numeric check ----
  cluster_num <- suppressWarnings(as.numeric(images$`Cluster ID`))
  bad_cluster <- which(!is.na(images$`Cluster ID`) & images$`Cluster ID` != "" & is.na(cluster_num))
  if (length(bad_cluster) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_cluster), " image(s) have non-numeric Cluster ID"))
  }
  
  # ---- Cross-check Site Names with deployments ----
  bad_site_match <- which(!images$`Site Name` %in% deployments$`Site Name`)
  if (length(bad_site_match) > 0) {
    issues <- c(issues, paste0("❌ ", length(bad_site_match), " image(s) have Site Names not found in deployments"))
  }
  
  # ---- Timestamp parsing ----
  ts_parsed <- parse_timestamp_robust(images$Timestamp)
  ts_missing_input <- !is.na(images$Timestamp) & trimws(as.character(images$Timestamp)) != ""
  bad_ts <- which(ts_missing_input & is.na(ts_parsed))
  if (length(bad_ts) > 0) {
    site_counts <- table(images$`Site Name`[bad_ts])
    for (sn in names(site_counts)) {
      issues <- c(issues, paste0("❌ Bad Timestamp format at site ", sn,
                                 " — ", site_counts[[sn]], " occurrence(s); use a recognizable date-time like yyyy-mm-dd HH:MM:SS or mm/dd/yyyy HH:MM."))
    }
  }
  images$Timestamp <- ts_parsed
  
  # ---- Timestamp vs deployment window check ----
  dt_deploy <- data.table::data.table(deployments)
  dt_deploy[, dep_start := as.Date(`Start Date`, format = "%m/%d/%Y")]
  dt_deploy[, dep_end   := as.Date(`End Date`,   format = "%m/%d/%Y")]
  
  dt_images <- data.table::data.table(images)
  dt_images[, ts := ts_parsed]
  dt_images[, img_date := as.Date(ts)]
  
  merged_dates <- merge(
    dt_images,
    dt_deploy[, .(`Site Name`, dep_start, dep_end)],
    by = "Site Name",
    all.x = TRUE
  )
  
  buffer_days <- 0L
  
  out_of_window <- merged_dates[
    !is.na(img_date) &
      (
        img_date < (dep_start - buffer_days) |
          img_date > (dep_end + buffer_days)
      )
  ]
  
  if (nrow(out_of_window) > 0) {
    first_img_out <- out_of_window[, .(first_img = min(img_date)), by = `Site Name`]
    
    for (i in seq_len(nrow(first_img_out))) {
      site <- first_img_out$`Site Name`[i]
      start_date <- dt_deploy$dep_start[dt_deploy$`Site Name` == site]
      first_img <- first_img_out$first_img[i]
      
      issues <- c(issues, paste0(
        "⚠️ Timestamp outside deployment window at site ", site, ": ",
        "Start Date = ", start_date, ", First image = ", first_img, ". ",
        "Check the tagging/export software and adjust the first photo timestamp if needed."
      ))
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
  
  # Parse dates with multiple format attempts (including 2-digit years)
  parse_date_safe <- function(x) {
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
    
    if (inherits(parsed, "POSIXct")) parsed <- as.Date(parsed)
    return(parsed)
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
# TRIM IMAGES TO A USER-DEFINED NUMBER OF DAYS
# -------------------------------------------------------------------

trim_images_to_days <- function(images, max_days = 56) {
  max_days <- as.integer(max_days[[1]])
  if (is.na(max_days) || max_days < 1) {
    stop("max_days must be a positive integer.", call. = FALSE)
  }

  # Ensure Timestamp is POSIXct
  if (!inherits(images$Timestamp, "POSIXt")) {
    images$Timestamp <- parse_timestamp_robust(images$Timestamp)
  }
  
  # Determine first timestamp per site
  first_times <- images %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarize(first_date = min(Timestamp, na.rm = TRUE), .groups = "drop")
  
  # Join to images to calculate days from first image
  images_check <- images %>%
    dplyr::left_join(first_times, by = "Site Name") %>%
    dplyr::mutate(days_from_start = as.numeric(difftime(Timestamp, first_date,
                                                        units = "days")))
  
  # Check if any rows exceed the requested deployment length
  rows_to_remove <- images_check %>% dplyr::filter(days_from_start > max_days)
  
  if (nrow(rows_to_remove) == 0) {
    message("✅ No images taken after ", max_days, " days. Dataset left unchanged.")
    return(images)
  } else {
    # Count rows removed per site
    removal_counts <- rows_to_remove %>%
      dplyr::group_by(`Site Name`) %>%
      dplyr::summarize(removed = n(), .groups = "drop")
    
    message("⚠️ Images beyond ", max_days, " days were removed:")
    for (i in seq_len(nrow(removal_counts))) {
      message(paste0("  Site ", removal_counts$`Site Name`[i], ": ",
                     removal_counts$removed[i], " row(s) removed"))
    }
    
    # Keep only rows within the requested deployment length
    images_trimmed <- images_check %>%
      dplyr::filter(days_from_start <= max_days) %>%
      dplyr::select(-first_date, -days_from_start)
    
    return(images_trimmed)
  }
}

trim_images_56days <- function(images) {
  trim_images_to_days(images, max_days = 56)
}

# -------------------------------------------------------------------
# SUMMARY HELPERS (deployment, images, deer)
# -------------------------------------------------------------------

summarize_deployments <- function(deployments) {
  deployments %>%
    dplyr::mutate(
      Start_Date = as.Date(`Start Date`, format = "%m/%d/%Y"),
      End_Date   = as.Date(`End Date`,   format = "%m/%d/%Y"),
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

standardize_deer_species <- function(images) {
  if (any(tolower(images$Species) %in% c("white-tailed deer", "white tailed deer"))) {
    images <- images %>%
      dplyr::mutate(
        Species = dplyr::if_else(
          tolower(Species) %in% c("white-tailed deer", "white tailed deer"),
          "Deer",
          Species
        )
      )
    message("🦌 Standardized all 'white-tailed deer' classifications to 'Deer' for consistency.")
  } else {
    message("✅ No 'white-tailed deer' entries found — all deer species names already standardized.")
  }
  images
}

deer_summary_per_site <- function(images) {
  species_summary_per_site(images, "Deer")
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

deer_counts_per_camera <- function(images) {
  species_counts_per_camera(images, "Deer")
}

first_finite_or_na <- function(x) {
  vals <- suppressWarnings(as.numeric(x))
  vals <- vals[is.finite(vals)]
  if (length(vals)) vals[1] else NA_real_
}

species_counts_per_camera <- function(images, species_name, deployments = NULL) {
  counts <- filter_species_rows(images, species_name) %>%
    dplyr::group_by(`Site Name`) %>%
    dplyr::summarise(
      total_detections = sum(as.numeric(`Sighting Count`), na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(deployments) &&
      all(c("Site Name", "Latitude", "Longitude") %in% names(deployments))) {
    coords <- deployments %>%
      dplyr::group_by(`Site Name`) %>%
      dplyr::summarise(
        Latitude = first_finite_or_na(Latitude),
        Longitude = first_finite_or_na(Longitude),
        .groups = "drop"
      )
    counts <- counts %>%
      dplyr::left_join(coords, by = "Site Name")
  } else if (all(c("Latitude", "Longitude") %in% names(images))) {
    coords <- images %>%
      dplyr::group_by(`Site Name`) %>%
      dplyr::summarise(
        Latitude = first_finite_or_na(Latitude),
        Longitude = first_finite_or_na(Longitude),
        .groups = "drop"
      )
    counts <- counts %>%
      dplyr::left_join(coords, by = "Site Name")
  }

  counts
}

deer_daily_detections <- function(images) {
  species_daily_detections(images, "Deer")
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
