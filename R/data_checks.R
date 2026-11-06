# R/data_checks.R
# Cleaning, quality-control (QC), analysis-window, and summary helpers for
# uploaded camera data.
#
# In the QC messages, "❌" marks a problem that blocks the upload, "⚠️" marks
# a warning that the user should review but that doesn't block it, and "ℹ️"
# marks a note about data the models will leave out.

# Range of Camera Detection Angles (degrees, inclusive) the app accepts: the
# angles TTE's average path length across the viewshed
# (tte_path_length_multiplier()) was fitted to. Uploaded angles outside this
# range get a warning, and the models use the default angle from Model
# settings for those cameras.
camera_angle_min_deg <- 5
camera_angle_max_deg <- 90

# -------------------------------------------------------------------
# Parsing helpers
# -------------------------------------------------------------------

# Parse image timestamps written as month/day/year or year-month-day, with an
# optional time (HH:MM or HH:MM:SS) and 2- or 4-digit years. Returns POSIXct;
# values that can't be parsed become NA.
parse_timestamp_robust <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXt") || inherits(x, "Date")) {
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
        "mdy HMS", "mdy HM", "mdy"
      ),
      tz = tz,
      quiet = TRUE
    )
  )

  as.POSIXct(parsed, tz = tz)
}

# Strict MM/DD/YYYY parser for deployment dates. as.Date(x, format = "%m/%d/%Y")
# on its own is too forgiving: "6/15/24" silently becomes the year 0024 and
# trailing text ("06/15/2024abc") is ignored. This requires the whole value to
# be m/d/yyyy (1-2 digit month and day, 4-digit year) and returns NA otherwise.
parse_mdY_strict <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(format(x, "%Y-%m-%d")))
  }
  x <- trimws(as.character(x))
  ok <- !is.na(x) & grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", x)
  out <- rep(as.Date(NA), length(x))
  out[ok] <- as.Date(x[ok], format = "%m/%d/%Y")
  out
}

# Rows of an uploaded table in which every column is blank (for example, rows
# of commas that a spreadsheet can leave at the end of a CSV file).
blank_rows <- function(df) {
  blank <- rep(TRUE, nrow(df))
  for (x in df) {
    blank <- blank & (is.na(x) | trimws(as.character(x)) == "")
  }
  blank
}

# Check-log note about the blank rows removed from an uploaded file.
blank_rows_note <- function(blank) {
  rows <- which(blank)
  paste0(
    "ℹ️ Removed ", length(rows), " blank row(s): row(s) ",
    paste(head(rows, 10), collapse = ", "),
    if (length(rows) > 10) ", ..." else "", "."
  )
}

# -------------------------------------------------------------------
# Pre-QC cleaning
# -------------------------------------------------------------------

clean_deployment_import <- function(deployments) {
  # Standardize column names: drop trailing dots, trim whitespace (including
  # non-breaking spaces), and collapse internal runs of whitespace.
  nm <- names(deployments)
  nm <- gsub("\\.+$", "", nm)
  nm <- gsub("[\\s\\p{Z}]+$", "", nm, perl = TRUE)
  nm <- gsub("^[\\s\\p{Z}]+", "", nm, perl = TRUE)
  nm <- gsub("[\\s\\p{Z}]+", " ", nm, perl = TRUE)
  names(deployments) <- trimws(nm)

  # Trim leading/trailing whitespace from all character columns.
  deployments[] <- lapply(deployments, function(x) {
    if (is.character(x)) {
      gsub("^[\\s\\p{Z}]+|[\\s\\p{Z}]+$", "", x, perl = TRUE)
    } else {
      x
    }
  })

  # Camera Functioning: accept common variants and convert to Yes / No.
  if ("Camera Functioning" %in% names(deployments)) {
    x <- deployments$`Camera Functioning`
    if (is.logical(x)) {
      deployments$`Camera Functioning` <- ifelse(x, "Yes", "No")
    } else {
      xc <- trimws(as.character(x))
      low <- tolower(xc)
      yn <- xc
      yn[low %in% c("yes", "y", "true", "t", "1")] <- "Yes"
      yn[low %in% c("no", "n", "false", "f", "0")] <- "No"
      deployments$`Camera Functioning` <- yn
    }
  }

  deployments
}

clean_images_import <- function(images) {
  if (is.null(images) || !is.data.frame(images)) {
    stop("The images file could not be read as a table.", call. = FALSE)
  }

  # Remove blank rows, keeping each row's position in the file so problems
  # found below (including after multi-species rows are split) are reported
  # by file row.
  blank <- blank_rows(images)
  images$.row <- seq_len(nrow(images))
  if (any(blank)) {
    message(blank_rows_note(blank))
    images <- images[!blank, , drop = FALSE]
  }

  if (nrow(images) == 0) {
    warning("⚠️ The images file has no rows.")
    images$.row <- NULL
    return(images)
  }
  if (length(names(images)) == 0) {
    stop("The images file has no column names.", call. = FALSE)
  }

  # Standardize column names.
  names(images) <- trimws(gsub("\\.+$", "", names(images)))

  # Parse Timestamp (dropping an explicit " UTC" suffix). The column is only
  # converted when every non-blank value parses; otherwise it is left as text
  # so check_images() can report the unparseable values as a bad format
  # instead of as missing.
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
  # for 1 deer and one row for 4 squirrels; all other columns are repeated).
  if (all(c("Species", "Sighting Count") %in% names(images))) {
    images$`Sighting Count` <- as.character(images$`Sighting Count`)

    # separate_rows() fails with an unhelpful message if a row's
    # Species and Sighting Count split into different numbers of pieces, so
    # check rows that use the pipe-delimited format first.
    has_pipe <- grepl("|", images$Species, fixed = TRUE) |
      grepl("|", images$`Sighting Count`, fixed = TRUE)
    has_pipe[is.na(has_pipe)] <- FALSE
    if (any(has_pipe)) {
      species_n <- lengths(strsplit(images$Species[has_pipe], "|", fixed = TRUE))
      count_n <- lengths(strsplit(images$`Sighting Count`[has_pipe], "|", fixed = TRUE))
      mismatch_rows <- images$.row[has_pipe][species_n != count_n]
      if (length(mismatch_rows) > 0) {
        stop(
          "❌ Species and Sighting Count have a different number of ",
          "pipe-delimited values in row(s) ", paste(mismatch_rows, collapse = ", "),
          ". Each pipe-delimited Species entry needs a matching Sighting Count ",
          "entry (for example, Species = 'Deer|Squirrel' needs Sighting Count = '1|4').",
          call. = FALSE
        )
      }
    }

    images <- separate_rows(images, Species, `Sighting Count`, sep = "\\|")

    # Every row with a Species needs a numeric Sighting Count and a Cluster ID,
    # since each event's count is the largest Sighting Count among its images.
    # Rows without a Species (for example, triggers with no animal) aren't
    # counted for any species, so they're left alone here; check_images()
    # reports them as a warning.
    has_species <- !is.na(images$Species) & trimws(images$Species) != ""
    count_text <- trimws(images$`Sighting Count`)
    count_blank <- is.na(count_text) | count_text == ""
    count_bad <- !count_blank & is.na(suppressWarnings(as.numeric(count_text)))
    cluster_blank <- if ("Cluster ID" %in% names(images)) {
      is.na(images$`Cluster ID`) | trimws(images$`Cluster ID`) == ""
    } else {
      rep(FALSE, nrow(images))
    }
    rows_with <- function(flag) {
      rows <- unique(images$.row[has_species & flag])
      paste0(
        paste(head(rows, 10), collapse = ", "),
        if (length(rows) > 10) ", ..." else ""
      )
    }
    problems <- c(
      if (any(has_species & count_blank)) {
        paste0("Sighting Count is missing in row(s) ", rows_with(count_blank), ".")
      },
      if (any(has_species & count_bad)) {
        paste0(
          "Sighting Count isn't a number in row(s) ", rows_with(count_bad), " (",
          paste(head(unique(count_text[has_species & count_bad]), 5), collapse = ", "),
          ")."
        )
      },
      if (any(has_species & cluster_blank)) {
        paste0("Cluster ID is missing in row(s) ", rows_with(cluster_blank), ".")
      }
    )
    if (length(problems) > 0) {
      stop(
        "❌ ", paste(problems, collapse = " "),
        " Every row with a Species needs a numeric Sighting Count (or pipe-delimited ",
        "numbers for multi-species rows, for example '1|4') and a Cluster ID. ",
        "Fix these rows, then re-upload.",
        call. = FALSE
      )
    }

    images$`Sighting Count` <- as.numeric(count_text)
  }

  images$.row <- NULL
  images
}

# -------------------------------------------------------------------
# QC: deployments
# -------------------------------------------------------------------

check_deployments <- function(deployment, hemisphere = "Western",
                              lat_hemisphere = "Northern") {
  issues <- character()
  hemisphere <- if (identical(hemisphere, "Eastern")) "Eastern" else "Western"
  lat_hemisphere <- if (identical(lat_hemisphere, "Southern")) "Southern" else "Northern"

  # Problems that block the upload are collected here and reported together
  # after the per-row checks: a missing Site Name; a missing, badly formatted,
  # or reversed Start Date/End Date; a missing or zero Latitude/Longitude; any
  # Camera Malfunction Date problem; and a missing, non-numeric, or
  # non-positive Detection Distance. Every model depends on the dates and
  # coordinates, which set each camera's analysis window and location; REM
  # and TTE also depend on Detection Distance.
  blocking_issues <- character()

  # ---- Column presence (hard stop) ----
  required_cols <- c(
    "Site Name",
    "Start Date", "Start Time", "End Date", "End Time",
    "Latitude", "Longitude", "Camera Functioning",
    "Camera Malfunction Date", "Detection Distance"
  )

  missing_cols <- setdiff(required_cols, names(deployment))
  if (length(missing_cols) > 0) {
    stop(
      "❌ Deployment file is missing required column(s): ",
      paste(missing_cols, collapse = ", "), ". ",
      "If these columns are truly missing from your dataset (for example, no ",
      "Start Time or End Time was ever recorded), add the column names to your ",
      "CSV file and leave the cells blank. If you believe these columns exist but ",
      "are labeled slightly differently (e.g. 'Start_Date' instead of 'Start Date'), ",
      "rename them to match exactly. After fixing the file, re-upload it.",
      call. = FALSE
    )
  }

  # ---- Site Name uniqueness (hard stop) ----
  # A duplicate Site Name would collide downstream (camera counts, coordinate
  # lookups, the camera array used by USCR), so it stops the upload right away.
  site_names <- trimws(deployment$`Site Name`)
  dup_sites <- unique(site_names[duplicated(site_names) & !is.na(site_names) & site_names != ""])
  if (length(dup_sites) > 0) {
    stop(
      "❌ Deployment file has duplicate Site Name value(s): ",
      paste(dup_sites, collapse = ", "), ". ",
      "Each row must have a unique Site Name. If two cameras were deployed at the ",
      "same physical site at different times, give them distinct Site Names (for ",
      "example, add a suffix). If this is a duplicate row from a data-entry mistake, ",
      "remove the extra row(s). After fixing the file, re-upload it.",
      call. = FALSE
    )
  }

  # ---- Blank rows (removed) ----
  # Rows in which every column is blank are skipped below and removed from
  # the checked data. Row numbers in the messages are rows of the file.
  blank <- blank_rows(deployment)
  if (any(blank)) message(blank_rows_note(blank))

  # ---- Per-row checks ----
  for (i in seq_len(nrow(deployment))) {
    if (blank[i]) next

    # --- Site Name (blocking) ---
    # Any naming scheme is accepted, but a row with deployment data and no
    # Site Name can't be linked to its images.
    sn <- deployment$`Site Name`[i]
    if (is.na(sn) || sn == "") {
      msg <- paste("❌ Site Name missing in row", i)
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    }

    # --- Dates (blocking) ---
    for (col in c("Start Date", "End Date")) {
      val <- deployment[[col]][i]
      if (is.na(val) || val == "") {
        msg <- paste("❌", col, "missing in row", i)
        issues <- c(issues, msg)
        blocking_issues <- c(blocking_issues, msg)
      } else if (is.na(parse_mdY_strict(val))) {
        msg <- paste0("❌ Bad date in ", col, " row ", i, ": ", val,
                      " (should be MM/DD/YYYY with a 4-digit year)")
        issues <- c(issues, msg)
        blocking_issues <- c(blocking_issues, msg)
      }
    }

    # End Date must not be before Start Date: a reversed range would give
    # that camera a collapsed or wrong analysis window.
    start_parsed <- parse_mdY_strict(deployment$`Start Date`[i])
    end_parsed <- parse_mdY_strict(deployment$`End Date`[i])
    if (!is.na(start_parsed) && !is.na(end_parsed) && end_parsed < start_parsed) {
      msg <- paste0("❌ End Date is before Start Date in row ", i, ": ",
                    deployment$`Start Date`[i], " → ", deployment$`End Date`[i])
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    }

    # --- Camera Malfunction Date (blocking) ---
    # Required whenever Camera Functioning = No. It is the FIRST day whose
    # data should be discarded: the camera's analysis window ends at midnight
    # at the start of it (see analysis_windows()). Every problem with it
    # blocks the upload, because otherwise images taken after the malfunction
    # (e.g. a camera knocked to a bad angle that kept firing) would be
    # included.
    cam_func <- deployment$`Camera Functioning`[i]
    if (!is.na(cam_func) && tolower(cam_func) == "no") {
      val <- deployment$`Camera Malfunction Date`[i]
      msg <- NULL
      if (is.na(val) || val == "") {
        msg <- paste0("❌ Camera Malfunction Date missing in row ", i,
                      " (required because Camera Functioning = No)")
      } else if (is.na(parse_mdY_strict(val))) {
        msg <- paste0("❌ Bad date in Camera Malfunction Date row ", i, ": ", val,
                      " (should be MM/DD/YYYY with a 4-digit year)")
      } else {
        mal_parsed <- parse_mdY_strict(val)
        if (!is.na(start_parsed) && mal_parsed <= start_parsed) {
          msg <- paste0(
            "❌ Camera Malfunction Date in row ", i, " (", val,
            ") is on or before Start Date (", deployment$`Start Date`[i],
            "). The malfunction date is the first day to discard, so this camera ",
            "would have no usable days. Correct the date or remove the row."
          )
        } else if (!is.na(end_parsed) && mal_parsed > end_parsed) {
          msg <- paste0(
            "❌ Camera Malfunction Date in row ", i, " (", val,
            ") is after End Date (", deployment$`End Date`[i], "). If the camera worked ",
            "until it was picked up, set Camera Functioning = Yes and leave ",
            "Camera Malfunction Date blank."
          )
        }
      }
      if (!is.null(msg)) {
        issues <- c(issues, msg)
        blocking_issues <- c(blocking_issues, msg)
      }
    } else if (!is.na(cam_func) && tolower(cam_func) == "yes") {
      val <- deployment$`Camera Malfunction Date`[i]
      if (!is.na(val) && val != "") {
        issues <- c(issues, paste0(
          "⚠️ Camera Malfunction Date given in row ", i, " (", val,
          ") but Camera Functioning = Yes, so the date is ignored. If the camera ",
          "malfunctioned, set Camera Functioning = No."
        ))
      }
    }

    # --- Times (warning) ---
    # A missing or unreadable time shortens the analysis window (see
    # analysis_windows()) rather than blocking the upload.
    fallback <- c(
      "Start Time" = "the analysis window starts at midnight at the end of Start Date",
      "End Time" = "the analysis window ends at midnight at the start of End Date"
    )
    for (col in c("Start Time", "End Time")) {
      val <- deployment[[col]][i]
      if (is.na(val) || val == "") {
        issues <- c(issues, paste0("⚠️ ", col, " missing in row ", i, "; ", fallback[[col]]))
      } else if (is.na(parse_time_of_day(val))) {
        issues <- c(issues, paste0("⚠️ Bad time in ", col, " row ", i, ": ", val,
                                   " (should be HH:MM or HH:MM:SS, 24-hour); ",
                                   fallback[[col]]))
      }
    }

    # --- Detection Distance (blocking) ---
    # REM and TTE need a positive detection radius for every camera.
    dd_val <- deployment$`Detection Distance`[i]
    dd_num <- suppressWarnings(as.numeric(dd_val))
    msg <- if (is.na(dd_val) || dd_val == "") {
      paste("❌ Detection Distance missing in row", i)
    } else if (is.na(dd_num)) {
      paste0("❌ Detection Distance isn't a number in row ", i, ": ", dd_val)
    } else if (dd_num <= 0) {
      paste0("❌ Detection Distance isn't positive in row ", i, ": ", dd_val)
    }
    if (!is.null(msg)) {
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    }

    # --- Other numeric fields (warning) ---

    if ("Camera Detection Angle" %in% names(deployment)) {
      theta_val <- deployment$`Camera Detection Angle`[i]
      theta_num <- suppressWarnings(as.numeric(theta_val))
      if (is.na(theta_val) || theta_val == "") {
        issues <- c(issues, paste0(
          "⚠️ Camera Detection Angle missing in row ", i,
          "; the default angle from Model settings is used instead"
        ))
      } else if (is.na(theta_num)) {
        issues <- c(issues, paste0(
          "⚠️ Non-numeric value in Camera Detection Angle row ", i, ": ", theta_val,
          "; the default angle from Model settings is used instead"
        ))
      } else if (theta_num < camera_angle_min_deg || theta_num > camera_angle_max_deg) {
        issues <- c(issues, paste0(
          "⚠️ Camera Detection Angle out of range in row ", i, ": ", theta_val,
          " (must be from ", camera_angle_min_deg, " to ", camera_angle_max_deg,
          " degrees); the default angle from Model settings is used instead"
        ))
      }
    }

    if ("Camera Height" %in% names(deployment)) {
      ch_val <- deployment$`Camera Height`[i]
      if (!is.na(ch_val) && ch_val != "" && is.na(suppressWarnings(as.numeric(ch_val)))) {
        issues <- c(issues, paste0("⚠️ Non-numeric value in Camera Height row ", i, ": ", ch_val))
      }
    }

    # --- Latitude (blocking if missing or zero; sign fixed to match the
    # hemisphere selected in Model settings) ---
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

    # --- Longitude (same rules as Latitude) ---
    long_val <- suppressWarnings(as.numeric(deployment$Longitude[i]))
    if (is.na(long_val) || long_val == 0) {
      msg <- paste("❌ Longitude missing or invalid (0) in row", i)
      issues <- c(issues, msg)
      blocking_issues <- c(blocking_issues, msg)
    } else {
      fixed_lon <- if (identical(hemisphere, "Eastern")) abs(long_val) else -abs(long_val)
      if (fixed_lon != long_val) {
        deployment$Longitude[i] <- fixed_lon
        message("🛠 Fixed Longitude in row ", i, ": ", long_val, " → ", fixed_lon)
      }
    }

    # --- Camera Orientation (warning) ---
    if ("Camera Orientation" %in% names(deployment)) {
      val <- deployment$`Camera Orientation`[i]
      valid_cardinals <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
      if (!is.na(val) && val != "") {
        val_numeric <- suppressWarnings(as.numeric(val))
        valid <- toupper(val) %in% valid_cardinals ||
          (!is.na(val_numeric) && val_numeric >= 0 && val_numeric <= 359)
        if (!valid) {
          issues <- c(issues, paste0(
            "⚠️ Invalid Camera Orientation in row ", i, ": ", val,
            " (must be N, NE, E, SE, S, SW, W, NW, or 0-359 degrees)"
          ))
        }
      }
    }

    # --- Camera Functioning (warning) ---
    val <- deployment$`Camera Functioning`[i]
    if (is.na(val) || val == "") {
      issues <- c(issues, paste("⚠️ Camera Functioning missing in row", i))
    } else if (!tolower(val) %in% c("yes", "no")) {
      issues <- c(issues, paste0("⚠️ Invalid Camera Functioning value in row ", i, ": ", val,
                                 " (must be Yes or No)"))
    }
  }

  # ---- Duplicate coordinates (warning) ----
  # Two different Site Names sharing the exact same Latitude AND Longitude
  # usually means a copy-paste error. Latitude or Longitude values may repeat
  # individually (e.g. cameras on a north-south transect); only a duplicated
  # pair is flagged. Runs after the per-row loop so it uses the sign-corrected
  # coordinates.
  lat_num_all <- suppressWarnings(as.numeric(deployment$Latitude))
  lon_num_all <- suppressWarnings(as.numeric(deployment$Longitude))
  has_coords <- !is.na(lat_num_all) & !is.na(lon_num_all)
  coord_key <- paste(lat_num_all, lon_num_all)
  dup_keys <- unique(coord_key[has_coords][duplicated(coord_key[has_coords])])
  for (k in dup_keys) {
    dup_rows <- which(has_coords & coord_key == k)
    issues <- c(issues, paste0(
      "⚠️ Duplicate coordinates in ", length(dup_rows), " row(s) (Site Name: ",
      paste(deployment$`Site Name`[dup_rows], collapse = ", "),
      "): same Latitude/Longitude; verify these are meant to be different camera locations."
    ))
  }

  # ---- Analysis window (blocking) ----
  # Every camera needs some time between the start and end of its analysis
  # window. Rows with dates that couldn't be read are already blocked above.
  windows <- analysis_windows(deployment)
  for (i in which(!is.na(windows$days) & windows$days <= 0)) {
    msg <- paste0(
      "❌ No deployment time in row ", i, ": the analysis window would start at ",
      format(windows$window_start[i], "%m/%d/%Y %H:%M"), " and end at ",
      format(windows$window_end[i], "%m/%d/%Y %H:%M"),
      ". Check Start/End Date and Time (and Camera Malfunction Date)."
    )
    issues <- c(issues, msg)
    blocking_issues <- c(blocking_issues, msg)
  }

  # ---- Hard stop for blocking problems ----
  if (length(blocking_issues) > 0) {
    stop(
      "Deployment file has missing or invalid value(s) that block upload in ",
      length(blocking_issues), " place(s):\n",
      paste0("  - ", blocking_issues, collapse = "\n"),
      "\nEvery row needs a Site Name. Start Date and End Date must be present, in ",
      "MM/DD/YYYY format, and End Date must be on or after Start Date for every row; ",
      "with Start Time and End Time, they set each camera's analysis window, which ",
      "must have some time in it. Latitude and Longitude must be ",
      "present and non-zero for every row; they're used to project camera locations ",
      "for every model. When Camera Functioning = No, Camera Malfunction Date must be ",
      "an MM/DD/YYYY date after Start Date and no later than End Date. ",
      "Detection Distance must be a positive number (meters) for every row; REM and ",
      "TTE depend on it. ",
      "After fixing the file, re-upload it.",
      call. = FALSE
    )
  }

  # ---- Output ----
  if (length(issues) == 0) {
    message("✅ Deployments file is formatted correctly!")
  } else {
    warning("⚠️ Issues found in deployments file:\n", paste(issues, collapse = "\n"))
  }

  deployment <- deployment[!blank, , drop = FALSE]
  deployment$`Detection Distance` <- suppressWarnings(as.numeric(deployment$`Detection Distance`))
  deployment
}

# -------------------------------------------------------------------
# QC: images
# -------------------------------------------------------------------

check_images <- function(images, deployments) {
  issues <- character()

  # ---- Required columns (hard stop) ----
  required_cols <- c("Site Name", "Timestamp", "Species", "Sighting Count", "Cluster ID")
  missing_cols <- setdiff(required_cols, names(images))
  if (length(missing_cols) > 0) {
    stop(
      "❌ Images file is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ". Please fix the images file and re-upload it.",
      call. = FALSE
    )
  }

  # ---- Trim whitespace from character columns ----
  char_cols <- names(images)[vapply(images, is.character, logical(1))]
  images[char_cols] <- lapply(images[char_cols], trimws)

  # ---- Species (warning) ----
  # Rows without a Species (for example, triggers with no animal) aren't
  # counted for any species. Sighting Count and Cluster ID are required on
  # every other row (checked in clean_images_import()).
  missing_species <- which(is.na(images$Species) | images$Species == "")
  if (length(missing_species) > 0) {
    issues <- c(issues, paste0(
      "⚠️ Species is blank in ", length(missing_species),
      " image(s); these aren't counted for any species"
    ))
  }

  # ---- Cluster ID at more than one site (warning) ----
  # Cluster ID may be any format, and multiple rows legitimately share one
  # (e.g. burst photos, or a multi-species image split into one row per
  # species). It should, however, identify detections at a single camera.
  # detection_events() groups by Site Name and Cluster ID together, so such
  # clusters are kept separate per camera, but this usually signals a problem
  # with how Cluster IDs were assigned upstream.
  cluster_site_counts <- images |>
    filter(!is.na(`Cluster ID`) & `Cluster ID` != "") |>
    group_by(`Cluster ID`) |>
    summarize(n_sites = n_distinct(`Site Name`), .groups = "drop")
  bad_clusters <- cluster_site_counts$`Cluster ID`[cluster_site_counts$n_sites > 1]
  if (length(bad_clusters) > 0) {
    issues <- c(issues, paste0(
      "⚠️ ", length(bad_clusters), " Cluster ID value(s) appear at more than one Site Name: ",
      paste(head(bad_clusters, 10), collapse = ", "),
      if (length(bad_clusters) > 10) ", ..." else ""
    ))
  }

  # ---- Site Names must match the deployment file (hard stop) ----
  # detection_events() matches images to cameras' analysis windows by Site
  # Name, so images whose Site Name has no deployment row would be left out
  # of the counts without notice. Blank Site Names are reported separately.
  missing_site_idx <- which(is.na(images$`Site Name`) | images$`Site Name` == "")
  bad_site_match <- setdiff(
    which(!images$`Site Name` %in% deployments$`Site Name`),
    missing_site_idx
  )
  if (length(missing_site_idx) > 0 || length(bad_site_match) > 0) {
    bad_site_names <- unique(images$`Site Name`[bad_site_match])
    stop(
      "❌ ",
      if (length(missing_site_idx) > 0) {
        paste0(
          "Images file has ", length(missing_site_idx), " row(s) with a blank Site Name ",
          "(first few rows: ", paste(head(missing_site_idx, 10), collapse = ", "), "). "
        )
      } else {
        ""
      },
      if (length(bad_site_match) > 0) {
        paste0(
          "Images file has ", length(bad_site_match), " row(s) with Site Name value(s) not found ",
          "in the deployment file: ", paste(head(bad_site_names, 10), collapse = ", "),
          if (length(bad_site_names) > 10) ", ..." else "", ". "
        )
      } else {
        ""
      },
      "Every Site Name in the images file must exactly match a Site Name in the ",
      "deployment file. Fix the blank or mismatched Site Name(s) (check for typos ",
      "or a missing deployment row), then re-upload.",
      call. = FALSE
    )
  }

  # ---- Timestamps (hard stop) ----
  # An image without a readable Timestamp can't be placed in any camera's
  # analysis window, so it would be dropped from the counts. Rows are
  # reported by site rather than row number, because multi-species rows have
  # already been split in clean_images_import().
  ts_text <- trimws(as.character(images$Timestamp))
  ts_parsed <- parse_timestamp_robust(images$Timestamp)
  has_ts <- !is.na(ts_text) & !ts_text %in% c("", "NA", "NaN")
  missing_ts_idx <- which(!has_ts)
  bad_ts_idx <- which(has_ts & is.na(ts_parsed))
  if (length(missing_ts_idx) > 0 || length(bad_ts_idx) > 0) {
    sites_of <- function(idx) {
      sites <- unique(images$`Site Name`[idx])
      paste0(
        paste(head(sites, 10), collapse = ", "),
        if (length(sites) > 10) ", ..." else ""
      )
    }
    stop(
      "❌ ",
      if (length(missing_ts_idx) > 0) {
        paste0(
          "Timestamp is missing in ", length(missing_ts_idx), " image(s) (Site Name: ",
          sites_of(missing_ts_idx), "). "
        )
      } else {
        ""
      },
      if (length(bad_ts_idx) > 0) {
        bad_values <- unique(ts_text[bad_ts_idx])
        paste0(
          "Timestamp can't be read in ", length(bad_ts_idx), " image(s) (Site Name: ",
          sites_of(bad_ts_idx), "; value(s): ",
          paste(head(bad_values, 5), collapse = ", "),
          if (length(bad_values) > 5) ", ..." else "", "). "
        )
      } else {
        ""
      },
      "Every image needs a Timestamp, written as month/day/year or year-month-day ",
      "with an optional 24-hour time (for example, MM/DD/YYYY HH:MM:SS or ",
      "YYYY-MM-DD HH:MM), so it can be placed in its camera's analysis window. ",
      "Fix or remove these rows, then re-upload.",
      call. = FALSE
    )
  }
  images$Timestamp <- ts_parsed

  # ---- Image dates vs deployment dates (warning) ----
  # Deployment Start/End Times are approximate (they're written down when
  # someone visits the camera), so an image a few hours past End Time is
  # normal. This check looks for timestamps that are genuinely wrong, usually
  # from a camera whose date or clock was set incorrectly.
  #
  # So the comparison uses calendar DATES only and ignores Start/End Time. An
  # image is flagged only if its date is more than a day (`grace_days`) before
  # Start Date or after the date the camera was picked up.
  #
  # Malfunctioning cameras often keep firing until they're picked up, so
  # images between the Camera Malfunction Date and End Date are expected and
  # aren't a clock problem. They're reported as information, since the models
  # leave them out (the analysis window ends at midnight at the start of the
  # malfunction date; see analysis_windows()).
  grace_days <- 1

  malfunction_date <- parse_mdY_strict(deployments$`Camera Malfunction Date`)
  malfunction_date[!tolower(deployments$`Camera Functioning`) %in% "no"] <- NA

  deployment_windows <- data.frame(
    `Site Name` = deployments$`Site Name`,
    dep_start = parse_mdY_strict(deployments$`Start Date`),
    dep_end = parse_mdY_strict(deployments$`End Date`),
    malfunction_date = malfunction_date,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Calendar date as shown on the timestamp (no time-zone shifting).
  img_dates <- images |>
    transmute(`Site Name`, img_date = as.Date(format(Timestamp, "%Y-%m-%d")))

  window_check <- img_dates |>
    inner_join(deployment_windows, by = "Site Name") |>
    group_by(`Site Name`) |>
    summarize(
      dep_start = first(dep_start),
      dep_end = first(dep_end),
      malfunction_date = first(malfunction_date),
      first_image = min(img_date),
      last_image = max(img_date),
      n_images = n(),
      n_before = sum(img_date < dep_start - grace_days, na.rm = TRUE),
      n_after = sum(img_date > dep_end + grace_days, na.rm = TRUE),
      # Images on or after the malfunction date but not past pickup (+ grace):
      # the models drop these. Zero for working cameras (malfunction_date NA).
      n_post_malfunction = sum(
        img_date >= malfunction_date & img_date <= dep_end + grace_days,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  fmt_d <- function(d) format(d, "%m/%d/%Y")

  for (i in seq_len(nrow(window_check))) {
    w <- window_check[i, ]

    if (w$n_before > 0) {
      issues <- c(issues, paste0(
        "⚠️ Site ", w$`Site Name`, ": ", w$n_before, " of ", w$n_images,
        " image(s) are dated more than a day before the deployment ",
        "Start Date (Start Date = ", fmt_d(w$dep_start), "; earliest image date = ",
        fmt_d(w$first_image), "). Check that the camera's date/time was set correctly."
      ))
    }

    if (w$n_after > 0) {
      issues <- c(issues, paste0(
        "⚠️ Site ", w$`Site Name`, ": ", w$n_after, " of ", w$n_images,
        " image(s) are dated more than a day after the camera was ",
        "picked up (End Date = ", fmt_d(w$dep_end), "; latest image date = ",
        fmt_d(w$last_image), "). Check that the camera's date/time was set correctly."
      ))
    }

    if (w$n_post_malfunction > 0) {
      message(
        "ℹ️ Site ", w$`Site Name`, ": ", w$n_post_malfunction,
        " image(s) are dated on or after the Camera Malfunction Date (",
        fmt_d(w$malfunction_date), ") and will be excluded from the analysis."
      )
    }
  }

  # ---- Detection events outside the analysis window (information) ----
  # The models count a detection event only if its first image falls inside
  # the camera's analysis window (see analysis_windows()). Events that start
  # within a day of the deployment dates but outside the window would be
  # dropped silently, and can point to a wrong Start/End Time. Events already
  # flagged by the date check above, and events after a Camera Malfunction
  # Date, are reported there instead. The optional trim isn't applied here.
  windows <- analysis_windows(deployments)
  outside <- images |>
    filter(!is.na(Species), Species != "") |>
    group_by(`Site Name`, `Cluster ID`) |>
    summarize(start = min(Timestamp), .groups = "drop") |>
    inner_join(windows, by = "Site Name") |>
    inner_join(deployment_windows, by = "Site Name") |>
    mutate(start_date = as.Date(format(start, "%Y-%m-%d"))) |>
    filter(
      start_date >= dep_start - grace_days,
      start_date <= dep_end + grace_days
    ) |>
    group_by(`Site Name`, window_start, window_end) |>
    summarize(
      n_early = sum(start < window_start),
      n_late = sum(start > window_end & !malfunction),
      .groups = "drop"
    )

  fmt_dt <- function(x) format(x, "%m/%d/%Y %H:%M")
  for (i in seq_len(nrow(outside))) {
    w <- outside[i, ]
    if (w$n_early > 0) {
      message(
        "ℹ️ Site ", w$`Site Name`, ": ", w$n_early, " detection event(s) start before ",
        "the analysis window (", fmt_dt(w$window_start), ") and won't be counted. ",
        "If the camera was already running, check Start Date and Start Time."
      )
    }
    if (w$n_late > 0) {
      message(
        "ℹ️ Site ", w$`Site Name`, ": ", w$n_late, " detection event(s) start after ",
        "the analysis window (", fmt_dt(w$window_end), ") and won't be counted. ",
        "If the camera was still running, check End Date and End Time."
      )
    }
  }

  # ---- Output ----
  if (length(issues) == 0) {
    message("✅ Images file is formatted correctly and all Site Names match deployments!")
  } else {
    warning("⚠️ Issues found in images file:\n", paste(issues, collapse = "\n"))
  }

  images
}

# -------------------------------------------------------------------
# Analysis windows (deployment date-times and optional trimming)
# -------------------------------------------------------------------

# Seconds after midnight for HH:MM or HH:MM:SS (24-hour) times; NA for blank
# or unreadable values.
parse_time_of_day <- function(x) {
  x <- trimws(as.character(x))
  m <- regmatches(x, regexec("^([01]?[0-9]|2[0-3]):([0-5][0-9])(?::([0-5][0-9]))?$", x))
  vapply(m, function(parts) {
    if (length(parts) == 0) return(NA_real_)
    secs <- if (nzchar(parts[4])) as.numeric(parts[4]) else 0
    as.numeric(parts[2]) * 3600 + as.numeric(parts[3]) * 60 + secs
  }, numeric(1))
}

# The window of time each camera contributes to the analysis, from its start
# date-time to its end date-time:
# - Start: Start Date at Start Time. If Start Time is missing or unreadable,
#   the window starts at midnight at the end of Start Date, so time before the
#   camera was set isn't counted.
# - End: End Date at End Time. If End Time is missing or unreadable, the
#   window ends at midnight at the start of End Date. For a malfunctioning
#   camera (Camera Functioning = No), the window ends at midnight at the start
#   of Camera Malfunction Date.
# - Trim: if max_days is set, the window is at most max_days x 24 hours.
# Date-times are built in UTC from the recorded clock values, the same way
# image timestamps are parsed, so the two are directly comparable.
analysis_windows <- function(deployments, max_days = NULL) {
  to_midnight <- function(d) as.POSIXct(format(d, "%Y-%m-%d"), tz = "UTC")

  start_date <- parse_mdY_strict(deployments$`Start Date`)
  end_date <- parse_mdY_strict(deployments$`End Date`)
  start_secs <- parse_time_of_day(deployments$`Start Time`)
  end_secs <- parse_time_of_day(deployments$`End Time`)
  malfunction_date <- parse_mdY_strict(deployments$`Camera Malfunction Date`)
  malfunction <- tolower(deployments$`Camera Functioning`) %in% "no" &
    !is.na(malfunction_date)

  window_start <- to_midnight(start_date) +
    ifelse(is.na(start_secs), 24 * 3600, start_secs)
  window_end <- to_midnight(end_date) + ifelse(is.na(end_secs), 0, end_secs)
  window_end[malfunction] <- to_midnight(malfunction_date[malfunction])

  trimmed <- rep(FALSE, length(window_start))
  if (!is.null(max_days) && is.finite(max_days) && max_days > 0) {
    trim_end <- window_start + as.integer(max_days[[1]]) * 24 * 3600
    trimmed <- !is.na(window_end) & window_end > trim_end
    window_end[trimmed] <- trim_end[trimmed]
  }

  data.frame(
    `Site Name` = deployments$`Site Name`,
    window_start = window_start,
    window_end = window_end,
    days = as.numeric(difftime(window_end, window_start, units = "days")),
    start_time_missing = is.na(start_secs),
    end_time_missing = is.na(end_secs) & !malfunction,
    malfunction = malfunction,
    trimmed = trimmed,
    check.names = FALSE
  )
}

# One row per detection event (Site Name + Cluster ID + Species), the unit the
# models count. Grouping by Site Name keeps a Cluster ID reused at two sites
# from merging. An event is in the camera's analysis window if its first image
# is; `count` is the largest Sighting Count among its images, and
# `days_in_window` is the part of its duration (first to last image) inside
# the window, which the models remove from the camera's effort. Rows without a
# Species aren't events.
detection_events <- function(images, deployments, max_days = NULL) {
  windows <- analysis_windows(deployments, max_days = max_days) |>
    select(`Site Name`, window_start, window_end)

  images |>
    filter(!is.na(Species), Species != "") |>
    group_by(`Site Name`, `Cluster ID`, Species) |>
    summarize(
      n_images = n(),
      count    = max(`Sighting Count`),
      start    = min(Timestamp),
      end      = max(Timestamp),
      .groups  = "drop"
    ) |>
    inner_join(windows, by = "Site Name") |>
    mutate(
      in_window = start >= window_start & start <= window_end,
      days_in_window = if_else(
        in_window,
        as.numeric(difftime(pmin(end, window_end), start, units = "days")),
        0
      )
    )
}

# -------------------------------------------------------------------
# Summary helpers (Data summary tab)
# -------------------------------------------------------------------

# Camera-days from each camera's analysis window (see analysis_windows()), so
# they match the windows the models use, including the optional trim. Time
# spent within detection events, which the models also remove, isn't
# subtracted here.
summarize_deployments <- function(deployments, max_days = NULL) {
  analysis_windows(deployments, max_days = max_days) |>
    summarize(
      n_cameras         = n(),
      mean_days         = mean(days),
      total_camera_days = sum(days)
    )
}

# Per species: images in the file, plus detection events, animals counted
# (sum of event counts), and cameras with events inside the analysis windows,
# all as the models count them.
summarize_species <- function(images, events) {
  image_counts <- images |>
    filter(!is.na(Species), Species != "") |>
    count(Species, name = "n_images")

  event_counts <- events |>
    filter(in_window) |>
    group_by(Species) |>
    summarize(
      n_events = n(),
      animals  = sum(count),
      cameras  = n_distinct(`Site Name`),
      .groups  = "drop"
    )

  image_counts |>
    left_join(event_counts, by = "Species") |>
    mutate(across(c(n_events, animals, cameras), \(x) coalesce(x, 0))) |>
    arrange(desc(animals), Species)
}

# Per camera, for one species: images in the file, plus detection events and
# animals counted inside the camera's analysis window. Every camera is listed,
# with zeros where the species wasn't detected.
species_by_camera <- function(images, events, deployments, species_name) {
  image_counts <- images |>
    filter(Species == species_name) |>
    count(`Site Name`, name = "n_images")

  event_counts <- events |>
    filter(Species == species_name, in_window) |>
    group_by(`Site Name`) |>
    summarize(n_events = n(), animals = sum(count), .groups = "drop")

  deployments |>
    transmute(
      `Site Name`,
      Latitude  = as.numeric(Latitude),
      Longitude = as.numeric(Longitude)
    ) |>
    left_join(image_counts, by = "Site Name") |>
    left_join(event_counts, by = "Site Name") |>
    mutate(across(c(n_images, n_events, animals), \(x) coalesce(x, 0)))
}

# Animals counted per camera per day, for one species, by the date of each
# event's first image (events inside the analysis windows only).
species_by_day <- function(events, species_name) {
  events |>
    filter(Species == species_name, in_window) |>
    mutate(Date = as.Date(format(start, "%Y-%m-%d"))) |>
    group_by(`Site Name`, Date) |>
    summarize(animals = sum(count), .groups = "drop")
}
