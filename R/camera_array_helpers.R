camera_array_supported_extensions <- function() {
  c("zip", "shp", "dbf", "shx", "prj", "cpg", "kml", "kmz", "gpkg")
}

stage_spatial_uploads <- function(upload_df) {
  if (is.null(upload_df) || !nrow(upload_df)) {
    return(NULL)
  }

  upload_dir <- file.path(
    tempdir(),
    paste0("camera_array_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sample.int(1e6, 1))
  )
  dir.create(upload_dir, recursive = TRUE, showWarnings = FALSE)

  staged_files <- character(nrow(upload_df))
  for (i in seq_len(nrow(upload_df))) {
    out_path <- file.path(upload_dir, upload_df$name[i])
    ok <- file.copy(upload_df$datapath[i], out_path, overwrite = TRUE)
    if (!isTRUE(ok)) {
      stop("Failed to stage uploaded file: ", upload_df$name[i], call. = FALSE)
    }
    staged_files[i] <- out_path
  }

  list(
    dir = upload_dir,
    files = staged_files
  )
}

read_spatial_entry <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "kmz") {
    kmz_dir <- file.path(tempdir(), paste0("kmz_", sample.int(1e6, 1)))
    dir.create(kmz_dir, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(path, exdir = kmz_dir)
    kml_files <- list.files(kmz_dir, pattern = "\\.kml$", recursive = TRUE, full.names = TRUE)
    if (!length(kml_files)) {
      stop("KMZ file did not contain a readable KML layer.", call. = FALSE)
    }
    return(read_spatial_entry(kml_files[1]))
  }

  if (ext == "kml") {
    layer_names <- tryCatch(sf::st_layers(path)$name, error = function(e) character())
    if (!length(layer_names)) {
      return(sf::st_read(path, quiet = TRUE))
    }
    return(
      dplyr::bind_rows(
        lapply(layer_names, function(layer_name) {
          sf::st_read(path, layer = layer_name, quiet = TRUE)
        })
      )
    )
  }

  sf::st_read(path, quiet = TRUE)
}

read_uploaded_spatial_group <- function(upload_df, label = "spatial layer", required = FALSE) {
  if (is.null(upload_df) || !nrow(upload_df)) {
    if (required) {
      stop("Please upload at least one ", label, " file.", call. = FALSE)
    }
    return(NULL)
  }

  stage <- stage_spatial_uploads(upload_df)
  staged_files <- stage$files
  staged_ext <- tolower(tools::file_ext(staged_files))

  valid_ext <- camera_array_supported_extensions()
  bad_ext <- setdiff(unique(staged_ext), valid_ext)
  if (length(bad_ext)) {
    stop(
      "Unsupported ", label, " file type(s): ",
      paste(bad_ext, collapse = ", "),
      ". Supported: ", paste(valid_ext, collapse = ", "),
      call. = FALSE
    )
  }

  zip_files <- staged_files[staged_ext == "zip"]
  direct_entries <- staged_files[staged_ext %in% c("shp", "kml", "kmz", "gpkg")]

  if (length(zip_files)) {
    for (zip_path in zip_files) {
      unzip_dir <- file.path(stage$dir, paste0(tools::file_path_sans_ext(basename(zip_path)), "_unzipped"))
      dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
      utils::unzip(zip_path, exdir = unzip_dir)
      unzip_entries <- list.files(
        unzip_dir,
        pattern = "\\.(shp|kml|kmz|gpkg)$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
      direct_entries <- c(direct_entries, unzip_entries)
    }
  }

  if (!length(direct_entries)) {
    if (required) {
      stop("No readable ", label, " files were found after upload.", call. = FALSE)
    }
    return(NULL)
  }

  layer_list <- lapply(direct_entries, function(one_path) {
    tryCatch(
      read_spatial_entry(one_path),
      error = function(e) {
        stop(
          "Could not read ", label, " file '", basename(one_path), "': ", e$message,
          call. = FALSE
        )
      }
    )
  })

  layer_list <- Filter(Negate(is.null), layer_list)
  if (!length(layer_list)) {
    if (required) {
      stop("No valid ", label, " layers could be read.", call. = FALSE)
    }
    return(NULL)
  }

  out <- dplyr::bind_rows(layer_list)
  out <- sf::st_make_valid(out)

  if (is.na(sf::st_crs(out))) {
    stop("The uploaded ", label, " must have a defined coordinate reference system (CRS).", call. = FALSE)
  }

  out
}

validate_camera_boundary <- function(boundary_sf) {
  if (!inherits(boundary_sf, "sf")) {
    stop("Boundary input is not a valid sf object.", call. = FALSE)
  }
  if (is.null(sf::st_geometry(boundary_sf))) {
    stop("Boundary input does not have a geometry column.", call. = FALSE)
  }

  geom_type <- unique(as.character(sf::st_geometry_type(boundary_sf)))
  if (!all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("Boundary files must contain polygon or multipolygon geometry.", call. = FALSE)
  }

  if (is.na(sf::st_crs(boundary_sf))) {
    stop("Boundary files must have a defined CRS.", call. = FALSE)
  }

  boundary_sf
}

determine_utm_crs <- function(boundary_sf) {
  boundary_ll <- sf::st_transform(boundary_sf, 4326)
  centroid <- suppressWarnings(sf::st_centroid(sf::st_union(sf::st_geometry(boundary_ll))))
  coords <- sf::st_coordinates(centroid)
  lon <- coords[1, 1]
  lat <- coords[1, 2]
  utm_zone <- floor((lon + 180) / 6) + 1

  if (lat >= 0) {
    paste0("EPSG:", 269, sprintf("%02d", utm_zone))
  } else {
    paste0("EPSG:", 327, sprintf("%02d", utm_zone))
  }
}

transform_to_boundary_crs <- function(layer_sf, boundary_crs) {
  if (is.null(layer_sf)) return(NULL)
  if (is.na(sf::st_crs(layer_sf))) {
    stop("Uploaded spatial layers must have a defined CRS.", call. = FALSE)
  }
  if (sf::st_crs(layer_sf) != sf::st_crs(boundary_crs)) {
    return(sf::st_transform(layer_sf, sf::st_crs(boundary_crs)))
  }
  layer_sf
}

build_camera_design_area <- function(boundary_sf,
                                     roads = NULL,
                                     trails = NULL,
                                     buildings = NULL,
                                     parking_lots = NULL,
                                     exclusion_areas = NULL,
                                     buffer_dist_m = 10) {
  boundary_sf <- validate_camera_boundary(boundary_sf)
  boundary_utm <- sf::st_transform(boundary_sf, determine_utm_crs(boundary_sf))
  boundary_utm <- sf::st_make_valid(boundary_utm)
  boundary_utm <- sf::st_sf(geometry = sf::st_union(boundary_utm))

  optional_layers <- list(
    roads = roads,
    trails = trails,
    buildings = buildings,
    parking_lots = parking_lots
  )
  optional_layers <- lapply(optional_layers, transform_to_boundary_crs, boundary_crs = boundary_utm)
  exclusion_areas <- transform_to_boundary_crs(exclusion_areas, boundary_crs = boundary_utm)

  exclusion_parts <- list()

  for (layer_name in names(optional_layers)) {
    layer_sf <- optional_layers[[layer_name]]
    if (!is.null(layer_sf) && nrow(layer_sf) > 0) {
      buffered <- sf::st_buffer(layer_sf, dist = buffer_dist_m)
      buffered <- suppressWarnings(sf::st_intersection(buffered, boundary_utm))
      if (nrow(buffered) > 0) {
        exclusion_parts[[layer_name]] <- sf::st_make_valid(buffered)
      }
    }
  }

  if (!is.null(exclusion_areas) && nrow(exclusion_areas) > 0) {
    exclusion_parts[["uploaded_exclusions"]] <- sf::st_make_valid(exclusion_areas)
  }

  exclusion_union <- NULL
  if (length(exclusion_parts)) {
    exclusion_union <- sf::st_sf(
      geometry = sf::st_union(do.call(c, lapply(exclusion_parts, sf::st_geometry)))
    )
    exclusion_union <- sf::st_make_valid(exclusion_union)
  }

  allowed_area <- boundary_utm
  if (!is.null(exclusion_union) && nrow(exclusion_union) > 0) {
    allowed_area <- suppressWarnings(sf::st_difference(boundary_utm, exclusion_union))
    allowed_area <- sf::st_make_valid(allowed_area)
  }

  allowed_area <- allowed_area[!sf::st_is_empty(allowed_area), , drop = FALSE]
  if (!nrow(allowed_area)) {
    stop("No allowable camera area remained after subtracting exclusions and buffers.", call. = FALSE)
  }

  list(
    boundary = boundary_utm,
    allowed_area = allowed_area,
    excluded_area = exclusion_union
  )
}

generate_camera_candidates <- function(allowed_area, spacing_m, camera_prefix = "CAM") {
  if (!is.numeric(spacing_m) || length(spacing_m) != 1 || is.na(spacing_m) || spacing_m <= 0) {
    stop("Camera spacing must be a positive number of meters.", call. = FALSE)
  }

  union_geom <- sf::st_union(allowed_area)
  grid_pts <- sf::st_make_grid(union_geom, cellsize = spacing_m, what = "centers", square = TRUE)
  inside <- lengths(sf::st_intersects(grid_pts, union_geom)) > 0
  grid_pts <- grid_pts[inside]

  if (!length(grid_pts)) {
    stop("No candidate camera locations were generated. Try a smaller spacing value.", call. = FALSE)
  }

  pts_sf <- sf::st_sf(geometry = grid_pts)
  coords <- sf::st_coordinates(pts_sf)
  pts_sf$x_utm <- coords[, 1]
  pts_sf$y_utm <- coords[, 2]
  pts_sf <- pts_sf[order(-pts_sf$y_utm, pts_sf$x_utm), , drop = FALSE]
  pts_sf$camera_id <- sprintf("%s_%03d", camera_prefix, seq_len(nrow(pts_sf)))
  pts_sf
}

nearest_neighbor_distances_m <- function(camera_sf) {
  if (is.null(camera_sf) || nrow(camera_sf) < 2) {
    return(rep(NA_real_, nrow(camera_sf)))
  }

  coords <- sf::st_coordinates(camera_sf)
  out <- rep(NA_real_, nrow(coords))
  for (i in seq_len(nrow(coords))) {
    deltas <- coords[-i, , drop = FALSE] - matrix(coords[i, ], nrow = nrow(coords) - 1L, ncol = 2L, byrow = TRUE)
    out[i] <- min(sqrt(rowSums(deltas^2)))
  }
  out
}

camera_grid_indices <- function(camera_sf) {
  if (is.null(camera_sf) || !nrow(camera_sf)) {
    return(data.frame(row_idx = integer(), col_idx = integer()))
  }

  coords <- if (all(c("x_utm", "y_utm") %in% names(camera_sf))) {
    cbind(camera_sf$x_utm, camera_sf$y_utm)
  } else {
    sf::st_coordinates(camera_sf)
  }

  x_vals <- sort(unique(round(coords[, 1], 6)))
  y_vals <- sort(unique(round(coords[, 2], 6)), decreasing = TRUE)

  data.frame(
    row_idx = match(round(coords[, 2], 6), y_vals),
    col_idx = match(round(coords[, 1], 6), x_vals)
  )
}

camera_edge_priority <- function(camera_sf) {
  idx <- camera_grid_indices(camera_sf)
  if (!nrow(idx)) {
    return(data.frame(edge_distance = numeric(), corner_distance = numeric()))
  }

  row_min <- min(idx$row_idx)
  row_max <- max(idx$row_idx)
  col_min <- min(idx$col_idx)
  col_max <- max(idx$col_idx)

  edge_distance <- pmin(
    idx$row_idx - row_min,
    row_max - idx$row_idx,
    idx$col_idx - col_min,
    col_max - idx$col_idx
  )

  corner_distance <- pmin(
    abs(idx$row_idx - row_min) + abs(idx$col_idx - col_min),
    abs(idx$row_idx - row_min) + abs(idx$col_idx - col_max),
    abs(idx$row_idx - row_max) + abs(idx$col_idx - col_min),
    abs(idx$row_idx - row_max) + abs(idx$col_idx - col_max)
  )

  data.frame(
    edge_distance = edge_distance,
    corner_distance = corner_distance
  )
}

select_camera_subset <- function(candidates_sf, budget = 0, n_alternates = 2, spacing_m = NA_real_) {
  total_n <- nrow(candidates_sf)
  if (!total_n) {
    stop("No candidate cameras were available for selection.", call. = FALSE)
  }

  budget <- as.integer(budget)
  n_alternates <- as.integer(max(0, n_alternates))
  prefix <- if ("camera_id" %in% names(candidates_sf)) {
    sub("_[0-9]+$", "", candidates_sf$camera_id[1])
  } else {
    "CAM"
  }
  candidates_sf$candidate_camera_id <- candidates_sf$camera_id

  if (is.na(budget) || budget <= 0 || budget >= total_n) {
    selected <- candidates_sf
    selected$status <- "final"
    selected$camera_id <- sprintf("%s_%03d", prefix, seq_len(nrow(selected)))
    nn_selected <- nearest_neighbor_distances_m(selected)
    mean_nn <- if (length(nn_selected) && any(is.finite(nn_selected))) mean(nn_selected, na.rm = TRUE) else NA_real_
    max_nn <- if (length(nn_selected) && any(is.finite(nn_selected))) max(nn_selected, na.rm = TRUE) else NA_real_
    attr(selected, "spacing_note") <- ""
    attr(selected, "mean_nearest_neighbor_m") <- mean_nn
    attr(selected, "max_nearest_neighbor_m") <- max_nn
    return(selected)
  }

  n_target <- min(total_n, budget + n_alternates)
  coords <- sf::st_coordinates(candidates_sf)
  center <- colMeans(coords)
  first_idx <- which.min(rowSums((coords - matrix(center, nrow(coords), 2, byrow = TRUE))^2))

  chosen <- first_idx
  while (length(chosen) < n_target) {
    remaining <- setdiff(seq_len(total_n), chosen)
    dist_to_selected <- sapply(remaining, function(idx) {
      min(sqrt(rowSums((coords[chosen, , drop = FALSE] - matrix(coords[idx, ], nrow = length(chosen), ncol = 2, byrow = TRUE))^2)))
    })
    chosen <- c(chosen, remaining[which.max(dist_to_selected)])
  }

  out <- candidates_sf[chosen, , drop = FALSE]
  helper_idx <- camera_grid_indices(out)
  edge_priority <- camera_edge_priority(out)
  out$row_idx_helper <- helper_idx$row_idx
  out$col_idx_helper <- helper_idx$col_idx
  out$edge_distance_helper <- edge_priority$edge_distance
  out$corner_distance_helper <- edge_priority$corner_distance
  out$status <- "final"

  n_final <- min(budget, nrow(out))
  n_alt <- max(0, nrow(out) - n_final)
  if (n_alt > 0) {
    alt_order <- order(
      out$corner_distance_helper,
      out$edge_distance_helper,
      out$row_idx_helper,
      out$col_idx_helper
    )
    out$status[alt_order[seq_len(n_alt)]] <- "alternate"
  }

  final_pts <- out[out$status == "final", , drop = FALSE]
  alt_pts <- out[out$status == "alternate", , drop = FALSE]

  if (nrow(final_pts)) {
    final_pts <- final_pts[order(final_pts$row_idx_helper, final_pts$col_idx_helper), , drop = FALSE]
  }
  if (nrow(alt_pts)) {
    alt_pts <- alt_pts[order(alt_pts$corner_distance_helper, alt_pts$edge_distance_helper,
                             alt_pts$row_idx_helper, alt_pts$col_idx_helper), , drop = FALSE]
  }

  if (nrow(final_pts)) {
    final_pts$camera_id <- sprintf("%s_%03d", prefix, seq_len(nrow(final_pts)))
  }
  if (nrow(alt_pts)) {
    alt_pts$camera_id <- sprintf("%s_%03d", prefix, 100L + seq_len(nrow(alt_pts)))
  }

  out <- dplyr::bind_rows(final_pts, alt_pts)
  out$row_idx_helper <- NULL
  out$col_idx_helper <- NULL
  out$edge_distance_helper <- NULL
  out$corner_distance_helper <- NULL

  nn_final <- nearest_neighbor_distances_m(final_pts)
  mean_nn <- if (length(nn_final) && any(is.finite(nn_final))) mean(nn_final, na.rm = TRUE) else NA_real_
  max_nn <- if (length(nn_final) && any(is.finite(nn_final))) max(nn_final, na.rm = TRUE) else NA_real_
  spacing_note <- ""
  if (is.finite(spacing_m) && is.finite(max_nn) && max_nn > spacing_m * 1.5) {
    spacing_note <- paste0(
      "The reduced camera count keeps a well-spread subset, but it does not preserve the original ",
      round(spacing_m),
      " m spacing target everywhere. The largest nearest-neighbor spacing among final cameras is about ",
      round(max_nn),
      " m."
    )
  }

  attr(out, "spacing_note") <- spacing_note
  attr(out, "mean_nearest_neighbor_m") <- mean_nn
  attr(out, "max_nearest_neighbor_m") <- max_nn
  out
}

camera_points_for_export <- function(camera_sf) {
  if (is.null(camera_sf) || !nrow(camera_sf)) return(NULL)

  camera_ll <- sf::st_transform(camera_sf, 4326)
  ll_coords <- sf::st_coordinates(camera_ll)
  utm_coords <- sf::st_coordinates(camera_sf)

  out <- camera_sf
  out$longitude <- ll_coords[, 1]
  out$latitude <- ll_coords[, 2]
  out$utm_easting_m <- utm_coords[, 1]
  out$utm_northing_m <- utm_coords[, 2]
  out
}
