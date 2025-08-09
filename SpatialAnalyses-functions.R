

#corrected rgbif functions (fixing bugs after updating some spatial dependencies)
{
  
  ## --------------------------------------------------------------
  ## 0. Load rgbif – we rely on its internal helpers
  ## --------------------------------------------------------------
  library(rgbif)
  
  ## --------------------------------------------------------------
  ## 1. Patched async helper (forces HTTP/1.1)
  ## --------------------------------------------------------------
  gbif_async_get <- function(urls, parse = FALSE, curlopts = list()) {
    message("✅ using patched gbif_async_get (HTTP/1.1)")
    cc <- crul::Async$new(
      urls    = urls,
      headers = rgbif:::rgbif_ual,
      opts    = c(list(http_version = 1L), curlopts)
    )
    rgbif:::process_async_get(cc$get(), parse = parse)
  }
  
  ## --------------------------------------------------------------
  ## 2. check_name_data – preserve raw, add cleaned copy for GBIF
  ## --------------------------------------------------------------
  check_name_data <- function(name_data) {
    if (is.null(name_data)) stop("name_data is missing")
    
    to_df <- function(x) {
      data.frame(
        name_raw = x,
        name     = gsub("_", " ", x),   # <-- for matching
        stringsAsFactors = FALSE
      )
    }
    
    ## a) vector input ------------------------------------------------
    if (is.vector(name_data)) {
      stopifnot(is.character(name_data))
      out <- to_df(name_data)
      out$index <- seq_len(nrow(out))
      return(out)
    }
    
    ## b) one-col data.frame -----------------------------------------
    if (ncol(name_data) == 1) {
      colnames(name_data) <- "name_raw"
      stopifnot(is.character(name_data$name_raw))
      name_data$name  <- gsub("_", " ", name_data$name_raw)
      name_data$index <- seq_len(nrow(name_data))
      return(name_data)
    }
    
    ## c) multi-col data.frame ---------------------------------------
    orig <- colnames(name_data)
    colnames(name_data) <- tolower(gsub("_", "", orig))
    
    alias <- c("scientificname","sciname","names","species","speciesname",
               "spname","taxonname")
    if (!"name_raw" %in% colnames(name_data) && any(alias %in% colnames(name_data))) {
      colnames(name_data)[which(colnames(name_data) %in% alias)[1]] <- "name_raw"
    }
    if (!"name_raw" %in% colnames(name_data))
      stop("No column recognised as species name")
    
    name_data$name  <- gsub("_", " ", name_data$name_raw)
    name_data$index <- seq_len(nrow(name_data))
    name_data
  }
  
  ## --------------------------------------------------------------
  ## 3. Bring other internal helpers into the workspace
  ## --------------------------------------------------------------
  bind_rows            <- rgbif:::bind_rows
  process_alternatives <- rgbif:::process_alternatives
  
  ## --------------------------------------------------------------
  ## 4. Original name_backbone_checklist (only namespaced tweaks)
  ## --------------------------------------------------------------
  name_backbone_checklist <- function(
    name_data = NULL,
    rank = NULL, kingdom = NULL, phylum = NULL, class = NULL,
    order = NULL, family = NULL, genus = NULL,
    strict = FALSE, verbose = FALSE, curlopts = list()
  ) {
    name_data <- check_name_data(name_data)
    if (!is.null(c(rank, kingdom, phylum, class, order, family, genus)))
      name_data <- rgbif:::default_value_handler(
        name_data, rank, kingdom, phylum, class, order, family, genus)
    
    ## build list for GBIF
    data_list <- lapply(
      data.table::transpose(name_data[ , !names(name_data) %in% "name_raw"]),
      function(x) stats::setNames(as.list(x),
                                  colnames(name_data)[!colnames(name_data) %in% "name_raw"])
    )
    urls         <- rgbif:::make_async_urls(data_list, verbose, strict)
    matched_list <- gbif_async_get(urls, curlopts = curlopts)
    
    ## merge with verbatim
    verbatim_list <- lapply(
      data.table::transpose(name_data[ , c("name_raw","index") ]),
      function(x) stats::setNames(list(x[[1]], x[[2]]),
                                  c("verbatim_name", "verbatim_index"))
    )
    mvl           <- mapply(function(x,y) c(x,y), verbatim_list, matched_list,
                            SIMPLIFY = FALSE)
    matched_names <- bind_rows(mvl)
    
    if (verbose) {
      alts <- lapply(mvl, function(x) x[["alternatives"]])
      matched_names <- bind_rows(list(
        matched_names,
        process_alternatives(alts, verbatim_list)))
    }
    
    matched_names$verbatim_index <- as.numeric(matched_names$verbatim_index)
    matched_names <- matched_names[order(matched_names$verbatim_index), ]
    matched_names <- matched_names[!names(matched_names) %in% c("alternatives","note")]
    
    col_idx <- grep("verbatim_", names(matched_names))
    matched_names <- unique(
      matched_names[, c((1:ncol(matched_names))[-col_idx], col_idx)]
    )
    ## -- no second re-sort here, keeping historical row order
    
    if (verbose && "is_alternative" %in% names(matched_names))
      matched_names$is_alternative[is.na(matched_names$is_alternative)] <- FALSE
    
    matched_names
  }
  
  
  
}

#function for merging polygons when including multiple ranges per speciesKey

#' Merge Overlapping Species Polygons (Spherical)
#'
#' Merges duplicate polygon geometries per species (identified by a key column)
#' using `sf` geometry operations under spherical (s2) mode. Each duplicate
#' species gets a single merged geometry written back to all of its rows.
merge_species_shapes <- function(data, speciesKey_col = "speciesKey", shape_col = "Shape") {
  library(sf)
  
  # Ensure all geometries in the dataset are valid
  data[[shape_col]] <- st_make_valid(data[[shape_col]])
  
  # Identify species with duplicate entries
  species_with_duplicates <- unique(data[[speciesKey_col]][duplicated(data[[speciesKey_col]])])
  message("Found ", length(species_with_duplicates), " species with duplicates to process...")
  
  # Iterate over each duplicated species key
  for (sp_key in species_with_duplicates) {
    message("--------------------------------------------------")
    message("Processing speciesKey: ", sp_key)
    
    # Extract row indices for the current speciesKey
    indices <- which(data[[speciesKey_col]] == sp_key)
    subset_shapes <- data[[shape_col]][indices]
    
    # Ensure all geometries in the subset are valid
    valid_geometries <- lapply(subset_shapes, st_make_valid)
    valid_collection <- st_sfc(valid_geometries, crs = st_crs(data))
    
    # Merge geometries using st_union with buffered tolerance
    merged_geom <- tryCatch({
      st_make_valid(st_union(st_buffer(valid_collection, dist = 0.0000001)))
    }, error = function(e) {
      message(" Error merging speciesKey=", sp_key, ": ", e$message)
      return(NULL)
    })
    
    # Validate the merged geometry
    if (is.null(merged_geom) || st_is_empty(merged_geom) || !st_is_valid(merged_geom)) {
      message(" Merged geometry for speciesKey=", sp_key, " is invalid or empty. Skipping.")
      next
    }
    
    # Replace each row's geometry for the speciesKey
    for (idx in indices) {
      data[[shape_col]][idx] <- merged_geom[[1]]
      message(" Replaced geometry for row index: ", idx)
    }
    
    message(" Successfully merged geometries for speciesKey: ", sp_key)
  }
  
  #fixing CRS
  data$Shape <- st_sfc(data$Shape, crs='WGS84')
  
  # Final validity check
  final_validities <- sapply(data[[shape_col]], st_is_valid)
  if (!all(final_validities)) {
    invalid_idx <- which(!final_validities)
    message("Some geometries remain invalid after merging:")
    print(invalid_idx)
    stop("Some geometries remain invalid. Please investigate.")
  }
  
  message("Merging complete. Returning updated dataset.")
  return(data)
}

#' Merge Overlapping Species Polygons (Planar / s2 Disabled)
#'
#' Same as [merge_species_shapes()], but temporarily disables s2 spherical geometry
#' (`sf_use_s2(FALSE)`) to run planar unions (useful when s2 intersects/union fails).
merge_species_shapes_planar <- function(data, speciesKey_col = "speciesKey", shape_col = "Shape") {
  library(sf)
  
  # Ensure all geometries in the dataset are valid
  data[[shape_col]] <- st_make_valid(data[[shape_col]])
  
  # Identify species with duplicate entries
  species_with_duplicates <- unique(data[[speciesKey_col]][duplicated(data[[speciesKey_col]])])
  message("Found ", length(species_with_duplicates), " species with duplicates to process...")
  
  # Temporarily switch off s2 spherical geometry
  sf_use_s2(FALSE)
  
  # Iterate over each duplicated species key
  for (sp_key in species_with_duplicates) {
    message("--------------------------------------------------")
    message("Processing speciesKey: ", sp_key)
    
    # Extract row indices for the current speciesKey
    indices <- which(data[[speciesKey_col]] == sp_key)
    subset_shapes <- data[[shape_col]][indices]
    
    # Ensure all geometries in the subset are valid
    valid_geometries <- lapply(subset_shapes, st_make_valid)
    valid_collection <- st_sfc(valid_geometries, crs = st_crs(data))
    
    # Merge geometries using st_union without spherical geometry
    merged_geom <- tryCatch({
      st_make_valid(st_union(valid_collection))
    }, error = function(e) {
      message(" Error merging speciesKey=", sp_key, ": ", e$message)
      return(NULL)
    })
    
    # Validate the merged geometry
    if (is.null(merged_geom) || st_is_empty(merged_geom) || !st_is_valid(merged_geom)) {
      message(" Merged geometry for speciesKey=", sp_key, " is invalid or empty. Skipping.")
      next
    }
    
    # Replace each row's geometry for the speciesKey
    for (idx in indices) {
      data[[shape_col]][idx] <- merged_geom[[1]]
      message(" Replaced geometry for row index: ", idx)
    }
    
    message(" Successfully merged geometries for speciesKey: ", sp_key)
  }
  
  # Restore default s2 spherical geometry behavior
  sf_use_s2(TRUE)
  
  #fixing CRS
  data$Shape <- st_sfc(data$Shape, crs='WGS84')
  
  # Final validity check
  final_validities <- sapply(data[[shape_col]], st_is_valid)
  if (!all(final_validities)) {
    invalid_idx <- which(!final_validities)
    message("Some geometries remain invalid after merging:")
    print(invalid_idx)
    stop("Some geometries remain invalid. Please investigate.")
  }
  
  message("Merging complete. Returning updated dataset.")
  return(data)
}

#' Parallel Merge of Species Polygons (Planar)
#'
#' Validates all geometries, merges duplicates per species key in parallel with
#' `pbmcapply::pbmclapply()` while s2 is disabled, then writes merged geometries
#' back to the dataset with a progress bar. Designed for large datasets.
merge_species_shapes_planar_parallel <- function(data, speciesKey_col = "speciesKey", shape_col = "Shape", workers = 1) {
  library(sf)
  library(pbmcapply)  # Use pbmclapply for parallel processing with a progress bar
  library(pbapply)    # Use pblapply for progress bar in Step 4
  
  #options(digits = 22)
  
  # Step 1: Validate all geometries
  message("Step 1: Validating all geometries in the dataset...")
  data[[shape_col]] <- pbmclapply(data[[shape_col]], st_make_valid, mc.cores = workers)
  data[[shape_col]] <- do.call(st_sfc, data[[shape_col]])
  message("Step 1 complete: All geometries validated.")
  
  # Step 2: Identify species with duplicate entries
  species_with_duplicates <- unique(data[[speciesKey_col]][duplicated(data[[speciesKey_col]])])
  #species_with_duplicates <- species_with_duplicates[!is.na(species_with_duplicates)] #filter out any NAs
  message("Step 2: Found ", length(species_with_duplicates), " species with duplicates to process.")
  
  # Temporarily switch off s2 spherical geometry
  sf_use_s2(FALSE)
  
  # Step 3: Process each duplicated species in parallel
  message("Step 3: Merging geometries for duplicated species...")
  updated_geometries <- pbmclapply(species_with_duplicates, function(sp_key) {
    # Extract row indices for the current speciesKey
    indices <- which(data[[speciesKey_col]] == sp_key)
    subset_shapes <- data[[shape_col]][indices]
    
    # Ensure all geometries in the subset are valid
    valid_geometries <- lapply(subset_shapes, st_make_valid)
    valid_collection <- st_sfc(valid_geometries, crs = st_crs(data))
    
    # Merge geometries using st_union without spherical geometry
    merged_geom <- tryCatch({
      st_make_valid(st_union(valid_collection))
    }, error = function(e) {
      return(NULL)
    })
    
    # Validate the merged geometry
    if (is.null(merged_geom) || st_is_empty(merged_geom) || !st_is_valid(merged_geom)) {
      return(NULL)
    }
    
    list(indices = indices, merged_geom = merged_geom)
  }, mc.cores = workers)
  
  # updated_geometries <- pbmclapply(species_with_duplicates, function(sp_key) {
  #   tryCatch({
  #     # Extract row indices for the current speciesKey
  #     indices <- which(data[[speciesKey_col]] == sp_key)
  #     subset_shapes <- data[[shape_col]][indices]
  # 
  #     # Ensure all geometries in the subset are valid
  #     valid_geometries <- lapply(subset_shapes, st_make_valid)
  #     valid_collection <- st_sfc(valid_geometries, crs = st_crs(data))
  # 
  #     # Merge geometries using st_union without spherical geometry
  #     merged_geom <- st_make_valid(st_union(valid_collection))
  # 
  #     # Validate the merged geometry
  #     if (is.null(merged_geom) || st_is_empty(merged_geom) || !st_is_valid(merged_geom)) {
  #       message("Merged geometry is invalid, empty, or contains NA for speciesKey=", sp_key)
  #       return(NULL)
  #     }
  # 
  #     list(indices = indices, merged_geom = merged_geom)
  #   }, error = function(e) {
  #     # Handle any error by returning NULL and logging the speciesKey
  #     message("Error encountered for speciesKey=", sp_key, ": ", e$message)
  #     return(NULL)
  #   })
  # }, mc.cores = workers)
  #return(updated_geometries)
  
  
  # Filter valid results and report invalid entries
  valid_results <- Filter(Negate(is.null), updated_geometries)
  
  # Identify and log invalid species keys
  invalid_species <- species_with_duplicates[sapply(updated_geometries, is.null)]
  
  # Report the number of valid and invalid entries
  message("Number of valid species processed: ", length(valid_results))
  message("Number of invalid species skipped: ", length(invalid_species))
  
  # Log the invalid species keys for debugging
  if (length(invalid_species) > 0) {
    message("Invalid species keys: ", paste(invalid_species, collapse = ", "))
  }
  
  message("Step 3 complete: Finished merging geometries for all duplicated species.")
  #return(updated_geometries)
  
  # # Step 4: Update the dataset sequentially
  # message("Step 4: Updating dataset with merged geometries...")
  # for (result in updated_geometries) {
  #   if (!is.null(result)) {
  #     for (idx in result$indices) {
  #       data[[shape_col]][idx] <- result$merged_geom[[1]]
  #     }
  #   }
  # }
  # message("Step 4 complete: Dataset updated with merged geometries.")
  
  # # Step 4: Update the dataset with progress bar
  # message("Step 4: Updating dataset with merged geometries...")
  # pblapply(updated_geometries, function(result) {
  #   if (!is.null(result)) {
  #     for (idx in result$indices) {
  #       data[[shape_col]][idx] <- result$merged_geom[[1]]
  #     }
  #   }
  # })
  # message("Step 4 complete: Dataset updated with merged geometries.")
  
  # # Step 4: Update the dataset in batches
  # message("Step 4: Updating dataset with merged geometries...")
  # 
  # # Create a copy of the Shape column for batch updates
  # updated_shape <- data[[shape_col]]
  # 
  # # Apply updates to the Shape column
  # for (result in updated_geometries) {
  #   if (!is.null(result)) {
  #     updated_shape[result$indices] <- list(result$merged_geom[[1]])
  #   }
  # }
  # 
  # # Assign the updated Shape column back to the dataset
  # data[[shape_col]] <- st_sfc(updated_shape, crs = st_crs(data))
  # 
  # message("Step 4 complete: Dataset updated with merged geometries.")
  
  
  # Step 4: Update the dataset with merged geometries using a dynamic progress bar
  {
    message("Step 4: Updating dataset with merged geometries...")
    
    # Create a copy of the Shape column for batch updates
    updated_shape <- data[[shape_col]]
    
    # Total number of updates
    total_updates <- length(updated_geometries)
    
    # Initialize progress bar
    cat("Progress: [", strrep(" ", 50), "] 0%", sep = "")
    cat("\rProgress: [", sep = "")  # Move to the beginning of the progress bar line
    
    # Apply updates to the Shape column
    for (i in seq_along(updated_geometries)) {
      result <- updated_geometries[[i]]
      
      if (!is.null(result)) {
        updated_shape[result$indices] <- list(result$merged_geom[[1]])
      }
      
      # Update progress bar
      percent_complete <- (i / total_updates) * 100
      bar_length <- round((i / total_updates) * 50)  # Progress bar length (out of 50)
      cat(strrep("=", bar_length), strrep(" ", 50 - bar_length), "] ", sprintf("%.0f", percent_complete), "%", sep = "")
      if (i != total_updates) {
        cat("\rProgress: [", sep = "")  # Move the cursor back to the beginning of the bar
      }
    }
    
    # Complete the progress bar
    cat("\n")  # Move to a new line after completion
    
    # Assign the updated Shape column back to the dataset
    data[[shape_col]] <- st_sfc(updated_shape, crs = st_crs(data))
    
    message("Step 4 complete: Dataset updated with merged geometries.")
  }
  # Restore default s2 spherical geometry behavior
  sf_use_s2(TRUE)
  
  #fixing CRS
  data$Shape <- st_sfc(data$Shape, crs='WGS84')
  
  # Step 5: Perform final validity check in parallel
  message("Step 5: Performing final validity check...")
  final_validities <- pbmclapply(data[[shape_col]], st_is_valid, mc.cores = workers)
  if (!all(unlist(final_validities))) {
    stop("Step 5 failed: Some geometries remain invalid. Please investigate.")
  }
  message("Step 5 complete: All geometries are valid.")
  
  # Final reporting
  message("Merging process complete. Returning updated dataset.")
  
  return(data)
}

#' Collapse Duplicate Species Rows While Preserving Attributes
#'
#' Reduces multiple rows per species key into a single row. The geometry is taken
#' from the first row for that key, and *all other columns* are collapsed into
#' list-columns containing the set of values across the grouped rows.
collapse_species_data <- function(data, speciesKey_col = "speciesKey", shape_col = "Shape", verbose = TRUE) {
  library(sf)
  
  # Ensure input is a data.frame or tibble
  stopifnot(is.data.frame(data))
  if (!all(c(speciesKey_col, shape_col) %in% names(data))) {
    stop("Data does not contain specified speciesKey or Shape columns.")
  }
  
  # Identify all columns to collapse (besides speciesKey & Shape)
  all_cols <- names(data)
  other_cols <- setdiff(all_cols, c(speciesKey_col, shape_col))
  
  # Get unique speciesKeys
  unique_keys <- unique(data[[speciesKey_col]])
  if (verbose) {
    message("Number of unique speciesKey values: ", length(unique_keys))
  }
  
  # Initialize list to store results and counters for validation
  collapsed_data <- list()
  total_invalid <- 0
  total_fixed <- 0
  
  # Process each speciesKey sequentially
  for (sp_key in unique_keys) {
    if (verbose) {
      message("Processing speciesKey: ", sp_key)
    }
    
    # Subset rows for the current speciesKey
    subset_rows <- data[data[[speciesKey_col]] == sp_key, ]
    if (verbose) {
      message("  Number of rows for this speciesKey: ", nrow(subset_rows))
    }
    
    # Take the first geometry from the Shape column
    shape_value <- subset_rows[[shape_col]][1]
    
    # # Validate and fix the geometry if needed
    # if (!st_is_valid(shape_value)) {
    #   total_invalid <- total_invalid + 1
    #   if (verbose) {
    #     message("  Geometry is invalid. Attempting to fix...")
    #   }
    #   shape_value <- st_make_valid(shape_value)
    #   if (!st_is_valid(shape_value)) {
    #     if (verbose) {
    #       message("  Failed to fix geometry for speciesKey: ", sp_key)
    #     }
    #     next  # Skip this speciesKey if the geometry remains invalid
    #   }
    #   total_fixed <- total_fixed + 1
    #   if (verbose) {
    #     message("  Geometry successfully fixed.")
    #   }
    # } else if (verbose) {
    #   message("  Geometry is valid.")
    # }
    
    # Construct a 1-row data frame for the speciesKey
    row_df <- data.frame(
      speciesKey = sp_key,
      stringsAsFactors = FALSE
    )
    
    # Add other columns as lists of all values from the grouped rows
    for (colname in other_cols) {
      column_values <- subset_rows[[colname]]
      if (is.factor(column_values)) {
        column_values <- as.character(column_values)  # Convert factors to characters
      }
      row_df[[colname]] <- list(column_values)  # Store all values in a list
    }
    
    # Add the geometry directly into the data frame
    row_df[[shape_col]] <- shape_value
    
    # Append the row data frame to the results list
    collapsed_data[[as.character(sp_key)]] <- row_df
  }
  
  # Combine the list of one-row data frames into a single data frame
  if (verbose) {
    message("Combining collapsed data into a single data frame...")
  }
  result <- rbindlist(collapsed_data, ignore.attr = TRUE)
  
  # Final summary of geometry validation
  if (verbose) {
    message("Summary of geometry validation:")
    message("  Total geometries processed: ", length(unique_keys))
    #message("  Total invalid geometries detected: ", total_invalid)
    #message("  Total geometries successfully fixed: ", total_fixed)
    #if (total_invalid > total_fixed) {
    #  message("  Note: Some geometries could not be fixed and were excluded from the output.")
    #}
  }
  
  # Final message
  if (verbose) {
    message("Collapsing process complete. Returning consolidated data frame.")
  }
  
  return(result)
}

#' Parallel Collapse of Duplicate Species Rows
#'
#' Parallel version of [collapse_species_data()] using `pbmcapply::pbmclapply()`.
#' Each species key is processed in a separate worker and results are combined.
collapse_species_data_parallel <- function(data, speciesKey_col = "speciesKey", shape_col = "Shape", workers = 1, verbose = TRUE) {
  library(sf)
  library(pbmcapply)  # Parallel processing with progress bar
  
  # Ensure input is a data.frame or tibble
  stopifnot(is.data.frame(data))
  if (!all(c(speciesKey_col, shape_col) %in% names(data))) {
    stop("Data does not contain specified speciesKey or Shape columns.")
  }
  
  # Identify all columns to collapse (besides speciesKey & Shape)
  all_cols <- names(data)
  other_cols <- setdiff(all_cols, c(speciesKey_col, shape_col))
  
  # Get unique speciesKeys
  unique_keys <- unique(data[[speciesKey_col]])
  if (verbose) {
    message("Number of unique speciesKey values: ", length(unique_keys))
  }
  
  # Function to process each speciesKey
  process_speciesKey <- function(sp_key) {
    if (verbose) {
      message("Processing speciesKey: ", sp_key)
    }
    
    # Subset rows for the current speciesKey
    subset_rows <- data[data[[speciesKey_col]] == sp_key, ]
    if (verbose) {
      message("  Number of rows for this speciesKey: ", nrow(subset_rows))
    }
    
    # Take the first geometry from the Shape column
    shape_value <- subset_rows[[shape_col]][1]
    
    # # Validate and fix the geometry if needed
    # if (!st_is_valid(shape_value)) {
    #   if (verbose) {
    #     message("  Geometry is invalid. Attempting to fix...")
    #   }
    #   shape_value <- st_make_valid(shape_value)
    #   if (!st_is_valid(shape_value)) {
    #     if (verbose) {
    #       message("  Failed to fix geometry for speciesKey: ", sp_key)
    #     }
    #     return(NULL)  # Skip this speciesKey if the geometry remains invalid
    #   }
    #   if (verbose) {
    #     message("  Geometry successfully fixed.")
    #   }
    # } else if (verbose) {
    #   message("  Geometry is valid.")
    # }
    
    # Construct a 1-row data frame for the speciesKey
    row_df <- data.frame(
      speciesKey = sp_key,
      stringsAsFactors = FALSE
    )
    
    # Add other columns as lists of all values from the grouped rows
    for (colname in other_cols) {
      column_values <- subset_rows[[colname]]
      if (is.factor(column_values)) {
        column_values <- as.character(column_values)  # Convert factors to characters
      }
      row_df[[colname]] <- list(column_values)  # Store all values in a list
    }
    
    # Add the geometry directly into the data frame
    row_df[[shape_col]] <- shape_value
    
    return(row_df)
  }
  
  # Process each speciesKey in parallel using pbmclapply
  if (verbose) {
    message("Processing speciesKeys in parallel...")
  }
  collapsed_data <- pbmclapply(unique_keys, process_speciesKey, mc.cores = workers)
  
  # Remove NULL entries (skipped speciesKeys)
  collapsed_data <- Filter(Negate(is.null), collapsed_data)
  
  # Combine the list of one-row data frames into a single data frame
  if (verbose) {
    message("Combining collapsed data into a single data frame...")
  }
  result <- rbindlist(collapsed_data, ignore.attr = TRUE)
  
  # Final summary of geometry validation
  if (verbose) {
    message("Summary of geometry validation:")
    message("  Total geometries processed: ", length(unique_keys))
    #message("  Total skipped due to invalid geometries: ", length(unique_keys) - length(collapsed_data))
  }
  
  # Final message
  if (verbose) {
    message("Collapsing process complete. Returning consolidated data frame.")
  }
  
  return(result)
}


#' Assign H3 Hexagon Cells to a Polygon with Fallback
#'
#' Attempts to assign H3 cell indices to the given geometry at a specified
#' resolution. If normal polygon-to-cells mapping fails (e.g., very small or
#' degenerate geometry), falls back to using the bounding box centroid and
#' returns the cell(s) containing that point.
# Define a fallback function to ensure coverage
assign_h3_approx <- function(geometry, resolution) {
  # Try to assign H3 cells normally
  geometry_sf <- sf::st_sf(geometry = st_sfc(geometry), crs = 4326)
  h3_cells <- h3jsr::polygon_to_cells(geometry_sf, res = resolution, simple = TRUE)
  
  # Check if the result is NA (no H3 cell centers)
  if (all(is.na(h3_cells))) {
    # Use the bounding box centroid as a fallback
    centroid <- st_centroid(geometry)
    centroid_h3 <- h3jsr::point_to_cell(st_coordinates(centroid), res = resolution)
    
    # Expand to nearby H3 cells using get_disk
    h3_cells <- h3jsr::get_disk(centroid_h3, ring_size = 0)  # Adjust `ring_size` for desired coverage #here we set to 0 to avoid adding more cells
  }
  
  return(h3_cells)
}

#' Find Row Indices Containing a Given H3 Cell
#'
#' Scans a list-column of nested H3 cell vectors and returns the row indices
#' where a specific H3 cell appears (anywhere within the nested structure).
find_indices_for_h3_cell <- function(h3_cell, h3_cells_nested) {
  which(sapply(h3_cells_nested, function(cells) {
    h3_cell %in% unlist(cells, recursive = TRUE, use.names = FALSE)
  }))
}

#' Conditionally Log a Message
#'
#' Convenience helper to print a message only when \code{verbose = TRUE}.
log_message <- function(msg, verbose) {
  if (verbose) cat(msg, "\n")
}

#' Aggregate Species Rates and Richness per H3 Cell (Area-Weighted)
#'
#' Aggregates species-level metrics to H3 spatial cells using inverse-area
#' weights (1/area by default) computed from each species' total range area
#' (sum of \code{Shape_Area} across records with the same \code{speciesKey}).
#' Produces range-weighted and unweighted summaries for tip/lineage rates and
#' an H3-by-species weighting table.
aggregate_h3_metrics <- function(resident_species, verbose = TRUE, use_sqrt = FALSE, cores = 1) {
  log_message("Step 1: Validating input data...", verbose)
  
  required_cols <- c("h3_cells", "speciesKey", "log_tip_rate", "log_lineage_rate", "log_tipDR")
  missing_cols  <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'resident_species': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  has_verbatim_name <- "verbatim_name" %in% names(resident_species)
  if (!has_verbatim_name) {
    log_message("Warning: 'verbatim_name' column not found. Species names won't be included.", verbose)
  }
  
  row_lengths <- sapply(resident_species$h3_cells, length)
  empty_rows  <- which(row_lengths == 0)
  if (length(empty_rows) > 0) {
    cat("Rows with empty h3_cells (should not happen by definition):\n")
    print(empty_rows)
    stop("Encountered empty h3_cells, cannot proceed.")
  } else {
    log_message("No rows have empty h3_cells.", verbose)
  }
  
  log_message("Step 2: Flattening H3 cells and checking for extraneous cells...", verbose)
  flat_h3_cells <- unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE)
  actual_cells <- unique(flat_h3_cells)
  cat("Number of unique H3 cells in resident_species:", length(actual_cells), "\n")
  
  bad_cells <- setdiff(flat_h3_cells, actual_cells)
  if (length(bad_cells) > 0) {
    cat("WARNING: The following H3 cells are not actually used by any row in resident_species:\n")
    print(bad_cells)
    flat_h3_cells <- setdiff(flat_h3_cells, bad_cells)
    cat("Removed extraneous cells. Proceeding with filtered all_h3_cells.\n")
  } else {
    cat("No extraneous cells in all_h3_cells. Good to proceed.\n")
  }
  
  all_h3_cells <- unique(flat_h3_cells)
  cat("Final count of unique H3 cells to process:", length(all_h3_cells), "\n")
  
  log_message("Step 3: Aggregating metrics for each H3 cell...", verbose)
  
  #this old block computed range sizes as the count of the h3 cells
  # # Precompute range sizes for all species to avoid redundant calculations
  # all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
  #   sp_rows <- resident_species[resident_species$speciesKey == key, ]
  #   sp_cells <- unlist(sp_rows$h3_cells, recursive = TRUE, use.names = FALSE)
  #   length(unique(sp_cells))
  # })
  # names(all_range_sizes) <- unique(resident_species$speciesKey)
  
  #this updated block computes area as the actual areas of the shape files
  {
    # Precompute range sizes using Shape_Area for all species
    all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
      sp_rows <- resident_species[resident_species$speciesKey == key, ]
      # Sum up the Shape_Area values for the species
      sum(unlist(sp_rows$Shape_Area, recursive = TRUE, use.names = FALSE), na.rm = TRUE)
    })
    names(all_range_sizes) <- unique(resident_species$speciesKey)
    
    # Compute inverse weights for all species
    # all_range_weights <- 1 / all_range_sizes
    #all_range_sizes <- sqrt(all_range_sizes)
  }
  
  results_list <- pbmclapply(seq_along(all_h3_cells), function(i) {
    cell <- all_h3_cells[i]
    log_message(sprintf("Processing H3 cell %d of %d: %s", i, length(all_h3_cells), cell), verbose)
    
    indices <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
    if (length(indices) == 0) {
      stop(sprintf("Encountered an H3 cell (%s) with zero species. Should not happen!", cell))
    }
    
    cell_data <- resident_species[indices, ]
    unique_species_keys <- unique(cell_data$speciesKey)
    
    deduped_species_keys <- vector("character", length(unique_species_keys))
    deduped_tip_rates <- numeric(length(unique_species_keys))
    deduped_lineage_rates <- numeric(length(unique_species_keys))
    deduped_log_tipDR <- numeric(length(unique_species_keys)) 
    deduped_names <- vector("character", length(unique_species_keys))
    
    for (j in seq_along(unique_species_keys)) {
      this_key <- unique_species_keys[j]
      sub_idx <- which(cell_data$speciesKey == this_key)
      deduped_tip_rates[j] <- mean(cell_data$log_tip_rate[sub_idx], na.rm = TRUE)
      deduped_lineage_rates[j] <- mean(cell_data$log_lineage_rate[sub_idx], na.rm = TRUE)
      deduped_log_tipDR[j] <- mean(cell_data$log_tipDR[sub_idx], na.rm = TRUE)  # Extract mean log_tipDR
      deduped_species_keys[j] <- this_key
      if (has_verbatim_name) {
        these_names <- cell_data$verbatim_name[sub_idx]
        if (all(is.na(these_names))) {
          deduped_names[j] <- NA
        } else {
          deduped_names[j] <- these_names[!is.na(these_names)][1]
        }
      } else {
        deduped_names[j] <- NA
      }
    }
    
    species_keys <- deduped_species_keys
    tip_rates <- deduped_tip_rates
    lineage_rates <- deduped_lineage_rates
    species_names <- deduped_names
    
    range_sizes <- if (use_sqrt) sqrt(all_range_sizes[species_keys]) else all_range_sizes[species_keys]
    
    unnormalized_weights <- ifelse(range_sizes > 0, 1 / range_sizes, 0)
    sum_weights <- sum(unnormalized_weights, na.rm = TRUE)
    if (sum_weights == 0) {
      stop(sprintf("Sum of weights is zero for H3 cell=%s. Check data integrity.", cell))
    }
    normalized_weights <- unnormalized_weights / sum_weights
    tip_rate_range_weighted <- sum(tip_rates * normalized_weights, na.rm = TRUE)
    lineage_rate_range_weighted <- sum(lineage_rates * normalized_weights, na.rm = TRUE)
    log_tipDR_range_weighted <- sum(deduped_log_tipDR * normalized_weights, na.rm = TRUE)  # New weighted mean
    mean_tip_rate <- mean(tip_rates, na.rm = TRUE)
    lineage_rate <- mean(lineage_rates, na.rm = TRUE)
    mean_log_tipDR <- mean(deduped_log_tipDR, na.rm = TRUE)  # Compute mean log_tipDR here
    tip_rate_variance <- if (length(tip_rates) > 1) var(tip_rates, na.rm = TRUE) else NA
    lineage_rate_variance <- if (length(lineage_rates) > 1) var(lineage_rates, na.rm = TRUE) else NA
    log_tipDR_variance <- if (length(deduped_log_tipDR) > 1) var(deduped_log_tipDR, na.rm = TRUE) else NA  # Compute variance for log_tipDR
    
    # Compute range-weighted species count
    species_key_count_range_weighted <- sum(unnormalized_weights, na.rm = TRUE)
    
    log_message(sprintf(
      "  Number of deduplicated species: %d\n  Mean phenotypic rate: %.6f\n  Weighted mean phenotypic rate: %.6f\n  Phenotypic rate variance: %.6f\n  Mean lineage rate: %.6f\n  Weighted mean lineage rate: %.6f\n  Lineage rate variance: %.6f",
      length(species_keys), mean_tip_rate, tip_rate_range_weighted, tip_rate_variance, lineage_rate, lineage_rate_range_weighted, lineage_rate_variance
    ), verbose)
    
    weighting_details_df <- data.frame(
      h3_cell             = rep(cell, length(species_keys)),
      speciesKey          = species_keys,
      range_size          = range_sizes,
      weight_unnormalized = unnormalized_weights,
      weight_normalized   = normalized_weights
    )
    
    h3_metric <- list(
      h3_cell                 = cell, 
      mean_tip_rate           = mean_tip_rate, 
      tip_rate_range_weighted_mean = tip_rate_range_weighted, 
      tip_rate_variance       = tip_rate_variance, 
      mean_lineage_rate       = lineage_rate, 
      mean_lineage_rate_range_weighted = lineage_rate_range_weighted, 
      lineage_rate_variance   = lineage_rate_variance, 
      mean_tipDR              = mean_log_tipDR,  # Use precomputed mean log_tipDR
      tipDR_range_weighted_mean = log_tipDR_range_weighted,  # Weighted mean log_tipDR
      tipDR_variance          = log_tipDR_variance,  # Use precomputed variance log_tipDR
      species_list            = paste(unique(species_names), collapse = ", "), 
      species_key_list        = paste(unique(species_keys), collapse = ", "), 
      species_key_count       = length(species_keys), 
      species_key_count_range_weighted = species_key_count_range_weighted 
    )
    
    list(h3_metric = h3_metric, weighting_detail = weighting_details_df)
  }, mc.cores = cores)
  
  h3_metrics <- lapply(results_list, `[[`, "h3_metric")
  weighting_details_list <- lapply(results_list, `[[`, "weighting_detail")
  
  log_message("Step 4: Combining results into final outputs...", verbose)
  h3_metrics_df <- do.call(rbind, lapply(h3_metrics, as.data.frame))
  weighting_df <- do.call(rbind, weighting_details_list)
  
  if (verbose) {
    cat("\n===== Final h3_metrics_df structure =====\n")
    str(h3_metrics_df)
    cat("\n===== Sample of h3_metrics_df =====\n")
    print(head(h3_metrics_df))
    
    cat("\n===== Weighting details (weighting_df) =====\n")
    str(weighting_df)
    cat("\n===== Sample of weighting_df =====\n")
    print(head(weighting_df))
  }
  
  cat("\n===== Output Columns Description =====\n")
  cat("metrics_df:\n")
  cat("  - h3_cell: Identifier for the H3 spatial cell.\n")
  cat("  - mean_tip_rate: Calculated as the average of 'log_tip_rate' for all species in the H3 cell.\n")
  cat("  - tip_rate_range_weighted_mean: Calculated as the range-weighted mean of 'log_tip_rate' using species range sizes as weights.\n")
  cat("  - tip_rate_variance: Variance of 'log_tip_rate' for species in the H3 cell, or NA if only one species exists.\n")
  cat("  - mean_lineage_rate: Calculated as the average of 'log_lineage_rate' for all species in the H3 cell.\n")
  cat("  - mean_lineage_rate_range_weighted: Calculated as the range-weighted mean of 'log_lineage_rate' using species range sizes as weights.\n")
  cat("  - lineage_rate_variance: Variance of 'log_lineage_rate' for species in the H3 cell, or NA if only one species exists.\n")
  cat("  - mean_log_tipDR: Calculated as the average of 'log_tipDR' for all species in the H3 cell.\n")  # New metric
  cat("  - tipDR_range_weighted_mean: Calculated as the range-weighted mean of 'log_tipDR' using species range sizes as weights.\n")  # New metric
  cat("  - tipDR_variance: Variance of 'log_tipDR' for species in the H3 cell, or NA if only one species exists.\n")  # New metric
  cat("  - species_list: Comma-separated list of unique species names in the H3 cell.\n")
  cat("  - species_key_list: Comma-separated list of unique species identifiers (speciesKey) in the H3 cell.\n")
  cat("  - species_key_count: Total count of unique species in the H3 cell.\n")
  cat("  - species_key_count_range_weighted: Range-weighted species count for the H3 cell. This is the sum of inverse range size weights for all species in the cell, where smaller-range species contribute more to the total.\n")
  cat("\nweighting_df:\n")
  cat("  - h3_cell: Identifier for the H3 spatial cell.\n")
  cat("  - species_id: Unique identifier for each species in the H3 cell.\n")
  cat("  - species_range_size: Number of H3 cells the species spans, used for weighting calculations.\n")
  cat("  - unnormalized_weight: Raw weight calculated as 1 divided by the species range size (1 / species_range_size).\n")
  cat("  - normalized_weight: Weight normalized to sum to 1 across all species in the H3 cell.\n")
  
  return(list(
    metrics_df   = h3_metrics_df,
    weighting_df = weighting_df
  ))
}

#' Aggregate Range-Weighted Species Richness per H3 Cell (Counts Only)
#'
#' Computes H3-cell summaries focusing on richness, using inverse-area weights
#' derived from species range area (\code{Shape_Area}). Returns both the simple
#' species count and a range-weighted species count, plus a log-scale richness
#' metric.
aggregate_h3_metrics_counts <- function(resident_species, verbose = TRUE, use_sqrt = FALSE, cores = 1) {
  log_message("Step 1: Validating input data...", verbose)
  
  required_cols <- c("h3_cells", "speciesKey")
  missing_cols  <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from 'resident_species': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  has_verbatim_name <- "verbatim_name" %in% names(resident_species)
  if (!has_verbatim_name) {
    log_message("Warning: 'verbatim_name' column not found. Species names won't be included.", verbose)
  }
  
  row_lengths <- sapply(resident_species$h3_cells, length)
  empty_rows  <- which(row_lengths == 0)
  if (length(empty_rows) > 0) {
    cat("Rows with empty h3_cells (should not happen by definition):\n")
    print(empty_rows)
    stop("Encountered empty h3_cells, cannot proceed.")
  } else {
    log_message("No rows have empty h3_cells.", verbose)
  }
  
  log_message("Step 2: Flattening H3 cells and checking for extraneous cells...", verbose)
  flat_h3_cells <- unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE)
  actual_cells <- unique(flat_h3_cells)
  cat("Number of unique H3 cells in resident_species:", length(actual_cells), "\n")
  
  bad_cells <- setdiff(flat_h3_cells, actual_cells)
  if (length(bad_cells) > 0) {
    cat("WARNING: The following H3 cells are not actually used by any row in resident_species:\n")
    print(bad_cells)
    flat_h3_cells <- setdiff(flat_h3_cells, bad_cells)
    cat("Removed extraneous cells. Proceeding with filtered all_h3_cells.\n")
  } else {
    cat("No extraneous cells in all_h3_cells. Good to proceed.\n")
  }
  
  all_h3_cells <- unique(flat_h3_cells)
  cat("Final count of unique H3 cells to process:", length(all_h3_cells), "\n")
  
  log_message("Step 3: Aggregating metrics for each H3 cell...", verbose)
  
  #this old block computed range sizes as the count of the h3 cells
  # # Precompute range sizes for all species to avoid redundant calculations
  # all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
  #   sp_rows <- resident_species[resident_species$speciesKey == key, ]
  #   sp_cells <- unlist(sp_rows$h3_cells, recursive = TRUE, use.names = FALSE)
  #   length(unique(sp_cells))
  # })
  # names(all_range_sizes) <- unique(resident_species$speciesKey)
  
  #this updated block computes area as the actual areas of the shape files
  # Precompute range sizes using Shape_Area for all species
  all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
    sp_rows <- resident_species[resident_species$speciesKey == key, ]
    # Sum up the Shape_Area values for the species
    sum(unlist(sp_rows$Shape_Area, recursive = TRUE, use.names = FALSE), na.rm = TRUE)
  })
  names(all_range_sizes) <- unique(resident_species$speciesKey)
  
  results_list <- pbmclapply(seq_along(all_h3_cells), function(i) {
    cell <- all_h3_cells[i]
    log_message(sprintf("Processing H3 cell %d of %d: %s", i, length(all_h3_cells), cell), verbose)
    
    indices <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
    if (length(indices) == 0) {
      stop(sprintf("Encountered an H3 cell (%s) with zero species. Should not happen!", cell))
    }
    
    cell_data <- resident_species[indices, ]
    unique_species_keys <- unique(cell_data$speciesKey)
    
    deduped_species_keys <- vector("character", length(unique_species_keys))
    deduped_names <- vector("character", length(unique_species_keys)) 
    
    for (j in seq_along(unique_species_keys)) {
      this_key <- unique_species_keys[j]
      deduped_species_keys[j] <- this_key
      if (has_verbatim_name) {
        these_names <- cell_data$verbatim_name[cell_data$speciesKey == this_key]
        if (all(is.na(these_names))) {
          deduped_names[j] <- NA
        } else {
          deduped_names[j] <- these_names[!is.na(these_names)][1]
        }
      } else {
        deduped_names[j] <- NA
      }
    }
    
    species_keys <- deduped_species_keys
    species_names <- deduped_names
    
    #range_sizes <- all_range_sizes[species_keys]
    range_sizes <- if (use_sqrt) sqrt(all_range_sizes[species_keys]) else all_range_sizes[species_keys]
    log_weighted_values <- ifelse(range_sizes > 0, -log(range_sizes), 0)
    log_weighted_sum <- sum(log_weighted_values, na.rm = TRUE)
    
    unnormalized_weights <- ifelse(range_sizes > 0, 1 / range_sizes, 0)
    sum_weights <- sum(unnormalized_weights, na.rm = TRUE)
    if (sum_weights == 0) {
      stop(sprintf("Sum of weights is zero for H3 cell=%s. Check data integrity.", cell))
    }
    normalized_weights <- unnormalized_weights / sum_weights
    
    # New: Compute range-weighted species count
    species_key_count_range_weighted <- sum(unnormalized_weights, na.rm = TRUE)
    
    log_message(sprintf(
      "  Number of deduplicated species: %d",
      length(species_keys)
    ), verbose)
    
    weighting_details_df <- data.frame(
      h3_cell             = rep(cell, length(species_keys)),
      speciesKey          = species_keys,
      range_size          = range_sizes,
      weight_unnormalized = unnormalized_weights,
      weight_normalized   = normalized_weights
    )
    
    h3_metric <- list(
      h3_cell                 = cell, #Identifier for the H3 spatial cell.
      species_list            = paste(unique(species_names), collapse = ", "),
      species_key_list        = paste(unique(species_keys), collapse = ", "),
      species_key_count       = length(species_keys),
      species_key_count_range_weighted = species_key_count_range_weighted, # New metric
      log_range_weighted_species_richness = log_weighted_sum # 🚀 New metric added!
    )
    
    list(h3_metric = h3_metric, weighting_detail = weighting_details_df)
  }, mc.cores = cores)
  
  h3_metrics <- lapply(results_list, `[[`, "h3_metric")
  weighting_details_list <- lapply(results_list, `[[`, "weighting_detail")
  
  log_message("Step 4: Combining results into final outputs...", verbose)
  h3_metrics_df <- do.call(rbind, lapply(h3_metrics, as.data.frame))
  weighting_df <- do.call(rbind, weighting_details_list)
  
  if (verbose) {
    cat("\n===== Final h3_metrics_df structure =====\n")
    str(h3_metrics_df)
    cat("\n===== Sample of h3_metrics_df =====\n")
    print(head(h3_metrics_df))
    
    cat("\n===== Weighting details (weighting_df) =====\n")
    str(weighting_df)
    cat("\n===== Sample of weighting_df =====\n")
    print(head(weighting_df))
  }
  
  cat("\n===== Output Columns Description =====\n")
  cat("metrics_df:\n")
  cat("  - h3_cell: Identifier for the H3 spatial cell.\n")
  cat("  - species_list: Comma-separated list of unique species names in the H3 cell.\n")
  cat("  - species_key_list: Comma-separated list of unique species identifiers (speciesKey) in the H3 cell.\n")
  cat("  - species_key_count: Total count of unique species in the H3 cell.\n")
  cat("  - species_key_count_range_weighted: Range-weighted species count for the H3 cell. This is the sum of inverse range size weights for all species in the cell, where smaller-range species contribute more to the total.\n")
  cat("  - log_weighted_species_richness: **New metric!** The sum of log-transformed inverse range weights (`sum(log(1/area))`), providing a log-scale measure of range-weighted species richness.\n")
  
  cat("\nweighting_df:\n")
  cat("  - h3_cell: Identifier for the H3 spatial cell.\n")
  cat("  - species_id: Unique identifier for each species in the H3 cell.\n")
  cat("  - species_range_size: Number of H3 cells the species spans, used for weighting calculations.\n")
  cat("  - unnormalized_weight: Raw weight calculated as 1 divided by the species range size (1 / species_range_size).\n")
  cat("  - normalized_weight: Weight normalized to sum to 1 across all species in the H3 cell.\n")
  
  return(list(
    metrics_df   = h3_metrics_df,
    weighting_df = weighting_df
  ))
}


#' Aggregate Pairwise Mahalanobis Disparity per H3 Cell
#'
#' Calculates morphological disparity within each H3 spatial cell as the mean,
#' median, and range-weighted mean pairwise Mahalanobis distance among species'
#' residual trait vectors. Can use either standard covariance inversion
#' (\code{"solve"}) or shrinkage estimation (\code{"shrink"}) to handle
#' singular or ill-conditioned covariance matrices.
aggregate_h3_disparity.mahalanobis <- function(resident_species,
                                               verbose = TRUE,
                                               cores = 1,
                                               method_order = c("shrink", "solve"),
                                               plan_strategy = "multisession") {
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  
  if ("shrink" %in% method_order && !requireNamespace("corpcor", quietly = TRUE)) {
    stop("The 'corpcor' package is required for the 'shrink' method. Please install it.")
  }
  
  future::plan(strategy = plan_strategy, workers = cores)
  
  required_cols <- c("h3_cells", "speciesKey", "residuals", "Shape_Area")
  missing_cols <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  all_h3_cells <- unique(unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE))
  cat("Processing", length(all_h3_cells), "unique H3 cells\n")
  
  all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
    sp_rows <- resident_species[resident_species$speciesKey == key, ]
    sum(unlist(sp_rows$Shape_Area, recursive = TRUE, use.names = FALSE), na.rm = TRUE)
  })
  names(all_range_sizes) <- unique(resident_species$speciesKey)
  
  find_indices_for_h3_cell <- function(h3_cell, h3_cells_nested) {
    which(sapply(h3_cells_nested, function(cells) h3_cell %in% unlist(cells, recursive = TRUE, use.names = FALSE)))
  }
  
  result_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(all_h3_cells))
    
    future.apply::future_lapply(seq_along(all_h3_cells), function(i) {
      cell <- all_h3_cells[i]
      p(sprintf("Processing H3 cell %d", i))
      
      result <- tryCatch({
        indices <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
        cell_data <- resident_species[indices, ]
        unique_species_keys <- unique(cell_data$speciesKey)
        
        species_residuals <- lapply(unique_species_keys, function(key) {
          sp_rows <- cell_data[cell_data$speciesKey == key, ]
          if (length(sp_rows$residuals) == 0) return(NULL)
          as.numeric(sp_rows$residuals[[1]])
        })
        names(species_residuals) <- unique_species_keys
        species_residuals <- species_residuals[!sapply(species_residuals, is.null)]
        species_ids <- names(species_residuals)
        
        species_names <- if ("verbatim_name" %in% names(cell_data)) {
          sapply(species_ids, function(k) {
            nms <- cell_data$verbatim_name[cell_data$speciesKey == k]
            nms <- nms[!is.na(nms)]
            if (length(nms) > 0) nms[1] else NA
          })
        } else rep(NA, length(species_ids))
        
        fail_flag <- FALSE
        method_used <- NA
        mean_mahalanobis <- NA_real_
        mean_mahalanobis_weighted <- NA_real_
        dists <- NA
        weighted_dists <- NA
        pairnames <- NA
        
        if (length(species_residuals) >= 2) {
          residual_matrix <- do.call(rbind, species_residuals)
          cov_matrix <- cov(residual_matrix)
          inv_cov <- NULL
          
          for (method in method_order) {
            if (method == "solve") {
              inv_cov <- tryCatch(solve(cov_matrix), error = function(e) NULL)
            } else if (method == "shrink") {
              inv_cov <- tryCatch(corpcor::invcov.shrink(residual_matrix), error = function(e) NULL)
            }
            if (!is.null(inv_cov)) {
              method_used <- method
              break
            }
          }
          
          if (!is.null(inv_cov)) {
            n_species <- nrow(residual_matrix)
            n_pairs <- choose(n_species, 2)
            dists <- numeric(n_pairs)
            weighted_dists <- numeric(n_pairs)
            weights <- numeric(n_pairs)
            pairnames <- character(n_pairs)
            cell_ranges <- all_range_sizes[species_ids]
            
            ix <- 1
            for (j in 1:(n_species - 1)) {
              for (k in (j + 1):n_species) {
                diff_vec <- residual_matrix[j, ] - residual_matrix[k, ]
                dist_val <- sqrt(t(diff_vec) %*% inv_cov %*% diff_vec)
                dists[ix] <- dist_val
                r1 <- cell_ranges[j]
                r2 <- cell_ranges[k]
                weights[ix] <- if (r1 > 0 && r2 > 0) (1 / r1) * (1 / r2) else 0
                weighted_dists[ix] <- dist_val * weights[ix]
                pairnames[ix] <- paste(species_ids[j], species_ids[k], sep = "_")
                ix <- ix + 1
              }
            }
            
            mean_mahalanobis <- mean(dists, na.rm = TRUE)
            median_mahalanobis <- median(dists, na.rm = TRUE)
            if (sum(weights, na.rm = TRUE) > 0) {
              mean_mahalanobis_weighted <- sum(weighted_dists, na.rm = TRUE) / sum(weights, na.rm = TRUE)
            }
          } else {
            fail_flag <- TRUE
          }
        } else {
          fail_flag <- TRUE
        }
        
        list(
          metric = list(
            h3_cell = cell,
            mean_mahalanobis_disparity = mean_mahalanobis,
            mean_mahalanobis_disparity_weighted = mean_mahalanobis_weighted,
            median_mahalanobis_disparity = median_mahalanobis,
            species_key_count = length(species_ids),
            species_key_list = paste(species_ids, collapse = ", "),
            species_list = paste(unique(species_names), collapse = ", "),
            method_used = method_used
          ),
          failed = if (fail_flag) cell else NULL,
          distance_vector = if (!fail_flag) dists else NA,
          weighted_distance_vector = if (!fail_flag) weighted_dists else NA,
          pairwise_names = if (!fail_flag) pairnames else NA
        )
      }, error = function(e) {
        list(
          metric = list(
            h3_cell = cell,
            mean_mahalanobis_disparity = NA_real_,
            mean_mahalanobis_disparity_weighted = NA_real_,
            median_mahalanobis_disparity = NA_real_,
            species_key_count = 0,
            species_key_list = NA,
            species_list = NA,
            method_used = NA
          ),
          failed = cell,
          distance_vector = NA,
          weighted_distance_vector = NA,
          pairwise_names = NA
        )
      })
      
      return(result)
    }, future.stdout = TRUE, future.conditions = "condition", future.seed = TRUE)
  })
  
  metrics_df <- do.call(rbind, lapply(result_list, function(x) as.data.frame(x$metric)))
  failed_cells <- unlist(lapply(result_list, function(x) x$failed))
  
  distance_vectors <- setNames(
    lapply(result_list, function(x) x$distance_vector),
    sapply(result_list, function(x) x$metric$h3_cell)
  )
  
  weighted_distance_vectors <- setNames(
    lapply(result_list, function(x) x$weighted_distance_vector),
    sapply(result_list, function(x) x$metric$h3_cell)
  )
  
  pairwise_names <- setNames(
    lapply(result_list, function(x) x$pairwise_names),
    sapply(result_list, function(x) x$metric$h3_cell)
  )
  
  return(list(
    metrics_df = metrics_df,
    failed_cells = unique(failed_cells),
    distance_vectors = distance_vectors,
    weighted_distance_vectors = weighted_distance_vectors,
    pairwise_names = pairwise_names
  ))
}


#' Mahalanobis Dispersion to Centroid per H3 Cell
#'
#' Computes the average dispersion of species to the multivariate centroid
#' within each H3 spatial cell, using Mahalanobis distance based on species'
#' residual trait vectors. Also reports a range-weighted version that down-weights
#' wide-ranged species by inverse range area.
aggregate_h3_dispersion_to_centroid <- function(resident_species,
                                                verbose = TRUE,
                                                cores = 1,
                                                method_order = c("shrink", "solve"),
                                                plan_strategy = "multisession") {
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  
  if ("shrink" %in% method_order && !requireNamespace("corpcor", quietly = TRUE)) {
    stop("The 'corpcor' package is required for the 'shrink' method. Please install it.")
  }
  
  future::plan(strategy = plan_strategy, workers = cores)
  
  required_cols <- c("h3_cells", "speciesKey", "residuals", "Shape_Area")
  missing_cols <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  all_h3_cells <- unique(unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE))
  cat("Processing", length(all_h3_cells), "unique H3 cells\n")
  
  all_range_sizes <- sapply(unique(resident_species$speciesKey), function(key) {
    sp_rows <- resident_species[resident_species$speciesKey == key, ]
    sum(unlist(sp_rows$Shape_Area, recursive = TRUE, use.names = FALSE), na.rm = TRUE)
  })
  names(all_range_sizes) <- unique(resident_species$speciesKey)
  
  find_indices_for_h3_cell <- function(h3_cell, h3_cells_nested) {
    which(sapply(h3_cells_nested, function(cells) h3_cell %in% unlist(cells, recursive = TRUE, use.names = FALSE)))
  }
  
  result_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(all_h3_cells))
    
    future.apply::future_lapply(seq_along(all_h3_cells), function(i) {
      cell <- all_h3_cells[i]
      p(sprintf("Processing H3 cell %d", i))
      
      result <- tryCatch({
        indices <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
        cell_data <- resident_species[indices, ]
        unique_species_keys <- unique(cell_data$speciesKey)
        
        species_residuals <- lapply(unique_species_keys, function(key) {
          sp_rows <- cell_data[cell_data$speciesKey == key, ]
          if (length(sp_rows$residuals) == 0) return(NULL)
          as.numeric(sp_rows$residuals[[1]])
        })
        names(species_residuals) <- unique_species_keys
        species_residuals <- species_residuals[!sapply(species_residuals, is.null)]
        species_ids <- names(species_residuals)
        
        fail_flag <- FALSE
        method_used <- NA
        mean_distance_to_centroid <- NA_real_
        weighted_distance_to_centroid <- NA_real_
        
        if (length(species_residuals) >= 2) {
          residual_matrix <- do.call(rbind, species_residuals)
          centroid <- colMeans(residual_matrix, na.rm = TRUE)
          cov_matrix <- cov(residual_matrix)
          inv_cov <- NULL
          
          for (method in method_order) {
            if (method == "solve") {
              inv_cov <- tryCatch(solve(cov_matrix), error = function(e) NULL)
            } else if (method == "shrink") {
              inv_cov <- tryCatch(corpcor::invcov.shrink(residual_matrix), error = function(e) NULL)
            }
            if (!is.null(inv_cov)) {
              method_used <- method
              break
            }
          }
          
          if (!is.null(inv_cov)) {
            dists <- apply(residual_matrix, 1, function(vec) {
              diff_vec <- vec - centroid
              sqrt(t(diff_vec) %*% inv_cov %*% diff_vec)
            })
            mean_distance_to_centroid <- mean(dists, na.rm = TRUE)
            
            cell_ranges <- all_range_sizes[species_ids]
            weights <- ifelse(cell_ranges > 0, 1 / cell_ranges, 0)
            weighted_dists <- dists * weights
            weighted_distance_to_centroid <- if (sum(weights, na.rm = TRUE) > 0) {
              sum(weighted_dists, na.rm = TRUE) / sum(weights, na.rm = TRUE)
            } else {
              NA_real_
            }
          } else {
            fail_flag <- TRUE
          }
        } else {
          fail_flag <- TRUE
        }
        
        list(
          metric = data.frame(
            h3_cell = cell,
            mean_mahalanobis_distance_to_centroid = mean_distance_to_centroid,
            weighted_mahalanobis_distance_to_centroid = weighted_distance_to_centroid,
            species_key_count = length(species_ids),
            method_used = method_used
          ),
          failed = if (fail_flag) cell else NULL
        )
      }, error = function(e) {
        list(
          metric = data.frame(
            h3_cell = cell,
            mean_mahalanobis_distance_to_centroid = NA_real_,
            weighted_mahalanobis_distance_to_centroid = NA_real_,
            species_key_count = 0,
            method_used = NA
          ),
          failed = cell
        )
      })
      
      return(result)
    }, future.seed = TRUE)
  })
  
  metrics_df <- do.call(rbind, lapply(result_list, function(x) x$metric))
  failed_cells <- unlist(lapply(result_list, function(x) x$failed))
  
  return(list(
    metrics_df = metrics_df,
    failed_cells = unique(failed_cells)
  ))
}


#' Pairwise Phylogenetic Divergence Times per H3 Cell
#'
#' Extracts pairwise divergence times (in the same units as the phylogeny's branch lengths)
#' for all species occurring in each H3 spatial cell. Divergence times are computed as
#' the total tree height minus the MRCA height for each pair of species tips.
aggregate_h3_divergence_times <- function(resident_species,
                                          phy,
                                          key_to_tip_map,
                                          cores = 1,
                                          plan_strategy = "multisession") {
  requireNamespace("phytools")
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  
  if (!all(c("speciesKey", "phylo") %in% names(key_to_tip_map))) {
    stop("Mapping object must contain 'speciesKey' and 'phylo' columns.")
  }
  
  required_cols <- c("h3_cells", "speciesKey")
  missing_cols <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  future::plan(strategy = plan_strategy, workers = cores)
  
  all_h3_cells <- unique(unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE))
  cat("Processing", length(all_h3_cells), "unique H3 cells\n")
  
  mrca_matrix <- phytools::findMRCA(phy, type = "height")
  tree_height <- max(phytools::nodeHeights(phy)[, 2])
  tip_names <- phy$tip.label
  tip_index <- setNames(seq_along(tip_names), tip_names)
  
  find_indices_for_h3_cell <- function(h3_cell, h3_cells_nested) {
    which(sapply(h3_cells_nested, function(cells) h3_cell %in% unlist(cells, recursive = TRUE, use.names = FALSE)))
  }
  
  result_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(all_h3_cells))
    
    future.apply::future_lapply(seq_along(all_h3_cells), function(i) {
      cell <- all_h3_cells[i]
      p(sprintf("Processing H3 cell %d", i))
      
      result <- tryCatch({
        indices <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
        cell_data <- resident_species[indices, ]
        species_keys <- as.character(unique(cell_data$speciesKey))
        key_to_tip_map$speciesKey <- as.character(key_to_tip_map$speciesKey)
        
        matched_tips <- merge(
          data.frame(speciesKey = species_keys, stringsAsFactors = FALSE),
          key_to_tip_map,
          by = "speciesKey"
        )
        
        present_pairs <- matched_tips[matched_tips$phylo %in% tip_names, ]
        present_tips <- unique(present_pairs$phylo)
        fail_flag <- FALSE
        
        max_div <- NA_real_
        mean_div <- NA_real_
        median_div <- NA_real_
        mode_div <- NA_real_
        divergence_vec <- NA
        pairnames <- NA
        
        if (length(present_tips) >= 2) {
          tip_ids <- tip_index[present_tips]
          height_matrix <- mrca_matrix[tip_ids, tip_ids]
          div_matrix <- tree_height - height_matrix
          pairwise_divs <- div_matrix[upper.tri(div_matrix)]
          
          # Get speciesKeys in tip order for naming
          tip_order <- names(tip_ids)
          matched_ordered <- present_pairs[match(tip_order, present_pairs$phylo), ]
          species_ids <- matched_ordered$speciesKey
          
          pairnames <- character(length(pairwise_divs))
          ix <- 1
          for (j in 1:(length(species_ids) - 1)) {
            for (k in (j + 1):length(species_ids)) {
              pairnames[ix] <- paste(species_ids[j], species_ids[k], sep = "_")
              ix <- ix + 1
            }
          }
          
          divergence_vec <- pairwise_divs
          mean_div <- mean(pairwise_divs, na.rm = TRUE)
          max_div <- max(pairwise_divs, na.rm = TRUE)
          median_div <- median(pairwise_divs, na.rm = TRUE)
          dens <- density(pairwise_divs, na.rm = TRUE)
          mode_div <- dens$x[which.max(dens$y)]
        } else {
          fail_flag <- TRUE
          divergence_vec <- NA
          pairnames <- NA
        }
        
        list(
          metric = list(
            h3_cell = cell,
            max_divergence_time = max_div,
            mean_divergence_time = mean_div,
            median_divergence_time = median_div,
            mode_divergence_time = mode_div,
            n_species = length(present_tips),
            species_key_list = if (length(species_keys) > 0) paste(species_keys, collapse = ", ") else NA,
            phylo_tip_list = if (length(present_tips) > 0) paste(present_tips, collapse = ", ") else NA
          ),
          failed = if (fail_flag) cell else NULL,
          divergence_vector = divergence_vec,
          pairwise_names = pairnames
        )
      }, error = function(e) {
        list(
          metric = list(
            h3_cell = cell,
            max_divergence_time = NA_real_,
            mean_divergence_time = NA_real_,
            median_divergence_time = NA_real_,
            mode_divergence_time = NA_real_,
            n_species = 0,
            species_key_list = NA,
            phylo_tip_list = NA
          ),
          failed = cell,
          divergence_vector = NA,
          pairwise_names = NA
        )
      })
      
      return(result)
    }, future.stdout = TRUE, future.conditions = "condition", future.seed = TRUE)
  })
  
  metrics_df <- do.call(rbind, lapply(result_list, function(x) as.data.frame(x$metric)))
  failed_cells <- unlist(lapply(result_list, function(x) x$failed))
  
  divergence_vectors <- setNames(
    lapply(result_list, function(x) x$divergence_vector),
    sapply(result_list, function(x) x$metric$h3_cell)
  )
  pairwise_names <- setNames(
    lapply(result_list, function(x) x$pairwise_names),
    sapply(result_list, function(x) x$metric$h3_cell)
  )
  
  return(list(
    metrics_df = metrics_df,
    failed_cells = unique(failed_cells),
    divergence_vectors = divergence_vectors,
    pairwise_names = pairwise_names
  ))
}


#' Summarize Disparity–Divergence Relationships via Regression Through the Origin
#'
#' Fits a no-intercept regression of pairwise morphological disparity
#' (Mahalanobis distances) against pairwise phylogenetic divergence times
#' for each H3 grid cell, and summarizes residual statistics as a measure of
#' “excess” disparity beyond that expected from divergence time.
summarize_disparity_vs_divergence <- function(
    mahalanobis_data,
    divergence_data,
    use_weighted   = FALSE,
    use_median     = FALSE,
    bootstrap      = FALSE,
    n_boot         = 1000,
    min_pairs      = 5,
    standardize    = FALSE,
    keep_model_fit = FALSE,
    cores          = 1,
    plan_strategy  = "multisession"
) {
  requireNamespace("fastmatch")
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  requireNamespace("e1071")  # for skewness
  
  # parallel plan
  future::plan(strategy = plan_strategy, workers = cores)
  
  distance_key <- if (use_weighted) "weighted_distance_vectors" else "distance_vectors"
  if (!distance_key %in% names(mahalanobis_data)) {
    stop("Missing expected key in mahalanobis_data: ", distance_key)
  }
  
  common_cells <- intersect(
    names(mahalanobis_data[[distance_key]]),
    names(divergence_data$divergence_vectors)
  )
  if (length(common_cells) == 0) {
    stop("No grid cells in common between distance and divergence data")
  }
  
  # setup a text progress bar, one tick per cell
  progressr::handlers("txtprogressbar")
  p <- progressr::progressor(steps = length(common_cells))
  
  result_list <- progressr::with_progress({
    future.apply::future_lapply(common_cells, function(cell) {
      p()  # advance
      
      # extract matched distance & divergence vectors
      dist_vec   <- mahalanobis_data[[distance_key]][[cell]]
      div_vec    <- divergence_data$divergence_vectors[[cell]]
      dist_names <- mahalanobis_data$pairwise_names[[cell]]
      div_names  <- divergence_data$pairwise_names[[cell]]
      
      if (is.null(dist_vec) || is.null(div_vec) ||
          is.null(dist_names) || is.null(div_names)) {
        return(NULL)
      }
      
      shared_pairs <- intersect(dist_names, div_names)
      if (length(shared_pairs) < min_pairs) return(NULL)
      
      # index matching
      dist_idx <- fastmatch::fmatch(shared_pairs, dist_names)
      div_idx  <- fastmatch::fmatch(shared_pairs, div_names)
      
      dist_matched <- dist_vec[dist_idx]
      div_matched  <- div_vec[div_idx]
      
      # ---- FIXED HERE: removed extra ')' ----
      if (anyNA(dist_matched) || anyNA(div_matched)) return(NULL)
      
      if (standardize) {
        dist_matched <- scale(dist_matched)[, 1]
        div_matched  <- scale(div_matched)[, 1]
      }
      
      # fit no-intercept OLS
      initial_model <- lm(dist_matched ~ div_matched -1)
      resids <- residuals(initial_model)
      
      mean_resid    <- mean(resids, na.rm = TRUE)
      sd_resid      <- sd(resids, na.rm = TRUE)
      std_resid     <- mean_resid / sd_resid
      prop_pos      <- mean(resids > 0, na.rm = TRUE)
      median_resid  <- median(resids, na.rm = TRUE)
      skew_resid    <- e1071::skewness(resids, na.rm = TRUE, type = 2)
      prop_extreme  <- mean(abs(resids) > 2 * sd_resid, na.rm = TRUE)
      prop_ext_pos  <- mean(resids > 2 * sd_resid, na.rm = TRUE)
      prop_ext_neg  <- mean(resids < -2 * sd_resid, na.rm = TRUE)
      
      if (bootstrap) {
        mean_resids <- numeric(n_boot)
        std_resids  <- numeric(n_boot)
        all_boot    <- vector("list", n_boot)
        
        for (i in seq_len(n_boot)) {
          ix       <- sample(seq_along(resids), replace = TRUE)
          boot_fit <- lm(dist_matched[ix] ~ div_matched[ix] -1)
          boot_res <- residuals(boot_fit)
          mean_resids[i] <- mean(boot_res, na.rm = TRUE)
          std_resids[i]  <- mean_resids[i] / sd(boot_res, na.rm = TRUE)
          all_boot[[i]]  <- boot_res
        }
        
        all_boot_resids <- unlist(all_boot)
        resid_ci        <- quantile(mean_resids, c(0.025, 0.975), na.rm = TRUE)
        std_ci          <- quantile(std_resids,  c(0.025, 0.975), na.rm = TRUE)
        
        mean_pooled     <- mean(all_boot_resids, na.rm = TRUE)
        sd_pooled       <- sd(all_boot_resids, na.rm = TRUE)
        std_pooled      <- mean_pooled / sd_pooled
        
        result <- list(
          h3_cell                  = cell,
          mean_residual            = mean(mean_resids, na.rm = TRUE),
          resid_ci_lower           = resid_ci[1],
          resid_ci_upper           = resid_ci[2],
          standardized_residual    = if (standardize) NA else mean(std_resids, na.rm = TRUE),
          std_resid_ci_lower       = if (standardize) NA else std_ci[1],
          std_resid_ci_upper       = if (standardize) NA else std_ci[2],
          mean_pooled_residual     = mean_pooled,
          sd_pooled_residual       = sd_pooled,
          standardized_pooled_res  = std_pooled,
          prop_positive_residuals  = mean(all_boot_resids > 0, na.rm = TRUE),
          median_residual          = median(all_boot_resids, na.rm = TRUE),
          skewness_residual        = e1071::skewness(all_boot_resids, na.rm = TRUE, type = 2),
          prop_extreme_residuals   = mean(abs(all_boot_resids) > 2 * sd_pooled, na.rm = TRUE),
          prop_extreme_pos         = mean(all_boot_resids > 2 * sd_pooled, na.rm = TRUE),
          prop_extreme_neg         = mean(all_boot_resids < -2 * sd_pooled, na.rm = TRUE),
          n_pairs                  = length(shared_pairs),
          weighted                 = use_weighted,
          standardized             = standardize
        )
      } else {
        result <- list(
          h3_cell                 = cell,
          mean_residual           = mean_resid,
          standardized_residual   = if (standardize) NA else std_resid,
          prop_positive_residuals = prop_pos,
          median_residual         = median_resid,
          skewness_residual       = skew_resid,
          prop_extreme_residuals  = prop_extreme,
          prop_extreme_pos        = prop_ext_pos,
          prop_extreme_neg        = prop_ext_neg,
          n_pairs                 = length(shared_pairs),
          weighted                = use_weighted,
          standardized            = standardize
        )
      }
      
      if (keep_model_fit) result$model_fit <- initial_model
      return(result)
    }, future.seed = TRUE)
  })
  
  # drop NULLs, assemble summary_df & model_fits
  clean   <- Filter(Negate(is.null), result_list)
  summary <- do.call(rbind, lapply(clean, as.data.frame))
  fits    <- if (keep_model_fit) setNames(lapply(clean, `[[`, "model_fit"), summary$h3_cell) else NULL
  
  future::plan(sequential)
  return(list(summary_df = summary, model_fits = fits))
}


#' Calculate Multivariate Phylogenetic Signal Metrics for Each H3 Grid Cell
#'
#' Computes phylogenetic signal statistics (Z, K, and Pagel's lambda) for
#' multivariate trait data within each H3 grid cell using \code{\link[geomorph]{physignal.z}}.
aggregate_phylo_signal <- function(
    resident_species,
    phy,
    key_to_tip_map,
    trait_column = "residuals",
    min_species = 3,
    iter = 0,
    lambda_method = "front",
    cores = 1,
    plan_strategy = "multisession",
    ...
) {
  requireNamespace("ape")
  requireNamespace("geomorph")
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  
  future::plan(strategy = plan_strategy, workers = cores)
  
  key_to_tip_map$speciesKey <- as.character(key_to_tip_map$speciesKey)
  resident_species$speciesKey <- as.character(resident_species$speciesKey)
  
  all_h3_cells <- unique(unlist(resident_species$h3_cells, recursive = TRUE, use.names = FALSE))
  cat("Computing multivariate phylogenetic signal (Z, K, lambda) for", length(all_h3_cells), "H3 cells\n")
  
  find_indices_for_h3_cell <- function(h3_cell, h3_cells_nested) {
    which(sapply(h3_cells_nested, function(cells) h3_cell %in% unlist(cells, recursive = TRUE, use.names = FALSE)))
  }
  
  result_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(all_h3_cells))
    
    future.apply::future_lapply(seq_along(all_h3_cells), function(i) {
      cell <- all_h3_cells[i]
      
      tryCatch({
        idx <- find_indices_for_h3_cell(cell, resident_species$h3_cells)
        cell_data <- resident_species[idx, ]
        species_keys <- unique(as.character(cell_data$speciesKey))
        
        matched_tips <- merge(
          data.frame(speciesKey = species_keys, stringsAsFactors = FALSE),
          key_to_tip_map,
          by = "speciesKey"
        )
        
        tip_names <- phy$tip.label
        matched_tips <- matched_tips[matched_tips$phylo %in% tip_names, ]
        present_tips <- unique(matched_tips$phylo)
        
        if (length(present_tips) < min_species) {
          p()
          return(data.frame(h3_cell = cell, Z = NA, K = NA, lambda = NA, p_value = NA, n_species = length(present_tips)))
        }
        
        tip_tree <- ape::drop.tip(phy, setdiff(tip_names, present_tips))
        matched_cell <- merge(cell_data, matched_tips, by = "speciesKey")
        
        trait_matrix_list <- lapply(tip_tree$tip.label, function(tip) {
          vals <- matched_cell[[trait_column]][matched_cell$phylo == tip]
          if (length(vals) == 0) return(NULL)
          mat <- matrix(unlist(vals), nrow = 1)
          rownames(mat) <- tip
          return(mat)
        })
        trait_matrix_list <- Filter(Negate(is.null), trait_matrix_list)
        
        if (length(trait_matrix_list) < min_species) {
          p()
          return(data.frame(h3_cell = cell, Z = NA, K = NA, lambda = NA, p_value = NA, n_species = length(trait_matrix_list)))
        }
        
        trait_matrix <- do.call(rbind, trait_matrix_list)
        
        # Match rows of trait_matrix exactly to tree tip labels
        trait_matrix <- trait_matrix[tip_tree$tip.label, , drop = FALSE]
        
        if (nrow(trait_matrix) != ape::Ntip(tip_tree)) {
          p()
          return(data.frame(h3_cell = cell, Z = NA, K = NA, lambda = NA, p_value = NA, n_species = nrow(trait_matrix)))
        }
        
        signal_z <- geomorph::physignal.z(
          A = trait_matrix,
          phy = tip_tree,
          lambda = lambda_method,
          iter = iter,
          ...
        )
        
        p()
        data.frame(
          h3_cell = cell,
          Z = signal_z$Z,
          K = signal_z$K,
          lambda = signal_z$lambda,
          p_value = signal_z$pvalue,
          n_species = nrow(trait_matrix)
        )
        
      }, error = function(e) {
        message(sprintf("Error in cell %s: %s", cell, e$message))
        p()
        return(data.frame(h3_cell = cell, Z = NA, K = NA, lambda = NA, p_value = NA, n_species = 0))
      })
    }, future.seed = TRUE)
  })
  
  final_df <- do.call(rbind, result_list)
  return(final_df)
}


#' Prepare latitude-annotated data frames for plotting
#'
#' Combines disparity results (required) with optional phylogenetic-signal (K/Z)
#' and dispersion results, and appends latitude (from each H3 cell center).
prepare_plot_data <- function(disparity_results,
                              kz_results = NULL,
                              add_kz_metrics = FALSE,
                              dispersion_results = NULL,
                              add_dispersion_metrics = FALSE) {
  df1 <- disparity_results$summary_df
  df1$h3_cell <- as.character(df1$h3_cell)
  df1$latitude <- sf::st_coordinates(h3jsr::cell_to_point(df1$h3_cell))[,2]
  
  df2 <- NULL
  if (add_kz_metrics && !is.null(kz_results)) {
    df2 <- kz_results
    df2$h3_cell <- as.character(df2$h3_cell)
    df2$latitude <- sf::st_coordinates(h3jsr::cell_to_point(df2$h3_cell))[,2]
  }
  
  df3 <- NULL
  if (add_dispersion_metrics && !is.null(dispersion_results)) {
    if (is.list(dispersion_results) && "metrics_df" %in% names(dispersion_results)) {
      df3 <- dispersion_results$metrics_df
    } else {
      df3 <- dispersion_results
    }
    df3$h3_cell <- as.character(df3$h3_cell)
    df3$latitude <- sf::st_coordinates(h3jsr::cell_to_point(df3$h3_cell))[,2]
  }
  
  list(df1 = df1, df2 = df2, df3 = df3)
}

#' Tag polygons with continent/region labels (planar unions)
#'
#' Intersects input polygons (e.g., H3 cells) with continent boundaries built
#' from \pkg{rnaturalearth} country polygons (unioned in planar mode) and adds a
#' \code{region} factor column. Central America is labeled explicitly.
prepare_continent_overlap <- function(polygons_sf) {
  # ensure we’re in WGS84 (as rnaturalearth output is)
  polygons_sf <- sf::st_transform(polygons_sf, crs = 4326)
  
  # Turn off s2 spherical to avoid topology errors
  sf::sf_use_s2(FALSE)
  
  # Define continent names
  continent_names <- c(
    "Africa", "Antarctica", "Asia", "Europe",
    "North America", "South America", "Oceania"
  )
  
  # Build planar union for each continent
  continent_polygons <- lapply(continent_names, function(cont) {
    # download countries
    countries <- rnaturalearth::ne_countries(
      continent = cont, returnclass = "sf"
    )
    # validify and union
    countries <- sf::st_make_valid(countries)
    sf::st_union(countries)
  })
  names(continent_polygons) <- continent_names
  
  # Manually union Central America
  central_countries <- c(
    "Belize", "Costa Rica", "El Salvador", "Guatemala",
    "Honduras", "Nicaragua", "Panama"
  )
  ca <- rnaturalearth::ne_countries(
    country = central_countries, returnclass = "sf"
  )
  central_poly <- sf::st_make_valid(ca) %>% sf::st_union()
  
  # Prepare output
  polygons_sf$region <- "Other"
  
  # Do intersects (planar)
  for (cont in continent_names) {
    hit <- sf::st_intersects(
      polygons_sf,
      continent_polygons[[cont]],
      sparse = FALSE
    )[,1]
    polygons_sf$region[hit] <- cont
  }
  
  # Override Central America
  hit_ca <- sf::st_intersects(polygons_sf, central_poly, sparse = FALSE)[,1]
  polygons_sf$region[hit_ca] <- "Central America"
  
  # Rename Oceania
  polygons_sf$region[polygons_sf$region == "Oceania"] <- "Australia/Oceania"
  
  # Factor with desired order
  polygons_sf$region <- factor(
    polygons_sf$region,
    levels = c(
      "Africa", "Antarctica", "Asia", "Europe",
      "North America", "Central America", "South America",
      "Australia/Oceania", "Other"
    )
  )
  
  # (Re‐)enable S2 if you like for other operations
  # sf::sf_use_s2(TRUE)
  
  polygons_sf
}

#' Build a latitude–metric plot with optional regional coloring and summaries
#'
#' Creates a \pkg{ggplot2} scatter of a chosen metric vs. latitude, optionally
#' coloring points by continent/region and adding reference lines, LOESS fits,
#' and sliding-window proportions (for K significance).
make_metric_plot <- function(df_plot,
                             metric,
                             polygons_sf,
                             point_size,
                             point_alpha,
                             add_loess,
                             zero_line_metrics,
                             one_line_metrics,
                             p_threshold,
                             p_adjust_method,
                             bin_size,
                             step_size,
                             max_lat,
                             color_by_region = FALSE,
                             legend_alpha = 1) {
  # --- setup -------------------------------------------------------
  # labels & palette
  metric_labels <- c(
    mean_residual = "Mean Residual",
    standardized_residual = "Standardized Residual",
    prop_positive_residuals = "Proportion Positive Residuals",
    median_residual = "Median Residual",
    skewness_residual = "Skewness of Residuals",
    prop_extreme_residuals = "Proportion Extreme Residuals",
    prop_extreme_pos = "Proportion Extreme Positive",
    prop_extreme_neg = "Proportion Extreme Negative",
    K = expression(K["mult"]),
    Z = "Z-Score (Phylogenetic Signal)",
    lambda = "Pagel's Lambda",
    p_value = "Phylogenetic Signal p-value",
    mean_mahalanobis_disparity   = "Mean Pairwise\nMahalanobis Distance",
    median_mahalanobis_disparity = "Median Pairwise\nMahalanobis Distance"
  )
  continent_colors <- c(
    "Africa"            = "#E76F00", "Antarctica"        = "#6BC4FF",
    "Asia"              = "#1FA87C", "Europe"            = "#4858D9",
    "North America"     = "#0085CA", "Central America"   = "#FFCC29",
    "South America"     = "#C72E29", "Australia/Oceania" = "#D845A6",
    "Other"             = "#9E9E9E"
  )
  
  # -- filter & join region ----------------------------------------
  df_metric <- df_plot %>%
    dplyr::filter(!is.na(.data[[metric]]))
  
  if (color_by_region || metric == "standardized_residual") {
    df_metric <- df_metric %>%
      dplyr::left_join(
        polygons_sf %>%
          sf::st_drop_geometry() %>%
          dplyr::select(h3_address, region),
        by = c("h3_cell" = "h3_address")
      )
  }
  
  # -- handle K-significance ---------------------------------------
  if (metric == "K" && "p_value" %in% names(df_metric)) {
    df_metric$significant <- p.adjust(df_metric$p_value, p_adjust_method) < p_threshold
  }
  
  
  # --- build aesthetic mapping -----------------------------------
  aes_map <- aes(x = latitude, y = .data[[metric]])
  if (metric == "K" && !color_by_region && "significant" %in% names(df_metric)) {
    aes_map <- modifyList(aes_map, aes(color = significant))
    scale_fun <- scale_color_manual(values = c("TRUE" = "black", "FALSE" = "darkgrey"))
    guide_fun <- guides(color = "none")
  } else if (color_by_region || metric == "standardized_residual") {
    aes_map <- modifyList(aes_map, aes(color = region))
    scale_fun <- scale_color_manual(values = continent_colors)
    guide_fun <- guides(
      color = guide_legend(
        title = "Region",
        override.aes = list(size = 4, alpha = legend_alpha)
      )
    )
  } else {
    scale_fun <- NULL
    guide_fun <- NULL
  }
  
  # --- start plot -------------------------------------------------
  p <- ggplot(df_metric, aes_map) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = point_size, alpha = point_alpha)
  
  # add color scales/guides
  if (!is.null(scale_fun)) p <- p + scale_fun
  if (!is.null(guide_fun)) p <- p + guide_fun
  
  # --- sliding window for K ---------------------------------------
  if (metric == "K" && "significant" %in% names(df_metric)) {
    bins       <- seq(-90, 90 - bin_size, by = step_size)
    prop_data  <- data.frame(
      bin_center = bins + bin_size / 2,
      prop_sig   = sapply(bins, function(b) {
        w <- df_metric$latitude >= b & df_metric$latitude < b + bin_size
        if (any(w)) mean(df_metric$significant[w], na.rm = TRUE) else NA
      })
    ) %>% tidyr::drop_na()
    
    p <- p + geom_line(
      data = prop_data,
      aes(x = bin_center, y = prop_sig),
      inherit.aes = FALSE,
      linewidth = 1, color = alpha("black", 0.75)
    )
  }
  
  # --- common layers ---------------------------------------------
  if (add_loess) {
    p <- p + geom_smooth(method = "loess", se = TRUE)
  }
  if (metric %in% zero_line_metrics) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  }
  if (metric %in% one_line_metrics) {
    p <- p + geom_hline(yintercept = 1, linetype = "dotted", color = "red")
  }
  
  p +
    scale_x_continuous(limits = c(-max_lat, max_lat), oob = scales::oob_squish) +
    scale_y_continuous(oob = scales::oob_squish) +
    labs(
      title = if (metric == "K") bquote(bold(K["mult"])) else (metric_labels[[metric]] %||% metric),
      x = "Latitude",
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_blank(),
      axis.line    = element_line(color = "black", linewidth = 0.7),
      axis.ticks   = element_line(color = "black", linewidth = 0.4),
      axis.text    = element_text(color = "black"),
      plot.title   = element_text(hjust = 0, face = "bold", size = 13),
      axis.text.x  = element_text(size = 11),
      plot.margin  = margin(5, 10, 5, 10)
    )
}

#' Plot multiple disparity / signal / dispersion metrics by latitude
#'
#' Wrapper that prepares inputs, attaches continent/region labels, and generates
#' a grid of latitude–metric plots using \pkg{patchwork}.
plot_disparity_metrics_by_latitude <- function(
    disparity_results,
    metrics = c(
      "mean_residual", "standardized_residual", "prop_positive_residuals", "median_residual", 
      "skewness_residual", "prop_extreme_residuals", "prop_extreme_pos", "prop_extreme_neg",
      "K", "Z",
      "mean_mahalanobis_disparity",    
      "median_mahalanobis_disparity"),
    kz_results = NULL,
    add_kz_metrics = FALSE,
    dispersion_results = NULL,
    add_dispersion_metrics = FALSE,
    zero_line_metrics = c("standardized_residual", "Z"),
    one_line_metrics  = c("K", "standardized_residual", "prop_positive_residuals"),
    point_size = 1.6,
    point_alpha = 0.6,
    add_loess = FALSE,
    bin_size = 5,
    step_size = 1,
    p_threshold = 0.05,
    p_adjust_method = "none",
    color_by_region = FALSE,
    legend_alpha = 1
) {
  library(dplyr); library(ggplot2); library(h3jsr)
  library(sf);   library(patchwork); library(scales)
  
  if (missing(add_kz_metrics))        add_kz_metrics        <- !is.null(kz_results)
  if (missing(add_dispersion_metrics)) add_dispersion_metrics <- !is.null(dispersion_results)
  
  # 1) Prepare data frames
  data_list <- prepare_plot_data(
    disparity_results, kz_results, add_kz_metrics,
    dispersion_results, add_dispersion_metrics
  )
  df1 <- data_list$df1; df2 <- data_list$df2; df3 <- data_list$df3
  
  # 2) Select metrics that actually exist
  all_names     <- union(names(df1), union(names(df2), names(df3)))
  valid_metrics <- intersect(metrics, all_names)
  if (!length(valid_metrics)) stop("No valid metrics to plot.")
  
  # 3) Compute axis limit
  max_lat <- ceiling(max(abs(c(
    df1$latitude,
    if (!is.null(df2)) df2$latitude,
    if (!is.null(df3)) df3$latitude
  )), na.rm = TRUE))
  
  # 4) **Always** rebuild region labels on your external polygons_sf.1
  polygons_for_plot <- prepare_continent_overlap(polygons_sf.1)
  
  # 5) Build each metric’s plot
  plots <- lapply(valid_metrics, function(metric) {
    df_plot <- if (!is.null(df2) && metric %in% names(df2)) {
      df2
    } else if (!is.null(df3) && metric %in% names(df3)) {
      df3
    } else {
      df1
    }
    make_metric_plot(
      df_plot, metric, polygons_for_plot,
      point_size, point_alpha, add_loess,
      zero_line_metrics, one_line_metrics,
      p_threshold, p_adjust_method,
      bin_size, step_size, max_lat,
      color_by_region = color_by_region, 
      legend_alpha = legend_alpha
    )
  })
  
  # 6) Combine and return
  patchwork::wrap_plots(plots, ncol = 3)
}


#' Estimate Local Trait Covariance Structure and Dimensionality
#'
#' Computes phylogenetically corrected and empirical trait covariance matrices 
#' for species within each H3 grid cell, and extracts multiple metrics describing 
#' the scale, shape, and evenness of trait space usage.
local_vcv <- function(resident_species,
                      phy,
                      key_to_tip_map,
                      n_min_species = 3,
                      verbose = TRUE,
                      plan_strategy = "multisession",
                      cores = 1) {
  
  # ---- Metric interpretation ----
  # det_Pinv:     Generalized variance of trait space (volume); sensitive to both scale and shape.
  # trace:        Total variance across all trait dimensions.
  # trace_per_species: Trace normalized by species richness in the cell.
  # lambda1:      Largest eigenvalue (dominant axis of trait variance).
  # lambda1_prop: Proportion of total variance on the leading axis.
  # condition_number: Ratio of max to min eigenvalue.
  # effective_dimensionality: exp(entropy) — how many trait axes are used.
  # evenness:     Normalized entropy (0–1); describes variance evenness.
  # effective_dim_per_species: Dimensionality per species — richness-normalized trait axis usage.
  # modularity:   Community structure in trait covariance (via Louvain clustering).
  
  requireNamespace("mvMORPH")
  requireNamespace("future")
  requireNamespace("future.apply")
  requireNamespace("progressr")
  requireNamespace("igraph")
  
  required_cols <- c("h3_cells", "speciesKey", "residuals")
  missing_cols  <- setdiff(required_cols, names(resident_species))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  all_h3_cells <- unique(unlist(resident_species$h3_cells, recursive = TRUE))
  key_to_tip_map$speciesKey <- as.character(key_to_tip_map$speciesKey)
  resident_species$speciesKey <- as.character(resident_species$speciesKey)
  
  future::plan(strategy = plan_strategy, workers = cores)
  
  find_idx <- function(cell, nested) {
    which(sapply(nested, function(v) cell %in% unlist(v, recursive = TRUE)))
  }
  
  result <- progressr::with_progress({
    p <- progressr::progressor(steps = length(all_h3_cells))
    
    future.apply::future_lapply(all_h3_cells, function(cell) {
      tryCatch({
        idx <- find_idx(cell, resident_species$h3_cells)
        cell_data <- resident_species[idx, ]
        species_k <- unique(cell_data$speciesKey)
        map_df <- key_to_tip_map[key_to_tip_map$speciesKey %in% species_k, ]
        tips <- intersect(map_df$phylo, phy$tip.label)
        
        if (length(tips) < n_min_species) return(list(h3_cell = cell, success = FALSE))
        
        sub_tree <- ape::drop.tip(phy, setdiff(phy$tip.label, tips))
        cell_data <- merge(cell_data, map_df, by = "speciesKey")
        
        tm_list <- lapply(sub_tree$tip.label, function(tp) {
          v <- cell_data$residuals[cell_data$phylo == tp]
          if (length(v) == 0) return(NULL)
          m <- matrix(unlist(v), nrow = 1); rownames(m) <- tp; m
        })
        tm_list <- Filter(Negate(is.null), tm_list)
        if (length(tm_list) < n_min_species) return(list(h3_cell = cell, success = FALSE))
        
        X <- do.call(rbind, tm_list)
        
        ## -- Empirical covariance metrics --
        empirical_cov <- cov(X)
        eig_emp <- eigen(empirical_cov, symmetric = TRUE)$values
        eig_emp <- eig_emp[eig_emp > 0]
        
        cor_emp <- abs(cor(X))
        graph_emp <- igraph::graph_from_adjacency_matrix((cor_emp + t(cor_emp)) / 2,
                                                         mode = "undirected", weighted = TRUE, diag = FALSE)
        modularity_emp <- igraph::modularity(igraph::cluster_louvain(graph_emp))
        #clust_emp <- igraph::cluster_leiden(graph_emp, objective_function = "modularity")
        #modularity_emp <- igraph::modularity(clust_emp)
        
        emp_metrics <- list(
          det = det(empirical_cov),
          trace = sum(eig_emp),
          trace_per_species = sum(eig_emp) / nrow(X),
          lambda1 = max(eig_emp),
          lambda1_prop = max(eig_emp) / sum(eig_emp),
          condition_number = log(max(eig_emp) / min(eig_emp)),
          effective_dim = exp(-sum((eig_emp / sum(eig_emp)) * log(eig_emp / sum(eig_emp) + 1e-12))),
          evenness = -sum((eig_emp / sum(eig_emp)) * log(eig_emp / sum(eig_emp) + 1e-12)) / log(length(eig_emp)),
          effective_dim_per_species = exp(-sum((eig_emp / sum(eig_emp)) * log(eig_emp / sum(eig_emp) + 1e-12))) / nrow(X),
          modularity = modularity_emp
        )
        
        ## -- Model-based covariance (mvGLS) --
        fit_gls <- tryCatch({
          mvMORPH::mvgls(Y ~ 1, data = list(Y = X), tree = sub_tree,
                         model = "BM", method = "LL", REML = TRUE, error = TRUE)
        }, error = function(e) NULL)
        
        if (!is.null(fit_gls)) {
          Pinv <- fit_gls$sigma$Pinv
          eig_vals <- tryCatch(eigen(Pinv, symmetric = TRUE)$values, error = function(e) rep(NA, ncol(X)))
          eig_vals <- eig_vals[eig_vals > 0]
          
          cor_gls <- abs(cov2cor(Pinv))
          graph_gls <- igraph::graph_from_adjacency_matrix((cor_gls + t(cor_gls)) / 2,
                                                           mode = "undirected", weighted = TRUE, diag = FALSE)
          modularity_gls <- igraph::modularity(igraph::cluster_louvain(graph_gls))
          #clust_gls <- igraph::cluster_leiden(graph_gls, objective_function = "modularity")
          #modularity_gls <- igraph::modularity(clust_gls)
          
          trace <- sum(eig_vals)
          trace_per_species <- trace / nrow(X)
          lambda1 <- max(eig_vals)
          lambda1_prop <- lambda1 / trace
          condition_number <- log(lambda1 / min(eig_vals))
          
          prop_var <- eig_vals / sum(eig_vals)
          entropy <- -sum(prop_var * log(prop_var + 1e-12))
          effective_dim <- exp(entropy)
          evenness <- entropy / log(length(eig_vals))
          
          p()
          return(list(
            h3_cell = cell,
            success = TRUE,
            Pinv = Pinv,
            det_Pinv = det(Pinv),
            trace = trace,
            trace_per_species = trace_per_species,
            lambda1 = lambda1,
            lambda1_prop = lambda1_prop,
            condition_number = condition_number,
            effective_dimensionality = effective_dim,
            evenness = evenness,
            effective_dim_per_species = effective_dim / nrow(X),
            modularity = modularity_gls,
            n_species = nrow(X),
            n_traits = ncol(X),
            empirical = emp_metrics
          ))
        } else {
          return(list(h3_cell = cell, success = FALSE))
        }
      }, error = function(e) list(h3_cell = cell, success = FALSE))
    }, future.seed = TRUE)
  })
  
  successes <- Filter(function(x) x$success, result)
  failures  <- sapply(Filter(function(x) !x$success, result), `[[`, "h3_cell")
  
  metrics_df <- do.call(rbind, lapply(successes, function(x) {
    data.frame(
      h3_cell = x$h3_cell,
      det_Pinv = x$det_Pinv,
      trace = x$trace,
      trace_per_species = x$trace_per_species,
      lambda1 = x$lambda1,
      lambda1_prop = x$lambda1_prop,
      condition_number = x$condition_number,
      effective_dimensionality = x$effective_dimensionality,
      evenness = x$evenness,
      effective_dim_per_species = x$effective_dim_per_species,
      modularity = x$modularity,
      n_species = x$n_species,
      n_traits = x$n_traits,
      stringsAsFactors = FALSE
    )
  }))
  metrics_df$Pinv <- lapply(successes, `[[`, "Pinv")
  
  empirical_metrics_df <- do.call(rbind, lapply(successes, function(x) {
    data.frame(
      h3_cell = x$h3_cell,
      det_emp = x$empirical$det,
      trace_emp = x$empirical$trace,
      trace_per_species_emp = x$empirical$trace_per_species,
      lambda1_emp = x$empirical$lambda1,
      lambda1_prop_emp = x$empirical$lambda1_prop,
      condition_number_emp = x$empirical$condition_number,
      effective_dimensionality_emp = x$empirical$effective_dim,
      evenness_emp = x$empirical$evenness,
      effective_dim_per_species_emp = x$empirical$effective_dim_per_species,
      modularity_emp = x$empirical$modularity,
      n_species = x$n_species,
      n_traits = x$n_traits,
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(
    metrics_df = metrics_df,
    empirical_metrics_df = empirical_metrics_df,
    failed_cells = failures
  ))
}


#plotting processing functions
{
  
  
  #' Remove Polygons Crossing the International Date Line
  #'
  #' Identifies and removes geometries whose bounding box longitude range exceeds 
  #' a given threshold (default 180 degrees), which typically indicates crossing 
  #' of the International Date Line.
  remove_dateline_crossing <- function(sf_obj) {
    # Compute the longitude range for each polygon's bounding box
    bbox_ranges <- sapply(st_geometry(sf_obj), function(geom) {
      bbox <- st_bbox(geom)
      abs(bbox["xmax"] - bbox["xmin"])  # Calculate longitude range
    })
    
    # Set a threshold: remove polygons with a longitude range > 180 degrees
    threshold <- 180
    removed_count <- sum(bbox_ranges > threshold)  # Count polygons to be removed
    sf_filtered <- sf_obj[bbox_ranges <= threshold, ]
    
    # Print the number of polygons removed
    cat(sprintf("Number of polygons removed: %d\n", removed_count))
    
    return(sf_filtered)
  }
  
  #' Convert H3 Indices to Spatial Polygons in Parallel
  #'
  #' Converts a vector of H3 cell indices into an \code{sf} object containing 
  #' polygon geometries, using parallel processing for efficiency.
  convert_h3_indices_to_polygons <- function(h3_indices, num_cores = NULL) {
    # Default to using all available cores minus 1 if num_cores is not provided
    if (is.null(num_cores)) {
      num_cores <- max(parallel::detectCores() - 1, 1)
    }
    
    # Internal function to convert a single H3 index to an sf polygon
    convert_h3_to_polygon <- function(h3_index) {
      if (is.na(h3_index)) {
        return(NULL)  # Handle missing indices gracefully
      }
      # Use h3jsr::cell_to_polygon to convert H3 index to geometry
      polygon <- h3jsr::cell_to_polygon(h3_index, simple = FALSE)
      return(polygon)
    }
    
    # Apply the conversion in parallel using pbmclapply
    h3_polygons_list <- pbmclapply(
      h3_indices,
      convert_h3_to_polygon,  # Function to apply
      mc.cores = num_cores    # Number of cores to use
    )
    
    # Filter out NULL results (if conversion failed for some H3 indices)
    h3_polygons_list <- h3_polygons_list[!sapply(h3_polygons_list, is.null)]
    
    # Combine the list of polygons into a single sf object
    if (length(h3_polygons_list) == 0) {
      stop("No valid H3 polygons were generated. Check your H3 indices.")
    }
    
    # Combine all polygons into a single sf object
    h3_polygons_sf <- do.call(rbind, h3_polygons_list)
    return(h3_polygons_sf)
  }
  
  #' Precompute H3 Cell Polygons and Attach to Metrics Data
  #'
  #' Generates spatial polygons for H3 cell indices from a metrics data frame,
  #' processes them in parallel, and applies geographic cleaning.
  precompute_h3_polygons <- function(results, num_cores = NULL) {
    # Extract metrics_df
    metrics_df <- results$metrics_df
    
    # Step 1: Check if the h3_cell column exists
    if (!"h3_cell" %in% names(metrics_df)) {
      stop("The `metrics_df` does not contain the required `h3_cell` column.")
    }
    
    # Step 2: Generate H3 polygons in parallel
    cat("Precomputing H3 polygons in parallel...\n")
    polygons_sf <- convert_h3_indices_to_polygons(metrics_df$h3_cell, num_cores)
    
    # Add H3 index as a column to polygons_sf
    polygons_sf$h3_index <- metrics_df$h3_cell
    
    # Step 3: Remove polygons crossing the International Date Line
    cat("Removing polygons crossing the International Date Line...\n")
    #polygons_sf <- remove_dateline_crossing(polygons_sf)
    polygons_sf <- st_wrap_dateline(polygons_sf)
    
    # Return the filtered polygons_sf
    return(polygons_sf)
  }
  

  #' Plot Global H3 Metrics
  #'
  #' Creates a global visualization of a selected metric for H3 grid cells, overlayed on a world map.
  plot_h3_metrics_global <- function(results, polygons_sf, metric_column, crs = "EPSG:4326", 
                                     ocean_color = "lightblue", border_thickness = 0.5,
                                     log_transform = FALSE) {
    library(sf)
    library(ggplot2)
    library(rnaturalearth)
    library(viridis)
    
    # Extract metrics_df
    metrics_df <- results$metrics_df
    
    # Step 1: Check if the metric_column exists in metrics_df
    if (!metric_column %in% names(metrics_df)) {
      stop(paste("The specified metric_column:", metric_column, "does not exist in the provided metrics_df."))
    }
    
    # Optional: Apply natural log transform
    if (log_transform) {
      metrics_df[[metric_column]] <- log(metrics_df[[metric_column]])
    }
    
    # Step 2: Dynamically generate the label and title based on metric_column
    metric_label <- gsub("_", " ", metric_column)  # Replace underscores with spaces
    metric_label <- tools::toTitleCase(metric_label)  # Capitalize words
    title <- paste("Global Visualization of", "\n", metric_label)
    
    # Step 3: Merge precomputed polygons with metrics_df
    merged_sf <- merge(polygons_sf, metrics_df[, c("h3_cell", metric_column)], 
                       by.x = "h3_index", by.y = "h3_cell")
    
    # Step 4: Load world map
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    
    # Step 5: Create an ocean polygon
    ocean <- st_polygon(list(cbind(
      c(seq(-180, 179, length.out = 100), rep(180, 100), 
        seq(179, -180, length.out = 100), rep(-180, 100)),
      c(rep(-90, 100), seq(-89, 89, length.out = 100),
        rep(90, 100), seq(89, -90, length.out = 100))
    ))) |>
      st_sfc(crs = "WGS84") |>
      st_as_sf()
    
    # Step 6: Transform ocean polygon to CRS
    ocean_transformed <- st_transform_proj(ocean, crs)
    world_transformed <- st_transform_proj(world, crs)
    merged_sf_transformed <- st_transform_proj(merged_sf, crs)
    
    # Step 7: Extract ocean border
    ocean_border <- st_cast(ocean_transformed, "MULTILINESTRING")
    
    # Step 8: Plot with the specified CRS
    ggplot() +
      geom_sf(data = ocean_transformed, fill = ocean_color, color = NA) +
      geom_sf(data = world_transformed, fill = "black", color = "black", size = 0.2) +
      geom_sf(data = merged_sf_transformed, aes_string(fill = metric_column), color = NA) +
      scale_fill_viridis_c(option = "C", name = metric_label, na.value = "white") +
      geom_sf(data = ocean_border, color = "black", linewidth = border_thickness) +
      coord_sf(crs = crs, expand = FALSE) +
      theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "gray75", linetype = "dashed", linewidth = 0.25),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = paste(title, "\n", "in CRS:", crs),
        x = "Longitude",
        y = "Latitude"
      ) 
  }
  
  
  #' Plot Global H3 Metrics with Custom Color Palette
  #'
  #' Creates a global visualization of a selected metric for H3 grid cells, using a
  #' fully custom color palette and optional graticule customization.
  plot_h3_metrics_global.custom <- function(results, polygons_sf, metric_column, custom_palette = NULL, crs = "EPSG:4326", 
                                            ocean_color = "lightblue", border_thickness = 0.5, 
                                            lat_graticule_spacing = 15, lon_graticule_spacing = 15, 
                                            graticule_linetype = 1, graticule_color = "gray50", 
                                            graticule_thickness = 0.5) {
    library(sf)
    library(ggplot2)
    library(rnaturalearth)
    
    # Extract metrics_df
    metrics_df <- results$metrics_df
    
    # Step 1: Check if the metric_column exists in metrics_df
    if (!metric_column %in% names(metrics_df)) {
      stop(paste("The specified metric_column:", metric_column, "does not exist in the provided metrics_df."))
    }
    
    # Step 2: Dynamically generate the label and title based on metric_column
    metric_label <- gsub("_", " ", metric_column)  # Replace underscores with spaces
    metric_label <- tools::toTitleCase(metric_label)  # Capitalize words
    title <- paste("Global Visualization of", "\n", metric_label)
    
    # Step 3: Merge precomputed polygons with metrics_df
    merged_sf <- merge(polygons_sf, metrics_df[, c("h3_cell", metric_column)], 
                       by.x = "h3_index", by.y = "h3_cell")
    
    # Step 4: Check for custom palette
    if (!is.null(custom_palette)) {
      # Ensure that the custom palette covers all H3 indices in the merged_sf
      merged_sf$custom_color <- custom_palette[merged_sf$h3_index]
      # Assign "white" to any H3 cells not present in the custom palette
      merged_sf$custom_color[is.na(merged_sf$custom_color)] <- "white"
    } else {
      stop("Custom palette is required but was not provided.")
    }
    
    # Step 5: Load world map
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    
    # Step 6: Create an ocean polygon
    ocean <- st_polygon(list(cbind(
      c(seq(-180, 179, length.out = 100), rep(180, 100), 
        seq(179, -180, length.out = 100), rep(-180, 100)),
      c(rep(-90, 100), seq(-89, 89, length.out = 100),
        rep(90, 100), seq(89, -90, length.out = 100))
    ))) |>
      st_sfc(crs = "WGS84") |>
      st_as_sf()
    
    # Step 7: Transform ocean polygon to CRS
    ocean_transformed <- st_transform_proj(ocean, crs)
    world_transformed <- st_transform_proj(world, crs)
    merged_sf_transformed <- st_transform_proj(merged_sf, crs)
    
    # Step 8: Extract ocean border
    ocean_border <- st_cast(ocean_transformed, "MULTILINESTRING")
    
    # Step 9: Validate user inputs for graticules
    if (!is.numeric(lat_graticule_spacing) || lat_graticule_spacing <= 0 || lat_graticule_spacing %% 2 != 0) {
      stop("Latitude graticule spacing must be a positive even integer.")
    }
    if (!is.numeric(lon_graticule_spacing) || lon_graticule_spacing <= 0 || lon_graticule_spacing %% 2 != 0) {
      stop("Longitude graticule spacing must be a positive even integer.")
    }
    if (!graticule_linetype %in% 1:6) {
      stop("Graticule line type (lty) must be an integer between 1 and 6.")
    }
    
    # Step 10: Compute graticules from 0 outward
    lat_ticks <- seq(0, 90, by = lat_graticule_spacing)
    lat_ticks <- c(-rev(lat_ticks), lat_ticks)  # Add north and south graticules
    lon_ticks <- seq(0, 180, by = lon_graticule_spacing)
    lon_ticks <- c(-rev(lon_ticks), lon_ticks)  # Add east and west graticules
    
    graticules <- st_graticule(lat = lat_ticks, lon = lon_ticks, crs = crs)
    
    # Step 11: Determine graticule labeling
    label_axes <- "--EN"  # Default: suppress latitude labels on top, but include longitude and latitude labels
    if (grepl("moll", crs)) {
      label_axes <- "--N"  # Suppress longitude labels completely for Mollweide or similar projections
    }
    
    # Step 12: Plot with the specified CRS
    ggplot() +
      geom_sf(data = ocean_transformed, fill = ocean_color, color = NA) +  # Ocean polygon with user-specified color
      geom_sf(data = world_transformed, fill = "black", color = "black", size = 0.2) +  # World map (black continents)
      geom_sf(data = merged_sf_transformed, aes(fill = custom_color), color = NA) +  # H3 polygons with custom colors
      scale_fill_identity(name = metric_label, na.value = "white") +  # Use identity scale for direct color values
      geom_sf(data = graticules, color = graticule_color, linetype = graticule_linetype, linewidth = graticule_thickness) +  # Graticules on top
      geom_sf(data = ocean_border, color = "black", linewidth = border_thickness) +  # Ocean border plotted last
      coord_sf(crs = crs, expand = FALSE, label_axes = label_axes) +  # Adjust graticule labels
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),  # Remove default graticules
        panel.border = element_blank(),  # Remove default plot border
        panel.background = element_blank(),  # Remove default panel background
        plot.background = element_rect(fill = "white", color = NA)  # Set plot background to white
      ) +
      labs(
        title = paste(title, "\n", "in CRS:", crs),
        x = "Longitude",
        y = "Latitude"
      ) 
  }
  
  #' Plot Global H3 Metrics with Custom Palette and Raster Basemap
  #'
  #' Creates a global map of H3-based metrics with user-specified colors, optional 
  #' raster basemap, and a mask applied to non-ocean areas.
  plot_h3_metrics_global.custom.raster <- function(results, polygons_sf, metric_column, custom_palette = NULL, 
                                                   crs = "EPSG:4326", ocean_color = "lightblue", 
                                                   border_thickness = 0.5, lat_graticule_spacing = 15, 
                                                   lon_graticule_spacing = 15, graticule_linetype = 1, 
                                                   graticule_color = "gray50", graticule_thickness = 0.5,
                                                   raster_basemap = NULL, raster_alpha = 1, mask_color = "white") {
    library(sf)
    library(ggplot2)
    library(rnaturalearth)
    library(tidyterra)
    
    # Extract metrics_df
    metrics_df <- results$metrics_df
    
    # Step 1: Check if the metric_column exists in metrics_df
    if (!metric_column %in% names(metrics_df)) {
      stop(paste("The specified metric_column:", metric_column, "does not exist in the provided metrics_df."))
    }
    
    # Step 2: Dynamically generate the label and title based on metric_column
    metric_label <- gsub("_", " ", metric_column)  # Replace underscores with spaces
    metric_label <- tools::toTitleCase(metric_label)  # Capitalize words
    title <- paste("Global Visualization of", "\n", metric_label)
    
    # Step 3: Merge precomputed polygons with metrics_df
    merged_sf <- merge(polygons_sf, metrics_df[, c("h3_cell", metric_column)], 
                       by.x = "h3_index", by.y = "h3_cell")
    
    # Step 4: Check for custom palette
    if (!is.null(custom_palette)) {
      # Ensure that the custom palette covers all H3 indices in the merged_sf
      merged_sf$custom_color <- custom_palette[merged_sf$h3_index]
      # Assign "white" to any H3 cells not present in the custom palette
      merged_sf$custom_color[is.na(merged_sf$custom_color)] <- "white"
    } else {
      stop("Custom palette is required but was not provided.")
    }
    
    # Step 5: Load world map
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    
    # Step 6: Create an ocean polygon
    ocean <- st_polygon(list(cbind(
      c(seq(-180, 179, length.out = 100), rep(180, 100), 
        seq(179, -180, length.out = 100), rep(-180, 100)),
      c(rep(-90, 100), seq(-89, 89, length.out = 100),
        rep(90, 100), seq(89, -90, length.out = 100))
    ))) |>
      st_sfc(crs = "WGS84") |>
      st_as_sf()
    
    # Step 7: Transform ocean polygon to CRS
    ocean_transformed <- st_transform_proj(ocean, crs)
    world_transformed <- st_transform_proj(world, crs)
    merged_sf_transformed <- st_transform_proj(merged_sf, crs)
    
    # Step 8: Extract ocean border
    ocean_border <- st_cast(ocean_transformed, "MULTILINESTRING")
    
    # Step 9: Validate user inputs for graticules
    if (!is.numeric(lat_graticule_spacing) || lat_graticule_spacing <= 0 || lat_graticule_spacing %% 1 != 0) {
      stop("Latitude graticule spacing must be a positive even integer.")
    } #also allow odd like 15
    if (!is.numeric(lon_graticule_spacing) || lon_graticule_spacing <= 0 || lon_graticule_spacing %% 1 != 0) {
      stop("Longitude graticule spacing must be a positive even integer.")
    } #also allow odd like 15
    if (!graticule_linetype %in% 1:6) {
      stop("Graticule line type (lty) must be an integer between 1 and 6.")
    }
    
    # Step 10: Compute graticules from 0 outward
    lat_ticks <- seq(0, 90, by = lat_graticule_spacing)
    lat_ticks <- c(-rev(lat_ticks), lat_ticks)  # Add north and south graticules
    lon_ticks <- seq(0, 180, by = lon_graticule_spacing)
    lon_ticks <- c(-rev(lon_ticks), lon_ticks)  # Add east and west graticules
    
    graticules <- st_graticule(lat = lat_ticks, lon = lon_ticks, crs = crs)
    
    # Step 11: Determine graticule labeling
    label_axes <- "--EN"  # Default: suppress latitude labels on top, but include longitude and latitude labels
    if (grepl("moll", crs)) {
      label_axes <- "--N"  # Suppress longitude labels completely for Mollweide or similar projections
    }
    
    # Step 12: Create a mask from the ocean polygon
    mask <- st_difference(st_as_sfc(st_bbox(ocean_transformed)), ocean_transformed)
    
    # Step 13: Plot with the specified CRS
    ggplot() +
      geom_sf(data = ocean_transformed, fill = ocean_color, color = NA) +  # Ocean polygon with user-specified color
      geom_sf(data = world_transformed, fill = "black", color = "black", size = 0.2) +  # World map (black continents)
      {
        if (!is.null(raster_basemap)) {
          # Add raster basemap above the ocean polygon
          geom_spatraster_rgb(data = raster_basemap, alpha = raster_alpha)
        }
      } +
      geom_sf(data = mask, fill = mask_color, color = NA) +  # Add mask layer
      geom_sf(data = merged_sf_transformed, aes(fill = custom_color), color = NA) +  # H3 polygons with custom colors
      scale_fill_identity(name = metric_label, na.value = "white") +  # Use identity scale for direct color values
      geom_sf(data = graticules, color = graticule_color, linetype = graticule_linetype, linewidth = graticule_thickness) +  # Graticules on top
      geom_sf(data = ocean_border, color = "black", linewidth = border_thickness) +  # Ocean border plotted last
      coord_sf(crs = crs, expand = FALSE, label_axes = label_axes) +  # Adjust graticule labels
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),  # Remove default graticules
        panel.border = element_blank(),  # Remove default plot border
        panel.background = element_blank(),  # Remove default panel background
        plot.background = element_rect(fill = "white", color = NA)  # Set plot background to white
      ) +
      labs(
        title = paste(title, "\n", "in CRS:", crs),
        x = "Longitude",
        y = "Latitude"
      ) 
  }
  
  
  #' Prepare VCV Data for Plotting
  #'
  #' Adds geographic latitude information to a VCV metrics data frame by converting 
  #' H3 cell indices to centroid points.
  prepare_vcv_plot_data <- function(vcv_df) {
    # Ensure h3_cell is character
    vcv_df$h3_cell <- as.character(vcv_df$h3_cell)
    
    # Compute centroids and extract latitude
    centroids <- h3jsr::cell_to_point(vcv_df$h3_cell)
    vcv_df$latitude <- sf::st_coordinates(centroids)[, 2]
    
    vcv_df
  }
  
  #' Plot VCV Metrics by Latitude
  #'
  #' Creates a scatterplot of a selected variance–covariance (VCV) metric against latitude,
  #' with optional regional coloring, smoothing, and reference lines.
  make_vcv_metric_plot <- function(df,
                                   metric,
                                   polygons_sf,
                                   point_size = 1.6,
                                   point_alpha = 0.6,
                                   add_loess = FALSE,
                                   zero_line_metrics = NULL,
                                   one_line_metrics = NULL,
                                   max_lat,
                                   color_by_region = FALSE,
                                   legend_alpha = 1) {
    # Metric display names
    metric_labels <- c(
      det_Pinv = "Determinant (Trait Volume)",
      log_det_Pinv = "Log Determinant (Trait Volume)",
      trace = "Trace (Total Variance)",
      trace_per_species = "Trace / Species (Dispersion)",
      lambda1 = "Largest Eigenvalue",
      lambda1_prop = "Proportion in PC1 (Axis Alignment)",
      condition_number = "Condition Number (Anisotropy)",
      log_condition_number = "Log Condition Number (Anisotropy)",
      effective_dimensionality = "Effective Dimensionality",
      evenness = "Variance Evenness",
      effective_dim_per_species = "Effective Dimensionality / log(Species)",
      modularity = "Modularity (Trait Covariance)"
    )
    
    # Remove NAs
    df_metric <- df %>% dplyr::filter(!is.na(.data[[metric]]))
    
    # Join region info if requested
    if (color_by_region) {
      df_metric <- df_metric %>%
        dplyr::left_join(
          polygons_sf %>%
            sf::st_drop_geometry() %>%
            dplyr::select(h3_address, region),
          by = c("h3_cell" = "h3_address")
        )
    }
    
    # Base plot: color by region or black
    if (color_by_region) {
      p <- ggplot2::ggplot(df_metric, aes(x = latitude, y = .data[[metric]], color = region)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +  # ADD THIS FIRST
        ggplot2::geom_point(size = point_size, alpha = point_alpha) +
        ggplot2::scale_color_manual(values = c(
          "Africa"            = "#E76F00",  # savanna orange
          "Antarctica"        = "#6BC4FF",  # iceberg blue
          "Asia"              = "#1FA87C",  # jade green
          "Europe"            = "#4858D9",  # slate indigo
          "North America"     = "#0085CA",  # bold azure
          "Central America"   = "#FFCC29",  # tropical yellow
          "South America"     = "#C72E29",  # rich crimson
          "Australia/Oceania" = "#D845A6",  # coral pink
          "Other"             = "#9E9E9E"   # neutral grey
        )) +
        ggplot2::guides(color = ggplot2::guide_legend(
          title = "Region",
          override.aes = list(size = 4, alpha = legend_alpha)
        ))
    } else {
      p <- ggplot2::ggplot(df_metric, aes(x = latitude, y = .data[[metric]])) +
        ggplot2::geom_point(color = "black", size = point_size, alpha = point_alpha)
    }
    
    # Optional loess smoothing
    if (add_loess) {
      p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE)
    }
    
    # Optional reference lines
    if (!is.null(zero_line_metrics) && metric %in% zero_line_metrics) {
      p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red")
    }
    if (!is.null(one_line_metrics) && metric %in% one_line_metrics) {
      p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dotted", color = "red")
    }
    
    # Final scales, labels, and theme
    p +
      ggplot2::scale_x_continuous(limits = c(-max_lat, max_lat), oob = scales::oob_squish) +
      ggplot2::scale_y_continuous(oob = scales::oob_squish) +
      ggplot2::labs(
        title = metric_labels[[metric]] %||% metric,
        x = "Latitude",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line  = element_line(color = "black", linewidth = 0.7),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        axis.text  = element_text(color = "black"),
        plot.title = element_text(hjust = 0, face = "bold", size = 13),
        axis.text.x= element_text(size = 11),
        plot.margin= margin(5, 10, 5, 10)
      )
  }
  
  #' Plot VCV Metrics by Latitude
  #'
  #' Creates a scatterplot of a selected variance–covariance (VCV) metric against latitude,
  #' with optional regional coloring, smoothing, and reference lines.
  plot_vcv_metrics_by_latitude_v2 <- function(
    results_vcv,
    metrics = c("trace_per_species", "lambda1_prop", "effective_dimensionality"),
    source = c("estimated", "empirical"),
    add_loess = FALSE,
    point_size = 1.6,
    point_alpha = 0.6,
    zero_line_metrics = NULL,
    one_line_metrics = NULL,
    color_by_region = FALSE,
    legend_alpha = 1) {
    library(dplyr)
    library(ggplot2)
    library(h3jsr)
    library(sf)
    library(patchwork)
    library(scales)
    
    source <- match.arg(source)
    
    # ---- Metric display names ----
    metric_labels <- c(
      det_Pinv = "Determinant (Trait Volume)",
      log_det_Pinv = "Log Determinant (Trait Volume)",
      trace = "Trace (Total Variance)",
      trace_per_species = "Trace / Species (Dispersion)",
      lambda1 = "Largest Eigenvalue",
      lambda1_prop = "Proportion in PC1 (Axis Alignment)",
      condition_number = "Condition Number (Anisotropy)",
      log_condition_number = "Log Condition Number (Anisotropy)",
      effective_dimensionality = "Effective Dimensionality",
      evenness = "Variance Evenness",
      effective_dim_per_species = "Effective Dimensionality / log(Species)",
      modularity = "Modularity (Trait Covariance)"
    )
    
    # pick the right dataframe and suffix
    suffix <- if (source == "empirical") "_emp" else ""
    df <- if (source == "empirical") results_vcv$empirical_metrics_df else results_vcv$metrics_df
    
    # add log columns if present
    if (paste0("det", suffix) %in% names(df)) {
      df[[paste0("log_det_Pinv", suffix)]] <- log(pmax(df[[paste0("det", suffix)]], 1e-12))
    }
    if (paste0("condition_number", suffix) %in% names(df)) {
      df[[paste0("log_condition_number", suffix)]] <- log(pmax(df[[paste0("condition_number", suffix)]], 1e-12))
    }
    
    # add effective_dim_per_species if possible
    if (paste0("effective_dimensionality", suffix) %in% names(df) &&
        "n_species" %in% names(df)) {
      df[[paste0("effective_dim_per_species", suffix)]] <-
        df[[paste0("effective_dimensionality", suffix)]] / log(df$n_species)
    }
    
    # prepare latitudes
    df <- prepare_vcv_plot_data(df)
    max_lat <- ceiling(max(abs(df$latitude), na.rm = TRUE))
    
    # build continent polygons once
    polygons_for_plot <- prepare_continent_overlap(polygons_sf.1)
    
    # suffix-aware metric names
    metrics_full <- paste0(metrics, suffix)
    valid_metrics <- intersect(metrics_full, names(df))
    if (length(valid_metrics) == 0) {
      stop("No valid metrics to plot. Available: ", paste(names(df), collapse = ", "))
    }
    
    # make each plot
    plots <- lapply(seq_along(valid_metrics), function(i) {
      raw <- metrics[i]
      col <- valid_metrics[i]
      # override the title in the metric labels
      old_label <- metric_labels[[raw]] %||% raw
      # call the metric plot
      make_vcv_metric_plot(
        df,
        col,
        polygons_for_plot,
        point_size    = point_size,
        point_alpha   = point_alpha,
        add_loess     = add_loess,
        zero_line_metrics = zero_line_metrics,
        one_line_metrics  = one_line_metrics,
        max_lat       = max_lat,
        color_by_region  = color_by_region
      ) +
        ggplot2::labs(title = old_label)
    })
    
    # stack and return
    patchwork::wrap_plots(plots, ncol = 3)
  }
  
  
  #' Load and Reproject a Raster
  #'
  #' Loads a raster file and reprojects it to a specified coordinate reference system (CRS),
  #' optionally adjusting the resolution.
  load_reproject_raster <- function(raster_path, target_crs, res = NULL, method = "bilinear") {
    library(terra)
    
    # Load the raster
    raster <- rast(raster_path)
    
    # Check if the CRS of the raster matches the target CRS
    if (!crs(raster) == target_crs) {
      # Reproject the raster, setting resolution if specified
      raster <- project(raster, crs = target_crs, res = res, method = method)
    } else {
      # If no reprojection is needed but resolution is specified, adjust resolution
      if (!is.null(res)) {
        raster <- resample(raster, resolution = res, method = method)
      }
    }
    
    return(raster)
  }
  
  
  
}
  


#' Process Bioclimatic Raster Data for H3 Polygons
#'
#' Extracts and summarizes bioclimatic raster values for a set of H3 polygons.
process_bioclim <- function(raster_layer, polygons_sf, summary_function = mean, na_rm = TRUE) {
  # Extract summary values for each polygon
  extracted_values <- extract(raster_layer, vect(polygons_sf), fun = summary_function, na.rm = na_rm, ID = FALSE)
  
  # Rename the column to include the raster layer name and summary metric
  raster_name <- gsub(".tif$", "", basename(sources(raster_layer)))  # Extract the base name
  colnames(extracted_values)[1] <- paste0(raster_name, "_", deparse(substitute(summary_function)))
  
  # Add H3 addresses to the extracted data
  extracted_values$h3 <- polygons_sf$h3_address
  
  # Return the processed data frame
  return(extracted_values)
}

#' Load and Process All Bioclimatic Raster Files in a Directory
#'
#' Iterates over all `.tif` bioclimatic raster files in a directory, extracts summary
#' statistics for each H3 polygon, and merges the results into a spatial dataset.
load_all_bioclim <- function(bioclim_dir, spatial_coords, polygons_sf, summary_function = mean, summary_function_name = NULL, na_rm = TRUE) {
  # List all bioclimatic raster files in the specified directory
  bioclim_files <- list.files(bioclim_dir, pattern = "\\.tif$", full.names = TRUE)
  
  # Print a summary of the files found
  cat("Bioclimatic files found:\n")
  print(bioclim_files)
  
  # Determine the summary function name
  function_name <- if (!is.null(summary_function_name)) {
    summary_function_name
  } else {
    deparse(substitute(summary_function))
  }
  
  # Iterate through each bioclimatic file
  for (file in bioclim_files) {
    cat(paste("Processing file:", file, '\n'))
    
    # Load the bioclimatic layer
    bioclim_layer <- tryCatch({
      rast(file)
    }, error = function(e) {
      cat(paste("Error loading", file, ":", e$message, '\n'))
      return(NULL)
    })
    
    # Skip to the next file if loading failed
    if (is.null(bioclim_layer)) next
    
    # Process the bioclim layer using the provided function
    extracted_data <- tryCatch({
      process_bioclim(bioclim_layer, polygons_sf, summary_function, na_rm)
    }, error = function(e) {
      cat(paste("Error processing", file, ":", e$message, '\n'))
      return(NULL)
    })
    
    # Skip to the next file if processing failed
    if (is.null(extracted_data)) next
    
    # Rename the column in case it wasn’t set properly in process_bioclim
    raster_name <- gsub(".tif$", "", basename(file))
    colnames(extracted_data)[1] <- paste0(raster_name, "_", function_name)
    
    # Merge the extracted data with spatial_coords
    cat(paste("Merging extracted data from", file, "into spatial_coords\n"))
    spatial_coords <- tryCatch({
      merge(spatial_coords, extracted_data, by = "h3", all.x = TRUE)
    }, error = function(e) {
      cat(paste("Error merging data from", file, ":", e$message, '\n'))
      return(spatial_coords)
    })
  }
  
  # Return the updated spatial_coords dataset
  return(spatial_coords)
}



#' Leave-One-Out Model Weights for Spatial Lag Models
#'
#' Performs leave-one-term-out model comparisons for spatial lag models 
#' (fitted with \code{\link[spatialreg]{lagsarlm}}) by refitting the model 
#' after removing each candidate predictor, and computing AIC/BIC-based weights.
loo_model_weights <- function(model, data, listw, 
                              control_terms = c("(Intercept)", 
                                                "asin(sqrt(sampling_frac))", 
                                                "tipDR_range_weighted_mean_z"),
                              response_var = "mean_lineage_rate_range_weighted",
                              n_terms = NULL,
                              recompute_full_model = FALSE) {
  # Build the design matrix from the original model and data
  mm <- model.matrix(formula(model), data = data)
  
  # Get expanded predictor names
  all_terms_expanded <- colnames(mm)
  
  # Define candidate columns: all columns except control terms
  candidate_cols <- setdiff(all_terms_expanded, control_terms)
  
  # Optionally limit candidate columns to first n_terms if specified
  if (!is.null(n_terms)) {
    candidate_cols <- candidate_cols[1:min(n_terms, length(candidate_cols))]
  }
  
  # Extract the response variable from the data
  resp <- data[[response_var]]
  
  # Optionally refit the full model to ensure consistency
  if (recompute_full_model) {
    # Construct full model data the same way as reduced models
    full_data <- as.data.frame(mm)  # Use design matrix
    full_data <- full_data[, setdiff(colnames(full_data), "(Intercept)"), drop = FALSE]  # Drop intercept
    full_data[[response_var]] <- resp  # Add response variable
    
    full_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                           data = full_data,
                           listw = listw,
                           zero.policy = TRUE,
                           method = "eigen")
    full_AIC <- AIC(full_model)
    full_BIC <- BIC(full_model)
  } else {
    full_AIC <- AIC(model)
    full_BIC <- BIC(model)
    full_model <- model  # Keep the original model if not recomputed
  }
  
  # Helper function: remove one column, refit the model, and return a list with AIC, BIC, and model object
  fit_after_removal <- function(term) {
    new_cols <- setdiff(colnames(mm), term)
    new_data <- as.data.frame(mm[, new_cols, drop = FALSE])
    # Remove the intercept column (it is added automatically by the formula)
    new_data <- new_data[, setdiff(colnames(new_data), "(Intercept)"), drop = FALSE]
    new_data[[response_var]] <- resp
    updated_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                              data = new_data,
                              listw = listw,
                              zero.policy = TRUE,
                              method = "eigen")
    list(AIC = AIC(updated_model), BIC = BIC(updated_model), model = updated_model)
  }
  
  # Use pblapply (parallel with progress bar) to iterate over candidate columns
  res_list <- pblapply(candidate_cols, fit_after_removal)
  
  # Extract AIC and BIC values for candidate models
  aic_vals <- sapply(res_list, function(x) x$AIC)
  bic_vals <- sapply(res_list, function(x) x$BIC)
  
  results <- data.frame(
    Term = candidate_cols,
    AIC = aic_vals,
    BIC = bic_vals,
    stringsAsFactors = FALSE
  )
  
  # Compute ΔAIC and ΔBIC relative to the full model's AIC and BIC
  results$DeltaAIC <- results$AIC - full_AIC
  results$DeltaBIC <- results$BIC - full_BIC
  
  # Compute weights using the aicw function (for both AIC and BIC)
  results$AIC_Weight <- sapply(results$AIC, function(candidate_AIC) 
    aicw(c(full_AIC, candidate_AIC))$aicweights[1]
  )
  results$BIC_Weight <- sapply(results$BIC, function(candidate_BIC) 
    aicw(c(full_BIC, candidate_BIC))$aicweights[1]
  )
  
  # Create a named list of the updated model objects
  model_list <- lapply(res_list, function(x) x$model)
  names(model_list) <- candidate_cols
  
  # Prepend the recomputed full model (if applicable)
  if (recompute_full_model) {
    model_list <- c(list("Full Model" = full_model), model_list)
  }
  
  # Add full model metrics to the results for reference
  full_metrics <- data.frame(
    Term = "Full Model",
    AIC = full_AIC,
    BIC = full_BIC,
    DeltaAIC = 0,
    DeltaBIC = 0,
    AIC_Weight = 1,  # Full model always gets weight 1 in its own comparison
    BIC_Weight = 1,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(full_metrics, results)
  results <- results[order(results$AIC), ]
  
  list(metrics = results, models = model_list)
}


#' Leave-One-Out Model Weights from Existing Spatial Lag Model Object
#'
#' Computes leave-one-term-out model comparisons for an existing spatial lag model 
#' object fitted with \code{\link[spatialreg]{lagsarlm}}, using its internal design 
#' matrix and response vector.
loo_model_weights_from_obj <- function(mod, listw, 
                                       control_terms = c("(Intercept)", 
                                                         "asin(sqrt(sampling_frac))", 
                                                         "tipDR_range_weighted_mean_z"),
                                       response_var = "y",  # assumes mod$y holds the response
                                       n_terms = NULL) {
  # Extract the design matrix and response from the model object
  X <- mod$X
  y <- mod$y
  
  # Get all column names from the design matrix
  all_terms <- colnames(X)
  
  # Define candidate columns by removing the control terms
  candidate_cols <- setdiff(all_terms, control_terms)
  
  # Optionally limit candidate columns to the first n_terms if specified
  if (!is.null(n_terms)) {
    candidate_cols <- candidate_cols[1:min(n_terms, length(candidate_cols))]
  }
  
  # Full model metrics (from the passed model object)
  full_AIC <- AIC(mod)
  full_BIC <- BIC(mod)
  full_model <- mod
  
  # Helper function: remove one candidate column, then refit the model
  fit_after_removal <- function(term) {
    new_cols <- setdiff(all_terms, term)
    new_data <- as.data.frame(X[, new_cols, drop = FALSE])
    # Remove the intercept column if present (it will be added automatically)
    new_data <- new_data[, setdiff(colnames(new_data), "(Intercept)"), drop = FALSE]
    new_data[[response_var]] <- y  # Add the response vector
    updated_model <- lagsarlm(as.formula(paste(response_var, "~ .")),
                              data = new_data,
                              listw = listw,
                              zero.policy = mod$zero.policy,
                              method = mod$method)
    list(AIC = AIC(updated_model), BIC = BIC(updated_model), model = updated_model)
  }
  
  # Use pblapply to iterate over candidate columns with a progress bar
  res_list <- pbapply::pblapply(candidate_cols, fit_after_removal)
  
  # Extract AIC and BIC values from each reduced model
  aic_vals <- sapply(res_list, function(x) x$AIC)
  bic_vals <- sapply(res_list, function(x) x$BIC)
  
  results <- data.frame(
    Term = candidate_cols,
    AIC = aic_vals,
    BIC = bic_vals,
    stringsAsFactors = FALSE
  )
  
  # Compute ΔAIC and ΔBIC relative to the full model's AIC/BIC
  results$DeltaAIC <- results$AIC - full_AIC
  results$DeltaBIC <- results$BIC - full_BIC
  
  # Compute weights using the aicw function (ensure aicw() is defined/loaded)
  results$AIC_Weight <- sapply(results$AIC, function(candidate_AIC) 
    aicw(c(full_AIC, candidate_AIC))$aicweights[1]
  )
  results$BIC_Weight <- sapply(results$BIC, function(candidate_BIC) 
    aicw(c(full_BIC, candidate_BIC))$aicweights[1]
  )
  
  # Create a named list of the updated model objects from the reduced fits
  model_list <- lapply(res_list, function(x) x$model)
  names(model_list) <- candidate_cols
  
  # Prepend the full model to the list of models
  model_list <- c(list("Full Model" = full_model), model_list)
  
  # Add full model metrics for reference
  full_metrics <- data.frame(
    Term = "Full Model",
    AIC = full_AIC,
    BIC = full_BIC,
    DeltaAIC = 0,
    DeltaBIC = 0,
    AIC_Weight = 1,
    BIC_Weight = 1,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(full_metrics, results)
  results <- results[order(results$AIC), ]
  
  list(metrics = results, models = model_list)
}





