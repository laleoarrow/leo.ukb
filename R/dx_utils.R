#' Internal robust dx runner
#' @keywords internal
#' @noRd
.dx_run <- function(args, intern = TRUE, ignore.stderr = FALSE, ...) {
  # Find dx
  dx_path <- Sys.which("dx")

  # macOS Homebrew fallback
  if (dx_path == "" && Sys.info()["sysname"] == "Darwin") {
    common_paths <- c("/opt/homebrew/bin/dx", "/usr/local/bin/dx")
    for (p in common_paths) {
      if (file.exists(p)) {
        dx_path <- p
        break
      }
    }
  }

  if (dx_path == "") dx_path <- "dx"
  
  # Log command for transparency
  # Quote args that have spaces for display clarity (system2 handles actual args)
  pretty_args <- sapply(args, function(x) if (grepl(" ", x)) shQuote(x) else x)
  leo.basic::leo_log("Executing: dx {paste(pretty_args, collapse = ' ')}", level = "info")

  # Use system2
  if (ignore.stderr) {
    res <- suppressWarnings(system2(dx_path, args = args, stdout = intern, stderr = FALSE, ...))
  } else {
    res <- system2(dx_path, args = args, stdout = intern, stderr = TRUE, ...)
  }

  # Handle exit status if intern=TRUE
  if (intern && !is.null(attr(res, "status")) && attr(res, "status") != 0 && !ignore.stderr) {
    # It might be better to return NULL or handle at call site, but we'll stick to original behavior
  }

  return(res)
}

#' Get Dataset Dictionary (Official Method)
#'
#' Downloads the full data dictionary for a specific dataset using `dx extract_dataset -ddd`.
#' This file contains the authoritative list of column names (e.g., `p21003_i0`).
#'
#' @param dataset The dataset ID or name (e.g. "app123.dataset" or "project-X:record-Y").
#' @return A data.table/data.frame containing the dictionary (must have a `name` column).
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Download and load dictionary for a specific dataset
#'   dict <- dx_get_dataset_dictionary("project-X:record-Y")
#'   
#'   # Uses cached version if available
#'   dict2 <- dx_get_dataset_dictionary("project-X:record-Y")
#' }
dx_get_dataset_dictionary <- function(dataset) {
  # Determine storage location: prefer "tmp" in project, otherwise tempdir()
  storage_dir <- if (dir.exists("tmp")) "tmp" else tempdir()
  
  # Optimization: Try to resolve the record name to check if the file already exists.
  # We use `dx describe` to get the Name.
  
  # Step 1: Resolve dataset ID if needed
  if (missing(dataset) || is.null(dataset)) stop("Dataset ID is required.")

  leo.basic::leo_log("Checking for info on {dataset}...", level = "info")
  
  # Get dataset name to predict filename
  # Command: dx describe <id> --json
  # We parse the "name" field.
  # Since we might not have jsonlite, we can try `dx describe` text output.
  # "Name" is usually the first line or labeled "Name".
  
  sys_cmd_args <- c("describe", dataset, "--json")
  desc_json <- tryCatch({
    .dx_run(sys_cmd_args, intern = TRUE, ignore.stderr = TRUE)
  }, error = function(e) NULL)
  
  record_name <- NULL
  if (!is.null(desc_json) && length(desc_json) > 0) {
      json_str <- paste(desc_json, collapse = " ")
      matches <- regmatches(json_str, regexec('"name"\\s*:\\s*"([^"]+)"', json_str))
      if (length(matches[[1]]) > 1) {
          record_name <- matches[[1]][2]
      }
  }
  
  target_file <- NULL
  if (!is.null(record_name)) {
      expected_file <- file.path(storage_dir, paste0(record_name, ".data_dictionary.csv"))
      if (file.exists(expected_file)) {
          leo.basic::leo_log("Using cached dictionary: {expected_file}", level = "success")
          return(data.table::fread(expected_file, data.table = FALSE))
      }
  }
  
  # If not cached or name resolution failed, proceed to download.
  leo.basic::leo_log("Fetching dictionary for {dataset}...", level = "info")
  
  # Use system2 for proper exit code handling
  # dx extract_dataset <id> -ddd --output <dir>/
  
  res <- tryCatch({
    system2("dx", 
            args = c("extract_dataset", dataset, "-ddd", "--output", paste0(storage_dir, "/")),
            stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    leo.basic::leo_log("Failed to run dx extract_dataset: {e$message}", level = "danger")
    return(NULL)
  })
  
  # Check exit status
  exit_code <- attr(res, "status")
  if (!is.null(exit_code) && exit_code != 0) {
    leo.basic::leo_log("dx extract_dataset failed with exit code {exit_code}", level = "danger")
    if (length(res) > 0) {
      leo.basic::leo_log("Output: {paste(res, collapse = '\\n')}", level = "info")
    }
    return(NULL)
  }
  
  if (is.null(res)) return(NULL)
  
  # Find the generated CSV file
  # We can't trust record_name 100% if regex failed, so we search directory again by time.
  # But we should search specifically for .data_dictionary.csv
  
  csv_files <- list.files(storage_dir, pattern = "\\.data_dictionary\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    leo.basic::leo_log("No dictionary CSV found after extraction.", level = "danger")
    return(NULL)
  }
  
  # Sort by mtime to pick the absolute newest one (the one we just downloaded)
  info <- file.info(csv_files)
  target_file <- rownames(info)[which.max(info$mtime)]
  
  leo.basic::leo_log("Dictionary downloaded to: {target_file}", level = "success")
  
  # Read it
  dict_df <- data.table::fread(target_file, data.table = FALSE)

  
  return(dict_df)
}


#' Get UKB Dataset Schema (Metadata)
#' 
#' Fetches official UK Biobank schema files (field.tsv, category.tsv).
#' Prioritizes:
#' 1. Local cache in `tmp/`
#' 2. DNAnexus project `Showcase metadata/` folder
#' 3. Official UKB Showcase website (public download)
#'
#' @param type type of schema: "field", "category" or "hierarchy"
#' @param force Logical. If TRUE, re-downloads the schema even if a local cache exists.
#' @return Path to the downloaded TSV file
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Fetch Field Schema (Schema 1)
#'   field_schema_path <- dx_get_schema("field")
#'   
#'   # Fetch Category Schema (Schema 2)
#'   cat_schema_path <- dx_get_schema("category")
#' }
dx_get_schema <- function(type = c("field", "category", "hierarchy"), force = FALSE) {
  type <- match.arg(type)
  
  # Schema IDs: 1=Field, 3=Category, 13=Hierarchy
  schema_id <- switch(type, "field" = 1, "category" = 3, "hierarchy" = 13)
  filename <- paste0(type, ".tsv")
  
  # Define paths
  storage_dir <- if (dir.exists("tmp")) "tmp" else tempdir()
  local_path <- file.path(storage_dir, filename)
  
  # 1. Check Local Cache
  if (!force && file.exists(local_path) && file.size(local_path) > 0) {
    # Optional: check age? For now assume stable.
    leo.basic::leo_log("Using cached schema: {local_path}", level = "info")
    return(local_path)
  }
  # 2. Try DNAnexus (official location in RAP)
  # Path: "/Showcase metadata/<type>.tsv"
  leo.basic::leo_log("Looking for {filename}...", level = "info")
  rap_path <- paste0("Showcase metadata/", filename)
  res <- tryCatch({
      .dx_run(c("download", rap_path, "-o", local_path, "-f"), intern = TRUE, ignore.stderr = TRUE)
      TRUE
  }, error = function(e) FALSE)
  
  if (file.exists(local_path) && file.size(local_path) > 0) {
      leo.basic::leo_log("Downloaded {filename} from RAP project.", level = "success")
      return(local_path)
  }
  # 3. Fallback: Download from UK Biobank Public Showcase
  # URL: https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id={1,2}
  leo.basic::leo_log("Metadata not found in project. Fetching from UKBiobank Showcase...", level = "info")
  url <- paste0("https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=", schema_id)
  tryCatch({
      utils::download.file(url, destfile = local_path, quiet = TRUE, method = "auto")
  }, error = function(e) {
      leo.basic::leo_log("Failed to download from UKB Showcase: {e$message}", level = "danger")
  })
  if (file.exists(local_path) && file.size(local_path) > 0) {
      leo.basic::leo_log("Downloaded {filename} from UKB Showcase.", level = "success")
      return(local_path)
  }
  return(NULL)
}

#' Find Valid Columns in Dataset (Official Regex)
#'
#' Uses the official UKB RAP regex strategy to find all columns belonging to a field ID.
#' Filters by entity to ensure only fields for the specified entity are returned.
#'
#' @param fields Vector of field IDs (e.g. `21003`, `p31`).
#' @param dictionary The dictionary dataframe from `dx_get_dataset_dictionary`.
#' @param entity The entity type to filter for (e.g. "participant", "hesin"). 
#'   Default is "participant". The dictionary contains multiple entities, so filtering is essential.
#' @return Character vector of valid column names found in the dictionary for the specified entity.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   dict <- dx_get_dataset_dictionary("project-X:record-Y")
#'   
#'   # Find columns for specific fields
#'   cols <- dx_find_columns(c("p31", "21003"), dict)
#' }
dx_find_columns <- function(fields, dictionary, entity = "participant") {
  if (is.null(dictionary) || nrow(dictionary) == 0) return(character(0))
  if (!("name" %in% names(dictionary))) {
    stop("Dictionary is missing 'name' column.")
  }
  
  # CRITICAL: Filter dictionary by entity first
  # The data_dictionary.csv contains ALL entities (participant, hesin, death, etc.)
  # We must only match fields for the requested entity
  if ("entity" %in% names(dictionary)) {
    dictionary <- dictionary[dictionary$entity == entity, ]
    if (nrow(dictionary) == 0) {
      leo.basic::leo_log("No fields found for entity '{entity}' in dictionary.", level = "warning")
      return(character(0))
    }
  } else {
    leo.basic::leo_log("Dictionary has no 'entity' column. Cannot filter by entity (risky).", level = "warning")
  }
  
  # Clean fields -> extract purely numeric ID
  # 21003 -> "21003", p31 -> "31"
  fields_char <- as.character(fields)
  clean_ids <- sub("^\\D*(\\d+).*", "\\1", fields_char)
  clean_ids <- unique(clean_ids[clean_ids != ""]) # Remove empties
  
  valid_cols <- character(0)
  all_names <- dictionary$name
  for (fid in clean_ids) {
    # Official tutorial: https://github.com/UK-Biobank/UKB-RAP-Notebooks-Access/blob/main/RStudio/A110_Export_participant_data.Rmd
    regex <- paste0("^p", fid, "(?![0-9])")
    matches <- grep(regex, all_names, value = TRUE, perl = TRUE) 
    if (length(matches) > 0) {
      valid_cols <- c(valid_cols, matches)
    } else {
      leo.basic::leo_log("Field {fid} (p{fid}) not found for entity '{entity}'.", level = "warning")
    }
  }
  
  return(unique(valid_cols))
}

#' Find Fields by Category ID
#'
#' Scans the official schema for fields belonging to a specific category ID.
#'
#' @param category_ids Vector of numeric Category IDs (e.g. `1712`).
#' @param entity Entity type to filter for (default: "participant"). Ignored if using schema.
#' @return Character vector of field IDs (e.g. `c("41270", "41271")`).
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Find fields for Category 1712 (First occurrences)
#'   # Automatically downloads schema if missing
#'   fields_1712 <- dx_find_fields_by_category(1712)
#'   
#'   # Use vector of categories
#'   fields_mix <- dx_find_fields_by_category(c(1712, 100))
#' }
#' @importFrom data.table fread 
dx_find_fields_by_category <- function(category_ids, entity = "participant") {
  if (is.null(category_ids) || length(category_ids) == 0) return(character(0))
  
  # 1. Expand categories via Hierarchy (resolve parents to children)
  target_cats <- suppressWarnings(as.integer(category_ids))
  target_cats <- target_cats[!is.na(target_cats)]
  
  h_path <- dx_get_schema("hierarchy")
  if (!is.null(h_path) && length(target_cats) > 0) {
    hier <- data.table::fread(h_path, select = c("parent_id", "child_id"))
    all_cats <- target_cats
    todo <- target_cats
    while(length(todo) > 0) {
      children <- hier$child_id[hier$parent_id %in% todo]
      new_children <- setdiff(children, all_cats)
      if (length(new_children) == 0) break
      all_cats <- c(all_cats, new_children)
      todo <- new_children
    }
    target_cats <- all_cats
  }
  
  # 2. Match fields in Schema
  f_path <- dx_get_schema("field")
  if (!is.null(f_path)) {
    field_schema <- data.table::fread(f_path, select = c("field_id", "main_category"))
    matched <- field_schema$field_id[field_schema$main_category %in% target_cats]
    if (length(matched) > 0) {
      leo.basic::leo_log("Found {length(matched)} fields in {length(category_ids)} categories via Official Schema.", level = "success")
      return(unique(as.character(matched)))
    }
  }
  
  return(character(0))
}

