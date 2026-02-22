#' Robust dx runner
#'
#' Executes dx CLI commands with automatic path detection and logging.
#' Supports macOS Homebrew installations and dry-run mode for testing.
#'
#' @param args Character vector of dx command arguments (e.g., `c("describe", "file-123")`).
#' @param intern Logical. If TRUE, capture stdout as character vector. If FALSE, return exit code.
#' @param ignore.stderr Logical. If TRUE, suppress stderr output.
#' @param dry_run Logical. If TRUE, only logs the command without executing it. Returns mock success.
#' @param verbose Logical. If TRUE (default), logs the dx command before execution.
#' @param ... Additional arguments passed to `system2()`.
#' @return If `intern=TRUE`, returns character vector of stdout. If `intern=FALSE`, returns exit code.
#'   In dry_run mode, returns `character(0)` if intern=TRUE, or `0L` if intern=FALSE.
#' @export
#' @examples
#' \dontrun{
#'   # Get current project info
#'   env_info <- dx_run(c("env"), intern = TRUE)
#'
#'   # Extract project ID from environment
#'   # Output format: "Current workspace\tproject-XXX"
#'   workspace_line <- grep("Current workspace\\t", env_info, value = TRUE)
#'   project_id <- sub(".*\\t(project-[^\\s]+).*", "\\1", workspace_line)
#'
#'   # Now use the extracted project ID
#'   result <- dx_run(c("describe", project_id), intern = TRUE)
#'
#'   # Find datasets in current project
#'   datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE)
#'
#'   # Upload file without capturing output (dry run)
#'   exit_code <- dx_run(c("upload", "file.txt", "--path", "/data/"), intern = FALSE, dry_run = TRUE)
#'
#'   # Test command without execution (dry run)
#'   dx_run(c("run", "app-table-exporter", "-y"), dry_run = TRUE)
#'
#'   # Ignore stderr for commands that may produce warnings
#'   dx_run(c("find", "data", "--brief"), ignore.stderr = TRUE)
#' }
dx_run <- function(args, intern = TRUE, ignore.stderr = FALSE, dry_run = FALSE, verbose = TRUE, ...) {
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

  if (dry_run) {
    if (verbose) {
      leo.basic::leo_log("[DRY RUN] Would execute: dx {paste(pretty_args, collapse = ' ')}", level = "warning")
    }
    # Return mock success for dry run
    return(if (intern) character(0) else 0L)
  }

  if (verbose) {
    leo.basic::leo_log("Executing: dx {paste(pretty_args, collapse = ' ')}")
  }

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

#' Get Dataset Dictionary
#'
#' Downloads the full data dictionary for a specific dataset using `dx extract_dataset -ddd`.
#' This file contains the authoritative list of column names (e.g., `p21003_i0`).
#' Uses caching to avoid repeated downloads.
#'
#' @param dataset The dataset ID or name (e.g. "app123.dataset" or "project-X:record-Y").
#' @param dry_run Logical. If TRUE, simulates the download without executing. Returns mock data.
#' @return A data.frame containing the dictionary with columns: `name`, `title`, etc.
#'   Returns NULL if download fails.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Find dataset record in current project
#'   # Step 1: Find all datasets in the project
#'   datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE)
#'
#'   # Step 2: Use the first dataset found (or select specific one)
#'   dataset_id <- datasets[1]  # e.g., "project-XXX:record-YYY"
#'
#'   # Step 3: Download and load dictionary for that dataset
#'   dict <- dx_get_dataset_dictionary(dataset_id)
#'   head(dict$name)  # Show first few column names
#'
#'   # Uses cached version if available (checks tmp/ directory)
#'   dict2 <- dx_get_dataset_dictionary(dataset_id)
#' }
dx_get_dataset_dictionary <- function(dataset) {
  # Determine storage location: prefer "tmp" in project, otherwise tempdir()
  storage_dir <- if (dir.exists("tmp")) "tmp" else tempdir()

  # Step 1: Resolve dataset ID if needed
  if (missing(dataset) || is.null(dataset)) stop("Dataset ID is required.")
  leo.basic::leo_log("Checking for info on {dataset}...")
  # Get dataset name to predict filename
  # Command: dx describe <id> --json
  # We parse the "name" field.
  # "Name" is usually the first line or labeled "Name".
  sys_cmd_args <- c("describe", dataset, "--json")
  desc_json <- tryCatch({
    dx_run(sys_cmd_args, intern = TRUE, ignore.stderr = TRUE)
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

  leo.basic::leo_log("Fetching dictionary for {dataset}...")

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
      leo.basic::leo_log("Output: {paste(res, collapse = '\\n')}")
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
  dict_df <- data.table::fread(target_file, data.table = FALSE)
  return(dict_df)
}

#' Get UKB Dataset Schema (Metadata)
#'
#' Fetches official UK Biobank schema files (field.tsv, category.tsv, hierarchy.tsv).
#' Prioritizes local cache, then DNAnexus project, then public UKB Showcase.
#'
#' @param type Schema type: "field", "category", or "hierarchy"
#'   - **"field"**: Field Schema (Schema 1) - Contains `field_id` and `main_category` mappings
#'   - **"category"**: Category Schema (Schema 3) - Category metadata (currently unused)
#'   - **"hierarchy"**: Hierarchy Schema (Schema 13) - Parent-child relationships between categories
#' @param force Logical. If TRUE, re-downloads even if cached locally.
#' @param dry_run Logical. If TRUE, simulates download without executing.
#' @return Path to the downloaded TSV file, or NULL if download fails.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Fetch Field Schema (Schema 1) - maps field IDs to categories
#'   field_path <- dx_get_schema("field")
#'   field_df <- data.table::fread(field_path)
#'   head(field_df[, c("field_id", "main_category")])
#'
#'   # Fetch Hierarchy Schema (Schema 13) - category parent-child relationships
#'   hier_path <- dx_get_schema("hierarchy")
#'   hier_df <- data.table::fread(hier_path)
#'   # Find all children of category 1712
#'   children <- hier_df$child_id[hier_df$parent_id == 1712]
#'
#'   # Force re-download (ignore cache)
#'   field_path <- dx_get_schema("field", force = TRUE)
#'
#'   # Test without actual download
#'   dx_get_schema("field", dry_run = TRUE)
#' }
dx_get_schema <- function(type = c("field", "category", "hierarchy"), force = FALSE, dry_run = FALSE) {
  type <- match.arg(type)

  # Schema IDs: 1=Field, 3=Category, 13=Hierarchy
  schema_id <- switch(type, "field" = 1, "category" = 3, "hierarchy" = 13)
  filename <- paste0(type, ".tsv")

  # Define paths
  if (!dir.exists("tmp")) dir.create("tmp")
  storage_dir <- "tmp"
  local_path <- file.path(storage_dir, filename)

  # 1. Check Local Cache
  if (!force && file.exists(local_path) && file.size(local_path) > 0) {
    # Optional: check age? For now assume stable.
    leo.basic::leo_log("Using cached schema: {local_path}")
    return(local_path)
  }

  if (dry_run) {
    leo.basic::leo_log("[DRY RUN] Would download schema: {filename}", level = "warning")
    leo.basic::leo_log("[DRY RUN] Would try: dx download Showcase metadata/{filename} -o {local_path} -f", level = "warning")
    leo.basic::leo_log("[DRY RUN] Or fallback to: {paste0('https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id=', schema_id)}", level = "warning")
    # Return existing file if available, otherwise return NULL
    if (file.exists(local_path)) {
      return(local_path)
    } else {
      return(NULL)
    }
  }

  # 2. Try DNAnexus (official location in RAP)
  # Path: "/Showcase metadata/<type>.tsv"
  leo.basic::leo_log("Looking for {filename}...")
  rap_path <- paste0("Showcase metadata/", filename)
  res <- tryCatch({
      dx_run(c("download", rap_path, "-o", local_path, "-f"), intern = TRUE, ignore.stderr = TRUE, dry_run = dry_run)
      TRUE
  }, error = function(e) FALSE)

  if (file.exists(local_path) && file.size(local_path) > 0) {
      leo.basic::leo_log("Downloaded {filename} from RAP project.", level = "success")
      return(local_path)
  }
  # 3. Fallback: Download from UK Biobank Public Showcase
  # URL: https://biobank.ndph.ox.ac.uk/ukb/scdown.cgi?fmt=txt&id={1,2}
  leo.basic::leo_log("Metadata not found in project. Fetching from UKBiobank Showcase...")
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

#' Expand Field IDs to Column Names
#'
#' Uses the official UKB RAP regex strategy to expand field IDs into all actual column names.
#' Filters by entity to ensure only fields for the specified entity are returned.
#'
#' @param fields Vector of field IDs (e.g. `c("21003", "p31", "41270")`).
#'   Accepts both numeric IDs and p-prefixed IDs.
#' @param dictionary The dictionary dataframe from `dx_get_dataset_dictionary()`.
#'   Must contain a `name` column with actual column names.
#' @param entity The entity type to filter for (e.g. "participant", "hesin", "death").
#'   Default is "participant". Essential for filtering multi-entity dictionaries.
#' @return Character vector of valid column names found in the dictionary for the specified entity.
#'   Returns empty character vector if no matches found.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Get dictionary first (find dataset record in project)
#'   datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE)
#'   dict <- dx_get_dataset_dictionary(datasets[1])
#'
#'   # Expand field IDs to column names
#'   cols <- dx_expand_field(c("p31", "21003"), dict)
#'   # Returns: c("p31", "p21003_i0", "p21003_i1", ...)
#'
#'   # Expand multiple fields
#'   cols <- dx_expand_field(c("41270", "41271", "41280"), dict)
#'
#'   # Filter by different entity
#'   hesin_cols <- dx_expand_field(c("41270"), dict, entity = "hesin")
#' }
dx_expand_field <- function(fields, dictionary, entity = "participant") {
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

  fields_char <- as.character(fields)
  all_names <- dictionary$name
  valid_cols <- character(0)

  for (f in fields_char) {
    # Strategy 1: Direct name match (works for non-participant entities like gp_scripts, hesin, death, etc.)
    if (f %in% all_names) {
      valid_cols <- c(valid_cols, f)
      next
    }

    # Strategy 2: Numeric UKB field ID expansion (participant entity convention)
    # 21003 -> "21003", p31 -> "31"
    fid <- sub("^\\D*(\\d+).*", "\\1", f)
    if (fid == "" || is.na(fid)) {
      leo.basic::leo_log("Field '{f}' not found for entity '{entity}'.", level = "warning")
      next
    }

    # Official tutorial: https://github.com/UK-Biobank/UKB-RAP-Notebooks-Access/blob/main/RStudio/A110_Export_participant_data.Rmd
    regex <- paste0("^p", fid, "(?![0-9])")
    matches <- grep(regex, all_names, value = TRUE, perl = TRUE)
    if (length(matches) > 0) {
      valid_cols <- c(valid_cols, matches)
    } else {
      leo.basic::leo_log("Field '{f}' (p{fid}) not found for entity '{entity}'.", level = "warning")
    }
  }

  return(unique(valid_cols))
}

#' Expand Category IDs to Field IDs
#'
#' Scans the official schema to expand category IDs into all field IDs belonging to those categories.
#' Only accepts numeric category IDs, not category names.
#'
#' **Logic**:
#' 1. Takes input category ID(s) (e.g., 1712)
#' 2. Recursively expands via hierarchy.tsv to find all child categories
#'    (e.g., 1712 → 2401, 2403, 2404, ...)
#' 3. Matches these categories against field.tsv's main_category column
#' 4. Returns all field_ids whose main_category matches any expanded category
#'
#' @param category_ids Vector of numeric Category IDs (e.g. `1712`, `c(1712, 100)`).
#'   Must be integers. String names like "First Occurrence" are NOT supported.
#' @param title_pattern_to_keep Optional regex. If provided, only keep fields
#'   whose `title` matches this pattern (case-insensitive).
#' @param title_pattern_to_exclude Optional regex. If provided, drop fields
#'   whose `title` matches this pattern (case-insensitive). Applied after keep.
#' @param entity Entity type to filter for (default: "participant").
#'   Currently unused but kept for future compatibility.
#' @return Character vector of field IDs (e.g. `c("41270", "41271", ...)`).
#'   Returns empty character vector if no fields found.
#' @keywords internal
#' @noRd
#' @examples
#' \dontrun{
#'   # Expand Category 1712 (First occurrences) to field IDs
#'   # Automatically downloads hierarchy and field schemas if missing
#'   fields_1712 <- dx_expand_field_by_category(1712)
#'   length(fields_1712)  # ~2330 fields
#'
#'   # Filter by title patterns (First occurrences)
#'   fields_no_source <- dx_expand_field_by_category(1712, title_pattern_to_exclude = "^Source")
#'   fields_date <- dx_expand_field_by_category(1712, title_pattern_to_keep = "^Date")
#'   setequal(fields_no_source, fields_date)
#'
#'   # Expand multiple categories
#'   fields_mix <- dx_expand_field_by_category(c(1712, 100))
#'   fields_mix_no_source <- dx_expand_field_by_category(c(1712, 100), title_pattern_to_exclude = "Source")
#' }
#' @importFrom data.table fread
dx_expand_field_by_category <- function(category_ids, 
                                        title_pattern_to_keep = NULL,
                                        title_pattern_to_exclude = NULL,
                                        entity = "participant") {
  if (is.null(category_ids) || length(category_ids) == 0) return(character(0))

  # 1. Expand categories via Hierarchy (resolve parents to children)
  target_cats <- suppressWarnings(as.integer(category_ids))
  target_cats <- target_cats[!is.na(target_cats)]

  if (length(target_cats) == 0) {
    leo.basic::leo_log("Invalid category_id format. Category IDs must be numeric integers (e.g., 1712, not 'First Occurrence').", level = "danger")
    return(character(0))
  }

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
    leo.basic::leo_log("Expanded {length(category_ids)} categories to {length(target_cats)} categories (including subcategories)")
  }

  # 2. Match fields in Schema
  f_path <- dx_get_schema("field")
  if (!is.null(f_path)) {
    field_schema <- data.table::fread(f_path, select = c("field_id", "main_category", "title"))
    matched_rows <- field_schema[field_schema$main_category %in% target_cats, .(field_id, title)]
    matched <- matched_rows$field_id
    if (length(matched) > 0) {
      leo.basic::leo_log("Found {length(matched)} fields in {length(category_ids)} categories via Official Schema.", level = "success")

      # Apply title filters if provided
      if (!is.null(title_pattern_to_keep) && nzchar(title_pattern_to_keep)) {
        keep_mask <- grepl(title_pattern_to_keep, matched_rows$title, ignore.case = TRUE)
        keep_mask[is.na(keep_mask)] <- FALSE
        matched_rows <- matched_rows[keep_mask, ]
        matched <- matched_rows$field_id
        leo.basic::leo_log("After applying title_pattern_to_keep filter, {length(matched)} fields remain.", level = "success")
      }
      if (!is.null(title_pattern_to_exclude) && nzchar(title_pattern_to_exclude)) {
        exclude_mask <- grepl(title_pattern_to_exclude, matched_rows$title, ignore.case = TRUE)
        exclude_mask[is.na(exclude_mask)] <- FALSE
        matched_rows <- matched_rows[!exclude_mask, ]
        matched <- matched_rows$field_id
        leo.basic::leo_log("After applying title_pattern_to_exclude filter, {length(matched)} fields remain.", level = "success")
      }
      return(unique(as.character(matched)))
    } else {
      leo.basic::leo_log("No fields found for category IDs: {paste(category_ids, collapse=', ')}", level = "warning")
    }
  }
  return(character(0))
}
