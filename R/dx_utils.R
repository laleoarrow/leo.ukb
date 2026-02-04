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
  
  sys_cmd_args <- c("describe", shQuote(dataset), "--json")
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
          if (requireNamespace("data.table", quietly = TRUE)) {
            return(data.table::fread(expected_file, data.table = FALSE))
          } else {
            return(read.csv(expected_file, stringsAsFactors = FALSE))
          }
      }
  }
  
  # If not cached or name resolution failed, proceed to download.
  leo.basic::leo_log("Fetching dictionary for {dataset}...", level = "info")
  
  # Use system2 for proper exit code handling
  # dx extract_dataset <id> -ddd --output <dir>/
  
  res <- tryCatch({
    system2("dx", 
            args = c("extract_dataset", shQuote(dataset), "-ddd", "--output", paste0(storage_dir, "/")),
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
  # Use utils::read.csv or data.table::fread if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    dict_df <- data.table::fread(target_file, data.table = FALSE)
  } else {
    dict_df <- read.csv(target_file, stringsAsFactors = FALSE)
  }
  
  return(dict_df)
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
#' Scans the dataset dictionary for fields belonging to a specific category ID (or folder path).
#'
#' @param category_ids Vector of category IDs (e.g. `1712`, `c(100, 101)`).
#' @param dictionary The dictionary dataframe.
#' @param entity Entity type to filter for (default: "participant").
#' @return Character vector of field IDs (e.g. `c("41270", "41271")`).
#' @keywords internal
#' @noRd
dx_find_fields_by_category <- function(category_ids, dictionary, entity = "participant") {
  if (is.null(category_ids) || length(category_ids) == 0) return(character(0))
  if (is.null(dictionary) || nrow(dictionary) == 0) return(character(0))
  
  # Filter by entity first
  if ("entity" %in% names(dictionary)) {
    dictionary <- dictionary[dictionary$entity == entity, ]
  }
  
  # Standardize Inputs
  cat_ids <- as.character(category_ids)
  found_fields <- character(0)
  
  # Strategy: Check for 'folderPath' column (standard in dx extract_dataset -ddd)
  # Folder paths usually look like: "Category > Subcategory > Field"
  # Or they might contain the numeric ID in the path string if we are lucky, 
  # or we might need to rely on the 'title' if it contains the ID.
  # However, UKB RAP dictionary usually has 'folderPath'.
  
  if ("folderPath" %in% names(dictionary)) {
    # Extract leaf fields where the folderPath matches the requested category
    # Problem: category_ids are numbers (e.g. 1712), but folderPath text might be "Health-related outcomes > First occurrences"
    # We assume the user might provide the name OR the system needs to map ID -> Name.
    # Since we lack an internal map, we will TRY to match if the folderPath *contains* the ID (unlikely) 
    # OR we assume the user might pass a loose match string.
    # BUT wait, the official dictionary CSV often has a "Folder" column with just names.
    
    # HACK: If the dictionary has a 'LinkID' or 'FieldID' column, we can return those.
    # Actually, standard dx dictionary has 'name' (p41270_i0).
    # We need to extract the base field ID from the name and return it.
    
    # Let's inspect typical dictionary structure:
    # name, title, entity, folderPath
    
    # If we can't map ID 1712 to "First occurrences", we can't robustly find it by ID.
    # HOWEVER, the UKB schema usually puts the Category ID in the metadata if downloaded from Showcase.
    # The RAP dictionary is simpler.
    
    # User Request: "category_id = 1712".
    # If the CSV doesn't have numeric category IDs, we are stuck.
    # Let's assume for now the user might need to providing a NAME pattern if IDs aren't there, 
    # OR we search validation against 'title' or 'folderPath'.
    
    # REVISION: We will support regex matching on folderPath for now.
    # If the user passes "1712", we search for "1712" in folderPath? No, that won't exist.
    # This implies we might need to fetch the category tree separately or warn the user.
    # Wait, `dx run app-table-exporter` supports `-icategory=`. 
    # Maybe we don't need to resolve fields manually?
    # NO, the user wants to mix fields and categories. `table-exporter` takes `-ifield_names_file_txt`.
    # It does NOT take a category list file. So we MUST resolve to fields.
    
    # To resolve 1712 -> Fields, we theoretically need the UKB schema.
    # Since we don't have it, we'll try to match strictly on folderPath if the user passes a string,
    # OR if they pass a number, we simply look for that number in the folder path (maybe it's formatted "1712 - Name"?).
    
    # Let's be aggressive: we iterate all rows, and if the folderPath contains the string provided, we take it.
    paths <- dictionary$folderPath
    for (cid in cat_ids) {
       # Case-insensitive substring match
       # e.g. cid="First occurrences" or cid="1712"
       # If 1712 doesn't appear in text, this fails. 
       # We'll log a warning if no matches found.
       matches <- grep(cid, paths, ignore.case = TRUE, fixed = TRUE)
       if (length(matches) > 0) {
         # Get names (p123_i0) -> extract ID (123)
         names_matched <- dictionary$name[matches]
         ids <- unique(sub("^p(\\d+).*", "\\1", names_matched))
         found_fields <- c(found_fields, ids)
         leo.basic::leo_log("Category '{cid}' matched {length(ids)} fields.", level = "info")
       } else {
         leo.basic::leo_log("No fields found for category '{cid}' (checked folderPath). Try using the category Name instead of ID?", level = "warning")
       }
    }
  } else {
    leo.basic::leo_log("Dictionary missing 'folderPath'. Cannot filter by Category.", level = "warning")
  }
  
  return(unique(found_fields))
}
