#' Log in to DNAnexus
#'
#' This is a wrapper for \code{dx login}. It supports both interactive login
#' (if no token provided) and token-based login.
#'
#' @param token Optional. A DNAnexus authentication token.
#' @return None (runs system command).
#' @export
dx_login <- function(token = NULL) {
  t0 <- Sys.time()

  # Get dx path for system() call
  dx_path <- Sys.which("dx")
  if (dx_path == "" && Sys.info()["sysname"] == "Darwin") {
    for (p in c("/opt/homebrew/bin/dx", "/usr/local/bin/dx")) {
      if (file.exists(p)) { dx_path <- p; break }
    }
  }
  if (dx_path == "") dx_path <- "dx"

  if (!is.null(token)) {
    # Token-based login (non-interactive)
    .dx_run(c("login", "--token", token), intern = FALSE)
  } else {
    # Interactive login: system() handles TTY prompts (keyboard input)
    # much better than system2() in some R consoles
    system(paste(shQuote(dx_path), "login"))
  }

  leo.basic::leo_time_elapsed(t0)
}

#' Extract UKB data using DNAnexus Table Exporter
#'
#' This function takes a file containing UKB field IDs and runs the Table Exporter
#' app on DNAnexus RAP to extract the data. The output will be saved to `data_files/` on the RAP project. Ref (https://github.com/UK-Biobank/UKB-RAP-Notebooks-Access/blob/main/RStudio/A110_Export_participant_data.Rmd)
#'
#' @param field_id Path to a file or a vector of Field IDs (e.g., `c("p31", "41270")`).
#' @param category_id Path to a file or a vector of numeric Category IDs (e.g., `1712`, `c(1712, 100)`).
#'   Category IDs must be integers and correspond to UKB data categories in the hierarchy schema.
#' @param output_prefix Output file prefix. If NULL, auto-generated based on fields/categories.
#' @param title_pattern_to_keep Optional regex to keep category-expanded fields by `title`
#'   (case-insensitive). Only applies when `category_id` is provided.
#' @param title_pattern_to_exclude Optional regex to drop category-expanded fields by `title`
#'   (case-insensitive). Applied after `title_pattern_to_keep`.
#' @param expand Logical. If TRUE (default), validates and expands field IDs (e.g. `21003`)
#'   into all instances/arrays present in the dataset (e.g. `p21003_i0`, `p21003_i1`)
#'   using the official dataset dictionary verification strategy.
#' @param entity Entity type. Options:
#'   \itemize{
#'     \item \code{"participant"} (default): Main phenotype/covariate data.
#'     \item \code{"hesin"}: Hospital episode statistics.
#'     \item \code{"death"}: Death registry records.
#'     \item \code{"image"}: Imaging metadata.
#'   }
#' @param coding_option How to encode categorical values. Options:
#'   \itemize{
#'     \item \code{"RAW"} (default): Return raw numeric codes (e.g., 1, 2, 3).
#'     \item \code{"REPLACE"}: Replace codes with human-readable labels (e.g., "Male", "Female").
#'   }
#' @param output_format Output file format. Options:
#'   \itemize{
#'     \item \code{"CSV"} (default): Comma-separated values.
#'     \item \code{"TSV"}: Tab-separated values.
#'   }
#' @param header_style Column header naming style. Options:
#'   \itemize{
#'     \item \code{"FIELD-NAME"} (default): Use field names like \code{p21003_i0}.
#'     \item \code{"FIELD-TITLE"}: Use descriptive titles like \code{"Age when attended assessment centre"}.
#'     \item \code{"UKB-FORMAT"}: Use raw UKB field IDs like \code{21003-0.0}.
#'   }
#' @param instance_type The DNAnexus instance type for the job. Default: \code{"mem1_ssd1_v2_x4"}.
#'   For large extracts, \code{"mem1_hdd1_v2_x8"} is recommended.
#' @param dataset The dataset ID/name to extract from (e.g., \code{"app12345_20240101.dataset"} or \code{"project-Gk2...:record-Fp3..."}).
#'   Defaults to the latest dataset in the project if NULL.
#' @param dry_run Logical. If TRUE, only skips the Table Exporter job submission.
#'   All other steps (dataset detection, dictionary validation, file upload) still run.
#'   Useful for testing command generation without launching jobs. Default is FALSE.
#' @return A list containing the job ID and expected output path (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' dx_login() # Make sure you have logged in to DNAnexus
#' dx_status() # Check the status of your login
#'
#' # 1. Basic extraction (extracts covariates defined in fields.txt)
#' dx_extract("tmp/fields.txt", output_prefix = "ukb_data", dry_run = TRUE)
#'
#' # 2. Extraction using vectors (Field & Category)
#' dx_extract(field_id = c("p53", "p191", "p40000", # basic demographics
#'                         "p41270", "p41271", "p41280", "p41281"), # ICD codes
#'            category_id = 1712, # First occurrences
#'            title_pattern_to_exclude = "^Source",
#'            output_prefix = "dignosis20260205")
#'
#' # 3. Just Category (numeric ID or vector of IDs)
#' dx_extract(category_id = 1712)  # Extract all fields in category 1712
#' dx_extract(category_id = c(1712, 100))  # Multiple categories
#'
#' # 3b. First occurrences (FO) title filters
#' # Keep only titles starting with "Date"
#' dx_extract(category_id = 1712, title_pattern_to_keep = "^Date")
#' # Exclude titles starting with "Source"
#' dx_extract(category_id = 1712, title_pattern_to_exclude = "^Source")
#'
#' # 4. Extract with human-readable labels and descriptive headers
#' dx_extract("tmp/fields.txt",
#'                output_prefix = "ukb_data_readable",
#'                coding_option = "REPLACE",
#'                header_style = "FIELD-TITLE")
#'
#' # 5. Explicitly specifying the dataset (Recommended for reproducibility)
#' dx_extract("tmp/fields.txt",
#'                dataset = "project-Gk2PzX0Jj1X4Y5Z6:record-Fp37890Qj9k2X1Y4")
#'
#' # 6. Dry run mode - Test without actually submitting the job
#' # Useful for debugging and verifying the command before execution
#' result <- dx_extract(
#'   field_id = c("p31", "p21003"),
#'   category_id = 1712,
#'   output_prefix = "test_run",
#'   dry_run = TRUE
#' )
#' # Shows all commands that would be executed
#' # Returns: list(job_id = "job-DRYRUN123456", output_path = "...", dry_run = TRUE)
#' }
dx_extract <- function(field_id = NULL,
                       category_id = NULL,
                       output_prefix = NULL,
                       title_pattern_to_keep = NULL,
                       title_pattern_to_exclude = NULL,
                       expand = TRUE,
                       entity = "participant",
                       coding_option = "RAW",
                       output_format = "CSV",
                       header_style = "FIELD-NAME",
                       instance_type = "mem1_ssd1_v2_x4",
                       dataset = NULL,
                       dry_run = FALSE
) {
  # Input validation
  # Support alias args for backward compatibility (though this is a breaking change version)
  # We focus on the new signature: field_id, category_id.

  if (is.null(field_id) && is.null(category_id)) {
    leo.basic::leo_log("At least one of 'field_id' or 'category_id' must be provided.", level = "danger")
    return(invisible(NULL))
  }

  process_input <- function(input, name) {
    if (is.null(input)) return(character(0))
    # Check if input is a single string that looks like a file path existing on disk
    if (length(input) == 1 && is.character(input) && file.exists(as.character(input))) {
      leo.basic::leo_log("Reading {name} from file: {input}")
      lines <- readLines(input, warn = FALSE)
      return(lines[trimws(lines) != ""])
    }
    # Otherwise treat as vector
    return(as.character(input))
  }

  fields <- process_input(field_id, "fields")
  cats <- process_input(category_id, "categories")

  # Resolve categories if provided
  if (length(cats) > 0) {
    # Validate that all category IDs are numeric
    cat_numeric <- suppressWarnings(as.numeric(cats))
    invalid_cats <- cats[is.na(cat_numeric)]

    if (length(invalid_cats) > 0) {
      leo.basic::leo_log("Invalid category IDs detected: {paste(invalid_cats, collapse=', ')}. Category IDs must be numeric (e.g., 1712, not 'First Occurrence').", level = "danger")
      return(invisible(NULL))
    }

    leo.basic::leo_log("Resolving {length(cats)} categories...")
    # dx_expand_field_by_category uses Official Schema (no dictionary needed)
    cat_fields <- dx_expand_field_by_category(
      cats,
      title_pattern_to_keep = title_pattern_to_keep,
      title_pattern_to_exclude = title_pattern_to_exclude,
      entity = entity
    )

    if (length(cat_fields) == 0) {
      leo.basic::leo_log("No fields found for the specified category IDs. Please verify the category IDs are correct.", level = "warning")
    }

    # Append found fields to the main list
    fields <- c(fields, cat_fields)
  }

  cli::cat_rule("Extracting Data using DNAnexus Table Exporter", col = "blue")

  # Display dx environment info
  tryCatch({
    dx_env <- .dx_run("env", intern = TRUE, ignore.stderr = TRUE)
    user_line <- grep("Current user", dx_env, value = TRUE)
    project_line <- grep("Current workspace name", dx_env, value = TRUE)
    if (length(user_line) > 0) leo.basic::leo_log(trimws(user_line))
    if (length(project_line) > 0) leo.basic::leo_log(trimws(project_line))
  }, error = function(e) {
    leo.basic::leo_log("Could not retrieve dx environment info", level = "warning")
  })

  # Auto-detect latest dataset if not specified
  if (is.null(dataset)) {
    datasets <- tryCatch({
      .dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE, ignore.stderr = TRUE)
    }, error = function(e) character(0))

    if (length(datasets) > 0) {
      # Take the first (or latest by timestamp if names include dates)
      # dx find data returns newest first by default
      dataset <- datasets[1]
      leo.basic::leo_log("Auto-detected dataset: {dataset}")
    } else {
      leo.basic::leo_log("No dataset found in project. Please specify the dataset parameter.", level = "danger")
      return(invisible(NULL))
    }
  }

  fields <- trimws(fields)
  fields <- fields[fields != ""]

  if (length(fields) == 0) {
    leo.basic::leo_log("No valid field IDs found.", level = "danger")
    return(invisible(NULL))
  }

  # Ensure eid is always the first field
  fields <- fields[!tolower(fields) %in% c("eid", "p_eid")]
  fields <- unique(fields)  # Remove duplicates
  fields <- c("eid", fields)
  leo.basic::leo_log("Loaded {length(fields)} field IDs (eid + {length(fields) - 1} fields)")

  # Expand fields
  # New Logic: Official Dataset Dictionary Validation - we fetch the actual dataset dictionary and filter.
  if (expand) {
    # If we haven't fetched dictionary yet (only fields provided), do it now
    if (!exists("ukb_dd")) {
      leo.basic::leo_log("Validating fields against dataset dictionary ({dataset})...")
      ukb_dd <- dx_get_dataset_dictionary(dataset, dry_run = FALSE)
    }

    if (!is.null(ukb_dd)) {
      # 2. Expand field IDs to actual column names (Official Regex), filtered by entity
      # Exclude 'eid' from lookup as it's added manually later
      valid_cols <- dx_expand_field(fields[fields != "eid"], ukb_dd, entity = entity)

      if (length(valid_cols) > 0) {
        fields <- c("eid", valid_cols)
        leo.basic::leo_log("Validated and expanded to {length(fields)} existing columns.", level = "success")
      } else {
        leo.basic::leo_log("No valid columns found for the requested fields in this dataset!", level = "danger")
        return(invisible(NULL))
      }
    } else {
      leo.basic::leo_log("Could not retrieve dataset dictionary. Skipping validation (risky).", level = "warning")
    }
  }

  # Generate output prefix provided
  if (is.null(output_prefix)) {
    # Smart naming based on inputs
    # Logic:
    # 1. Gather all explicit inputs: categories (c) and original fields (f)
    # 2. If distinct items < 3, SHOW ALL (e.g. c1712_f41270_f41271)
    # 3. If string length <= 20, SHOW ALL
    # 4. Otherwise use summary: c@{num_cats}_f@{num_fields}

    # Clean inputs for naming
    # Remove 'p' prefix if present for cleaner names
    clean_fields <- sub("^p", "", unique(process_input(field_id, "")))
    clean_cats <- unique(cats)

    # Build candidate parts
    parts <- character(0)
    if (length(clean_cats) > 0) parts <- c(parts, paste0("c", clean_cats))
    if (length(clean_fields) > 0) parts <- c(parts, paste0("f", clean_fields))

    candidate_base <- paste(parts, collapse = "_")
    total_items <- length(clean_cats) + length(clean_fields)

    final_base <- ""
    if (total_items < 3) {
      final_base <- candidate_base
    } else if (nchar(candidate_base) <= 20) {
      final_base <- candidate_base
    } else {
      # Fallback summary
      # Note: If 0 fields provided, we omit f@0 to be clean? — If 0, then do not use f@0, just put c@xxx
      # Or user wants consistent f@M_c@N format?
      # "use f@[fieldID count]_c@[cate count]"
      sum_parts <- character(0)
      if (length(clean_fields) > 0) sum_parts <- c(sum_parts, paste0("f@", length(clean_fields)))
      if (length(clean_cats) > 0) sum_parts <- c(sum_parts, paste0("c@", length(clean_cats)))
      final_base <- paste(sum_parts, collapse = "_")
    }

    # Always append total column count
    output_prefix <- paste0("extract_", final_base, "_n", length(fields))
  }

  # Save fields to file (locally)
  # If tmp/ exists, save there for user inspection (and do not delete).
  # Otherwise use tempdir().

  if (dir.exists("tmp")) {
    fields_file <- file.path("tmp", paste0(output_prefix, "_fields.txt"))
    writeLines(fields, fields_file)
    leo.basic::leo_log("Saved expanded fields locally: {fields_file}", level = "success")
    # Do NOT unlink if in project tmp
  } else {
    fields_file <- tempfile(fileext = ".txt")
    writeLines(fields, fields_file)
    on.exit(unlink(fields_file), add = TRUE)
  }

  # Upload the field names file to fields_files/ folder
  dx_file_name <- paste0(output_prefix, "_fields.txt")
  dx_file_path <- paste0("/fields_files/", dx_file_name)
  leo.basic::leo_log("Uploading field names file to: {dx_file_path}")

  # Delete existing file if present (to avoid duplicates)
  .dx_run(c("rm", "-f", dx_file_path), ignore.stderr = TRUE, dry_run = FALSE)

  upload_exit <- tryCatch({
    res <- .dx_run(c("upload", normalizePath(fields_file, mustWork = TRUE), "--path", dx_file_path, "--parents"), intern = FALSE, ignore.stderr = TRUE, dry_run = FALSE)
    if (is.null(res)) 1L else res
  }, error = function(e) {
    leo.basic::leo_log(paste0("Upload error: ", e$message), level = "danger")
    return(1L)
  })

  if (upload_exit != 0L && !dry_run) {
    leo.basic::leo_log("Failed to upload field names file", level = "danger")
    return(invisible(NULL))
  }

  leo.basic::leo_log("Upload successful", level = "success")

  # Run Table Exporter app
  # dx run app-table-exporter -ientity=participant -ioutput=<prefix> -icoding_option=RAW
  #   -idataset_or_cohort_or_dashboard=<dataset> -ifield_names_file_txt=<file>
  #   -ioutput_format=CSV -iheader_style=FIELD-NAME --destination data_files/
  args <- c(
    "run", "app-table-exporter",
    paste0("-ientity=", entity),
    paste0("-ioutput=", output_prefix),
    paste0("-icoding_option=", coding_option),
    paste0("-idataset_or_cohort_or_dashboard=", dataset),
    paste0("-ifield_names_file_txt=", dx_file_path),
    paste0("-ioutput_format=", output_format),
    paste0("-iheader_style=", header_style),
    "--destination", "/data_files/",
    "--priority", "high",
    "--instance-type", instance_type,
    "-y"  # auto-confirm
  )

  if (dry_run) {
    leo.basic::leo_log("[DRY RUN] Would run Table Exporter with command:", level = "warning")
    leo.basic::leo_log("dx {paste(args, collapse = ' ')}")
    leo.basic::leo_log("[DRY RUN] Skipping actual job submission", level = "warning")
    return(invisible(list(
      job_id = "job-DRYRUN123456",
      output_path = paste0("project:/data_files/", output_prefix, ".", tolower(output_format)),
      dry_run = TRUE
    )))
  }

  leo.basic::leo_log("Running Table Exporter...")
  leo.basic::leo_log("dx {paste(args, collapse = ' ')}")

  res <- .dx_run(args, intern = TRUE, ignore.stderr = FALSE, dry_run = dry_run)
  exit_code <- attr(res, "status")
  if (is.null(exit_code)) exit_code <- 0L

  if (exit_code != 0L) {
    leo.basic::leo_log("Table Exporter failed with exit code {exit_code}", level = "danger")
    leo.basic::leo_log("Output: {paste(res, collapse = '\\n')}")
    return(invisible(NULL))
  }

  # Parse job ID from output (format: "Job ID: job-XXXX")
  job_line <- grep("Job ID:", res, value = TRUE)
  job_id <- if (length(job_line) > 0) {
    sub(".*Job ID:\\s*", "", job_line[1])
  } else {
    NA_character_
  }

  output_path <- paste0("project:/data_files/", output_prefix, ".", tolower(output_format))

  leo.basic::leo_log("Job submitted: {job_id}", level = "success")
  leo.basic::leo_log("Output will be saved to: {output_path}")

  return(invisible(list(job_id = job_id, output_path = output_path)))
}

#' Check DNAnexus Job Status & Dashboard
#'
#' Checks the status of a specific job or displays a dashboard of recent jobs.
#'
#' When called without a \code{job_id}, it displays a "htop"-style dashboard of the
#' last several jobs.
#'
#' @param all_projects Logical. If TRUE, shows jobs from all projects in the account (account-wide).
#'   Default is FALSE (only current project).
#' @param job_id Optional. A string containing the job ID (e.g., "job-Gk2...") or a list
#'   returned by \code{\link{dx_extract}}. If provided, \code{all_projects} is ignored.
#' @param limit Numeric. Number of recent jobs to show in the dashboard (default: 3).
#' @return The state of the job as a string (if \code{job_id} provided) or a data frame
#'   of recent jobs (if \code{job_id} is NULL), invisibly.
#' @export
#' @examples
#' \dontrun{
#' # 1. View dashboard of recent jobs in current project
#' dx_status()
#'
#' # 2. View dashboard for ONLY the current project
#' dx_status(all_projects = FALSE)
#'
#' # 3. Check status of a specific job
#' dx_status(job_id = "job-Gk2PzX0Jj1X4Y5Z6")
#' }
dx_status <- function(all_projects = TRUE, job_id = NULL, limit = 5) {
  t0 <- Sys.time()
  # --- 1. Single Job Mode (If job_id is provided) ---
  if (!is.null(job_id)) {
    if (is.list(job_id) && !is.null(job_id$job_id)) job_id <- job_id$job_id
    if (is.null(job_id) || is.na(job_id) || !is.character(job_id) || !grepl("^job-", job_id)) {
      leo.basic::leo_log("Invalid Job ID provided.", level = "danger")
      return(NA_character_)
    }

    cmd_args <- c("describe", job_id, "--json")
    res <- tryCatch({.dx_run(cmd_args, intern = TRUE, ignore.stderr = TRUE, verbose = FALSE)}, error = function(e) NULL)
    if (is.null(res) || length(res) == 0) return(NA_character_)

    json_str <- paste(res, collapse = " ")
    matches <- regmatches(json_str, regexec('"state"\\s*:\\s*"([^"]+)"', json_str))
    state <- if (length(matches[[1]]) > 1) matches[[1]][2] else NA_character_

    if (!is.na(state)) {
      color_map <- list("done"="success", "failed"="danger", "running"="info", "runnable"="warning", "terminated"="danger")
      lvl <- if (!is.null(color_map[[state]])) color_map[[state]] else "info"
      leo.basic::leo_log("Job {job_id} is: {state}", level = lvl)
    }
    return(state)
  }

  # --- 2. Dashboard Mode (If job_id is NULL) ---
  status_id <- cli::cli_status("Fetching DNAnexus data...")
  on.exit(cli::cli_status_clear(status_id), add = TRUE)

  # Get environment info (contains user and project)
  env_info <- tryCatch(.dx_run("env", intern = TRUE, ignore.stderr = TRUE, verbose = FALSE), error = function(e) NULL)
  if (is.null(env_info) || length(env_info) == 0) {
    leo.basic::leo_log("Could not get DNAnexus environment. Please run 'dx login' in terminal.", level = "danger")
    return(invisible(NULL))
  }

  # Extract user and current project info
  user_line <- grep("Current user", env_info, value = TRUE)
  user <- if (length(user_line) > 0) sub("^Current user\\s+", "", user_line[1]) else "Unknown"

  curr_proj_line <- grep("Current workspace name", env_info, value = TRUE)
  curr_proj_id_line <- grep("Current workspace", env_info, value = TRUE)
  curr_proj_name <- if (length(curr_proj_line) > 0) sub(".*Current workspace name\\s+", "", curr_proj_line[1]) else "Unknown"
  curr_proj_id <- if (length(curr_proj_id_line) > 0) {
    val <- curr_proj_id_line[grep("^Current workspace\\s+project-", curr_proj_id_line)]
    if(length(val) == 0) val <- curr_proj_id_line[grep("Current workspace", curr_proj_id_line)]
    sub(".*Current workspace\\s+", "", val[1])
  } else ""

  # Pre-cache current project
  if (curr_proj_id != "" && curr_proj_name != "Unknown" && exists(".dx_cache")) {
    assign(curr_proj_id, curr_proj_name, envir = .dx_cache)
  }

  cat(cli::style_italic(glue::glue(" User: {cli::col_yellow(user)} | Project: {cli::col_cyan(curr_proj_name)} ({curr_proj_id})")), "\n\n")

  # Fetch recent jobs
  cmd_args <- c("find", "jobs", if(all_projects) "--all-projects" else character(0), "--num-results", as.character(limit), "--json")
  res <- tryCatch({
    .dx_run(cmd_args, intern = TRUE, ignore.stderr = TRUE, verbose = FALSE)
  }, error = function(e) NULL)

  if (is.null(res) || length(res) == 0) {
    leo.basic::leo_log("Login success, but no jobs found.")
    return(invisible(NULL))
  }

  # Parse JSON
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jobs_raw <- jsonlite::fromJSON(paste(res, collapse = " "))
    if (length(jobs_raw) == 0) {
      leo.basic::leo_log("No jobs found.")
      return(invisible(NULL))
    }
    jobs <- as.data.frame(jobs_raw)
  } else {
    leo.basic::leo_log("jsonlite not installed. Reverting to base dx display.", level = "warning")
    .dx_run(cmd_args, intern = FALSE)
    return(invisible(NULL))
  }

  # Project Mapping for dashboard view
  proj_map <- list()
  unique_proj_ids <- unique(jobs$project)
  for (pid in unique_proj_ids) {
    if (pid == curr_proj_id) {
      proj_map[[pid]] <- curr_proj_name
    } else if (all_projects) {
      p_info <- tryCatch(.dx_run(c("describe", pid, "--json"), intern = TRUE, ignore.stderr = TRUE), error = function(e) NULL)
      if (!is.null(p_info)) {
        p_meta <- jsonlite::fromJSON(paste(p_info, collapse = " "))
        proj_map[[pid]] <- p_meta$name
      } else {
        proj_map[[pid]] <- pid
      }
    } else {
      proj_map[[pid]] <- pid
    }
  }

  # Timing & Metadata Processing
  now_ms <- as.numeric(Sys.time()) * 1000

  # Duration calculation
  # startedRunning, stoppedRunning are in ms
  jobs$start <- ifelse(is.na(jobs$startedRunning), jobs$created, jobs$startedRunning)
  jobs$end <- ifelse(is.na(jobs$stoppedRunning), now_ms, jobs$stoppedRunning)
  jobs$duration_ms <- jobs$end - jobs$start

  format_dur <- function(ms) {
    s <- ms / 1000
    if (s < 60) return(paste0(round(s), "s"))
    m <- s / 60
    if (m < 60) return(paste0(round(m), "m"))
    h <- m / 60
    return(paste0(round(h, 1), "h"))
  }

  jobs$dur_str <- sapply(jobs$duration_ms, format_dur)
  jobs$age_str <- sapply(difftime(Sys.time(), as.POSIXct(jobs$created/1000, origin="1970-01-01")), function(age) {
    val_m <- as.numeric(age, units="mins")
    if (val_m < 60) return(paste0(round(val_m), "m"))
    val_h <- as.numeric(age, units="hours")
    if (val_h < 24) return(paste0(round(val_h, 1), "h"))
    return(paste0(round(as.numeric(age, units="days")), "d"))
  })

  # Clear status right before printing the table
  cli::cli_status_clear(status_id)

  # Header
  cols <- c("Job ID", "Status", "Project", "Cost", "Owner", "Since", "Duration", "Name")

  # Custom alignment function for Multi-byte (Chinese) and ANSI strings
  align_ansi <- function(x, width, side = "left") {
    clean_x <- cli::ansi_strip(as.character(x))
    n <- nchar(clean_x, type = "width") # Use visual display width
    pad_len <- max(0, width - n)
    pad <- strrep(" ", pad_len)
    if (side == "left") paste0(x, pad) else paste0(pad, x)
  }

  # Unify Header and Data column widths
  w <- c(id=28, status=12, project=20, cost=8, owner=8, since=7, dur=8)

  header_str <- paste(
    align_ansi(cols[1], w["id"]),
    align_ansi(cols[2], w["status"]),
    align_ansi(cols[3], w["project"]),
    align_ansi(cols[4], w["cost"]),
    align_ansi(cols[5], w["owner"]),
    align_ansi(cols[6], w["since"]),
    align_ansi(cols[7], w["dur"]),
    cols[8],
    sep = " | "
  )
  cat(cli::style_bold(header_str), "\n")
  cat(paste(rep("-", 140), collapse = ""), "\n")

  for (i in seq_len(nrow(jobs))) {
    s <- jobs$state[i]
    st_color <- switch(s,
                       "done" = cli::col_green(paste0("✔ ", s)),
                       "failed" = cli::col_red(paste0("✖ ", s)),
                       "running" = cli::col_blue(paste0("● ", s)),
                       "runnable" = cli::col_yellow(paste0("○ ", s)),
                       "waiting" = cli::col_yellow(paste0("○ ", s)),
                       "terminated" = cli::col_red(paste0("⚑ ", s)),
                       s
    )

    owner <- sub("user-", "", jobs$launchedBy[i])
    pname <- if (!is.null(proj_map[[jobs$project[i]]])) proj_map[[jobs$project[i]]] else jobs$project[i]

    # Cost string processing
    cost_val <- if ("totalPrice" %in% names(jobs)) jobs$totalPrice[i] else 0

    # Determine currency symbol for this specific job (Defensive Logic)
    row_currency <- "£"
    # Case 1: Nested data frame (non-flattened fromJSON)
    if ("currency" %in% names(jobs) && is.data.frame(jobs$currency)) {
      sym <- jobs$currency$symbol[i]
      if (!is.null(sym) && !is.na(sym)) row_currency <- sym
    }
    # Case 2: Flattened column name
    else if ("currency.symbol" %in% names(jobs)) {
      sym <- jobs$currency.symbol[i]
      if (!is.null(sym) && !is.na(sym)) row_currency <- sym
    }

    if (s %in% c("running", "runnable", "waiting") && (is.na(cost_val) || cost_val == 0)) {
      cost_str <- "Pending"
      p_color <- cli::style_italic
    } else {
      if (is.na(cost_val)) cost_val <- 0
      cost_str <- sprintf("%s%.2f", row_currency, cost_val)
      p_color <- cli::col_green
    }

    cat(align_ansi(jobs$id[i], w["id"]), " | ", sep="")
    cat(align_ansi(st_color, w["status"]), " | ", sep="")
    cat(align_ansi(cli::ansi_strtrim(pname, w["project"]), w["project"]), " | ", sep="")
    cat(align_ansi(p_color(cost_str), w["cost"]), " | ", sep="")
    cat(align_ansi(cli::ansi_strtrim(owner, w["owner"]), w["owner"]), " | ", sep="")

    age_styled <- cli::style_italic(jobs$age_str[i])
    dur_styled <- cli::col_cyan(jobs$dur_str[i])

    cat(align_ansi(age_styled, w["since"]), " | ", sep="")
    cat(align_ansi(dur_styled, w["dur"]), " | ", sep="")
    cat(cli::ansi_strtrim(jobs$name[i], 35))
    cat("\n")
  }
  cat(paste(rep("-", 140), collapse = ""), "\n")

  leo.basic::leo_time_elapsed(t0)
  return(invisible(jobs))
}
