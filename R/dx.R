#' Extract UKB data using DNAnexus Table Exporter
#'
#' This function takes a file containing UKB field IDs and runs the Table Exporter
#' app on DNAnexus RAP to extract the data. The output will be saved to `data_files/` on the RAP project.
#'
#' @param file Path to a text file containing UKB field IDs (one per line).
#' @param output_prefix Output file prefix (default: derived from input file name).
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
#' @param dataset The dataset ID/name to extract from (e.g., \code{"app12345_20240101.dataset"} or \code{"project-Gk2...:record-Fp3..."}). 
#'   Defaults to the latest dataset in the project if NULL.
#' @param instance_type The DNAnexus instance type for the job. Default: \code{"mem1_ssd1_v2_x4"}.
#'   For large extracts, \code{"mem1_hdd1_v2_x8"} is recommended.
#' @return A list containing the job ID and expected output path (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Basic extraction (extracts covariates defined in fields.txt)
#' leo_dx_extract("tmp/fields.txt", output_prefix = "ukb_data")
#'
#' # 2. Extract with human-readable labels and descriptive headers
#' leo_dx_extract("tmp/fields.txt", 
#'                output_prefix = "ukb_data_readable",
#'                coding_option = "REPLACE",
#'                header_style = "FIELD-TITLE")
#'                
#' # 3. Explicitly specifying the dataset (Recommended for reproducibility)
#' leo_dx_extract("tmp/fields.txt", 
#'                dataset = "project-Gk2PzX0Jj1X4Y5Z6:record-Fp37890Qj9k2X1Y4")
#' }
leo_dx_extract <- function(file,
                           output_prefix = NULL,
                           expand = TRUE,
                           entity = "participant",
                           coding_option = "RAW",
                           output_format = "CSV",
                           header_style = "FIELD-NAME",
                           instance_type = "mem1_ssd1_v2_x4",
                           dataset = NULL) {
  # Check if file exists
  if (!file.exists(file)) {
    leo.basic::leo_log("File not found: {file}", level = "danger")
    return(invisible(NULL))
  }
  cli::cat_rule("Extracting Data using DNAnexus Table Exporter", col = "blue")
  
  # Display dx environment info
  tryCatch({
    dx_env <- system("dx env", intern = TRUE, ignore.stderr = TRUE)
    user_line <- grep("Current user", dx_env, value = TRUE)
    project_line <- grep("Current workspace name", dx_env, value = TRUE)
    if (length(user_line) > 0) leo.basic::leo_log(trimws(user_line), level = "info")
    if (length(project_line) > 0) leo.basic::leo_log(trimws(project_line), level = "info")
  }, error = function(e) {
    leo.basic::leo_log("Could not retrieve dx environment info", level = "warning")
  })

  # Auto-detect latest dataset if not specified
  # Use `dx find data` instead of `dx ls` for reliability:
  # - dx ls is non-recursive and shell glob can be consumed by zsh
  # - dx find data searches the entire project
  if (is.null(dataset)) {
    datasets <- tryCatch({
      # --brief returns just the ID, --name filters by pattern
      system("dx find data --name '*.dataset' --brief", intern = TRUE, ignore.stderr = TRUE)
    }, error = function(e) character(0))

    if (length(datasets) > 0) {
      # Take the first (or latest by timestamp if names include dates)
      # dx find data returns newest first by default
      dataset <- datasets[1]
      leo.basic::leo_log("Auto-detected dataset: {dataset}", level = "info")
    } else {
      leo.basic::leo_log("No dataset found in project. Please specify the dataset parameter.", level = "danger")
      return(invisible(NULL))
    }
  }

  # Read and validate field IDs
  fields <- readLines(file, warn = FALSE)
  fields <- trimws(fields)
  fields <- fields[fields != ""]

  if (length(fields) == 0) {
    leo.basic::leo_log("No valid field IDs found in {file}", level = "danger")
    return(invisible(NULL))
  }

  # Ensure eid is always the first field
  fields <- fields[!tolower(fields) %in% c("eid", "p_eid")]
  fields <- unique(fields)  # Remove duplicates
  fields <- c("eid", fields)
  leo.basic::leo_log("Loaded {length(fields)} field IDs (eid + {length(fields) - 1} fields)", level = "info")

  # Expand fields
  # New Logic: Official Dataset Dictionary Validation - we fetch the actual dataset dictionary and filter.
  if (expand) {
    leo.basic::leo_log("Validating fields against dataset dictionary ({dataset})...", level = "info")
    
    # 1. Get Dictionary (Official)
    ukb_dd <- leo_dx_get_dataset_dictionary(dataset)
    
    if (!is.null(ukb_dd)) {
        # 2. Find valid columns (Official Regex), filtered by entity
        # Exclude 'eid' from lookup as it's added manually later
        valid_cols <- leo_dx_find_columns(fields[fields != "eid"], ukb_dd, entity = entity)
        
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

  # Generate output prefix from file name if not provided
  if (is.null(output_prefix)) {
    output_prefix <- tools::file_path_sans_ext(basename(file))
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
  dx_file_path <- paste0("fields_files/", dx_file_name)
  leo.basic::leo_log("Uploading field names file to: {dx_file_path}", level = "info")

  # Delete existing file if present (to avoid duplicates)
  system(paste0("dx rm -f ", dx_file_path), ignore.stdout = TRUE, ignore.stderr = TRUE)

  upload_cmd <- paste0("dx upload \"", fields_file, "\" --path ", dx_file_path)
  upload_exit <- tryCatch({
    system(upload_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }, error = function(e) {
    leo.basic::leo_log(paste0("Upload error: ", e$message), level = "danger")
    return(1L)
  })

  if (upload_exit != 0L) {
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
    "--destination", "data_files/",
    "--priority", "high",
    "--instance-type", instance_type,
    "-y"  # auto-confirm
 )

  leo.basic::leo_log("Running Table Exporter...", level = "info")
  leo.basic::leo_log("dx {paste(args, collapse = ' ')}", level = "info")

  res <- system2("dx", args = args, stdout = TRUE, stderr = TRUE)
  exit_code <- attr(res, "status")
  if (is.null(exit_code)) exit_code <- 0L

  if (exit_code != 0L) {
    leo.basic::leo_log("Table Exporter failed with exit code {exit_code}", level = "danger")
    leo.basic::leo_log("Output: {paste(res, collapse = '\\n')}", level = "info")
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
  leo.basic::leo_log("Output will be saved to: {output_path}", level = "info")

  return(invisible(list(job_id = job_id, output_path = output_path)))
}

#' Check DNAnexus Job Status
#'
#' Checks the status of a DNAnexus job (e.g. Table Exporter) using \`dx describe\`.
#'
#' @param job_id A string containing the job ID (e.g., "job-Gk2...") or a list 
#'   returned by \code{\link{leo_dx_extract}} (which contains a \`\$job_id\` element).
#' @return The state of the job as a string (e.g., "idle", "runnable", "running", 
#'   "done", "failed", "terminated"). Returns NA if the status could not be determined.
#' @export
#' @examples
#' \dontrun{
#' # Check status using a job ID
#' leo_dx_status("job-Gk2PzX0Jj1X4Y5Z6")
#'
#' # Check status from extraction result
#' res <- leo_dx_extract("fields.txt")
#' leo_dx_status(res)
#' }
leo_dx_status <- function(job_id) {
  # Handle list input (from leo_dx_extract)
  if (is.list(job_id) && !is.null(job_id$job_id)) {
    job_id <- job_id$job_id
  }
  
  if (is.null(job_id) || is.na(job_id) || !is.character(job_id) || !grepl("^job-", job_id)) {
    leo.basic::leo_log("Invalid Job ID provided.", level = "danger")
    return(NA_character_)
  }
  
  # Run dx describe
  # dx describe job-xxxx --json
  # We look for the "state" field.
  
  cmd <- glue::glue("dx describe \"{job_id}\" --json")
  res <- tryCatch({
    system(cmd, intern = TRUE, ignore.stderr = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(res) || length(res) == 0) {
    leo.basic::leo_log("Could not check status for {job_id}. (Is dx toolkit authenticated?)", level = "warning")
    return(NA_character_)
  }
  
  json_str <- paste(res, collapse = " ")
  
  # Simple regex for "state": "done"
  matches <- regmatches(json_str, regexec('"state"\\s*:\\s*"([^"]+)"', json_str))
  
  state <- NA_character_
  if (length(matches[[1]]) > 1) {
    state <- matches[[1]][2]
  }
  
  # Pretty print status
  if (!is.na(state)) {
    color_map <- list(
      "done" = "success",
      "failed" = "danger",
      "running" = "info",
      "idle" = "warning",
      "runnable" = "warning",
      "terminated" = "danger"
    )
    lvl <- if (!is.null(color_map[[state]])) color_map[[state]] else "info"
    
    leo.basic::leo_log("Job {job_id} is: {state}", level = lvl)
  }
  
  return(state)
}