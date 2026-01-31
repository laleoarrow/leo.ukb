#' Extract UKB data using DNAnexus Table Exporter
#'
#' This function takes a file containing UKB field IDs and runs the Table Exporter
#' app on DNAnexus RAP to extract the data. The output will be saved to `data_files/`.
#'
#' @param file Path to a text file containing UKB field IDs (one per line).
#' @param output_prefix Output file prefix (default: derived from input file name).
#' @param entity Entity type. Options:
#'   \itemize{
#'     \item \code{"participant"} (default): All UK Biobank participants. This is the main phenotype/covariate data.
#'     \item Other entity types depend on the dataset structure (e.g., "sample" for genetic data).
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
#'     \item \code{"UKB-FIELD-ID"}: Use raw UKB field IDs like \code{21003-0.0}.
#'   }
#' @param dataset The dataset name to extract from. Defaults to the latest dataset in the project.
#'   You can find available datasets with \code{dx ls *.dataset} in the terminal.
#'
#' @param instance_type The DNAnexus instance type for the job (default: \code{"mem1_ssd1_v2_x4"}).
#' @return A list containing the job ID and expected output path (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract covariates with default settings
#' leo_dx_extract("tmp/cov20260131.txt", output_prefix = "cov20260131")
#'
#' # Use human-readable labels and descriptive column names
#' leo_dx_extract("tmp/cov20260131.txt",
#'                output_prefix = "cov20260131",
#'                coding_option = "REPLACE",
#'                header_style = "FIELD-TITLE")
#' }
leo_dx_extract <- function(file,
                           output_prefix = NULL,
                           entity = "participant",
                           coding_option = "RAW",
                           output_format = "CSV",
                           header_style = "FIELD-NAME",
                           instance_type = "mem1_ssd1_v2_x4",
                           dataset = NULL) {

  # Display dx environment info
  cli::cat_rule("DNAnexus Environment", col = "blue")
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
  if (is.null(dataset)) {
    datasets <- tryCatch({
      system("dx ls *.dataset", intern = TRUE, ignore.stderr = TRUE)
    }, error = function(e) character(0))

    if (length(datasets) > 0) {
      # Sort by name (which includes timestamp) and take the last one
      dataset <- sort(datasets, decreasing = TRUE)[1]
      leo.basic::leo_log("Auto-detected dataset: {dataset}", level = "info")
    } else {
      leo.basic::leo_log("No dataset found in project. Please specify the dataset parameter.", level = "danger")
      return(invisible(NULL))
    }
  }

  if (!file.exists(file)) {
    leo.basic::leo_log("File not found: {file}", level = "danger")
    return(invisible(NULL))
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

  # Generate output prefix from file name if not provided
  if (is.null(output_prefix)) {
    output_prefix <- tools::file_path_sans_ext(basename(file))
  }

  # Create temp file with eid as first field
  tmp_file <- tempfile(fileext = ".txt")
  writeLines(fields, tmp_file)
  on.exit(unlink(tmp_file), add = TRUE)

  # Upload the field names file to fields_files/ folder
  dx_file_name <- paste0(output_prefix, "_fields.txt")
  dx_file_path <- paste0("fields_files/", dx_file_name)
  leo.basic::leo_log("Uploading field names file to: {dx_file_path}", level = "info")

  # Delete existing file if present (to avoid duplicates)
  system(paste0("dx rm -f ", dx_file_path), ignore.stdout = TRUE, ignore.stderr = TRUE)

  upload_cmd <- paste0("dx upload \"", tmp_file, "\" --path ", dx_file_path)
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
