# Analysis until tools ----
#' Basic imputation function for missing values
#'
#' This function provides multiple methods for imputing missing values in a dataset,
#' including mean, median, random forest, KNN, and multiple imputation.
#'
#' @param data A data frame containing the data with missing values
#' @param method Imputation method: "mean", "median", "rf" (random forest),
#'               "knn" (K-nearest neighbors), or "mice" (multiple imputation)
#' @param ... Additional arguments passed to the specific imputation function
#'
#' @return A data frame with missing values imputed
#' @export
#' @importFrom dplyr across mutate where
#' @importFrom stats median
#' @importFrom VIM kNN
#' @importFrom missRanger missRanger
#' @importFrom mice mice complete
#'
#' @examples
#' # Create sample data with missing values
#' set.seed(123)
#' sample_data <- data.frame(
#'   age = c(25, 30, NA, 40, 45),
#'   score = c(80, NA, 90, 85, NA),
#'   group = factor(c("A", "B", "A", NA, "B"))
#' )
#'
#' # Mean imputation
#' leo_impute_na(sample_data, method = "mean")
#'
#' # Random forest imputation
#' leo_impute_na(sample_data, method = "rf")
leo_impute_na <- function(data, method = "mean", ...) {
  if (!method %in% c("mean", "median", "rf", "knn", "mice")) {
    stop("Invalid method. Choose 'mean', 'median', 'rf', 'knn', or 'mice'")
  }

  # Calculate missing values per column
  na_before <- colSums(is.na(data))
  na_columns <- na_before[na_before > 0]

  if (length(na_columns) > 0) {
    leo_log("Missing values per column (count/percentage):")
    for (col_name in names(na_columns)) {
      na_count <- na_columns[col_name]
      na_pct <- round(na_count / nrow(data) * 100, 1)
      leo_log("  {col_name}: {na_count}/{nrow(data)} ({na_pct}%)")
    }
  } else {
    leo_log("No missing values found")
    return(data)
  }

  # Apply selected imputation method
  if (method == "mean") {
    imputed_data <- data %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  } else if (method == "median") {
    imputed_data <- data %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ifelse(is.na(.), stats::median(., na.rm = TRUE), .)))
  } else if (method == "rf") {
    if (!requireNamespace("missRanger", quietly = TRUE)) {
      stop("Package 'missRanger' required for random forest imputation")
    }
    imputed_data <- missRanger::missRanger(data, ...)
  } else if (method == "knn") {
    if (!requireNamespace("VIM", quietly = TRUE)) {
      stop("Package 'VIM' required for KNN imputation")
    }
    imputed_data <- VIM::kNN(data, ...)
  } else if (method == "mice") {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' required for multiple imputation")
    }
    imputed <- mice::mice(data, ...)
    imputed_data <- mice::complete(imputed)
  }

  # Check remaining missing values
  na_after <- colSums(is.na(imputed_data))
  remaining_na <- na_after[na_after > 0]

  if (length(remaining_na) > 0) {
    leo_log("Remaining missing values after imputation:", level = "warning")
    for (col_name in names(remaining_na)) {
      leo_log("  {col_name}: {remaining_na[col_name]}", level = "warning")
    }
  } else {
    leo_log("Imputation completed using {method} method")
  }

  return(imputed_data)
}
#' Bin a continuous variable by quantiles
#'
#' A lightweight wrapper of `cut()` using quantile-based breaks (equal-frequency).
#'
#' @param x Numeric/character vector.
#' @param bins Integer; number of quantile bins (default 4).
#' @param labels Character vector; default uses "Q1..Qk".
#' @param right Logical; passed to `cut()`.
#' @param include_lowest Logical; passed to `cut()`.
#' @param ... Reserved for future use.
#'
#' @return Character vector of quantile bins; NA preserved.
#' @export
#' @importFrom cli cli_alert_warning
#' @examples
#' x <- c(1:10, NA_integer_)
#' q4 <- leo_quantile(x, bins = 4); q4
#' table(q4, useNA = "ifany")
#'
#' x_tie <- c(rep(1L, 6), rep(2L, 6))
#' q4_tie <- leo_quantile(x_tie, bins = 4); q4_tie
#' table(q4_tie, useNA = "ifany")
leo_quantile <- function(x, bins = 4, labels = NULL, right = TRUE, include_lowest = TRUE, ...) {
  x_num <- suppressWarnings(as.numeric(x))
  probs <- seq(0, 1, length.out = bins + 1)
  breaks_use <- stats::quantile(x_num, probs = probs, na.rm = TRUE, type = 7, names = FALSE)

  breaks_use <- breaks_use[is.finite(breaks_use)]
  if (length(breaks_use) < 2) {
    cli::cli_alert_warning("leo_quantile(): insufficient finite break points; returning all NA.")
    return(rep(NA_character_, length(x_num)))
  }

  if (anyDuplicated(breaks_use)) {
    cli::cli_alert_warning("leo_quantile(): duplicated break points detected (ties). Collapsing duplicates.")
    breaks_use <- unique(breaks_use)
    if (length(breaks_use) < 2) return(rep(NA_character_, length(x_num)))
  }

  n_bins <- length(breaks_use) - 1
  if (is.null(labels)) labels <- paste0("Q", seq_len(n_bins))

  out <- cut(x_num, breaks = breaks_use, labels = labels, right = right, include.lowest = include_lowest)
  as.character(out)
}

#' @title Calculate P-value for Heterogeneity from Subgroup Summary Statistics
#' @description This function performs a test for heterogeneity (interaction)
#'              based on summary statistics (Hazard Ratios and P-values) from
#'              a subgroup analysis of a Cox model. It uses Cochran's Q test.
#' @param hrs A numeric vector of Hazard Ratios (HR) for each subgroup.
#' @param p_values A numeric vector of the corresponding P-values for each HR.
#' @param subgroup_names An optional character vector of subgroup names for display.
#'
#' @return A list containing the Q statistic, degrees of freedom (df),
#'         the P-value for heterogeneity, and the I-squared statistic.
#'         It also prints a summary to the console.
#' @examples
#' # Example: Subgroup analysis by smoking status (Smoker vs. Non-smoker)
#' hrs_example <- c(3.576652, 2.657723366)
#' p_values_example <- c(0.012254053, 0.002129672)
#' names_example <- c("Age>65", "Age<=65")
#' leo_heterogeneity_p(
#'   hrs = hrs_example,
#'   p_values = p_values_example,
#'   subgroup_names = names_example
#' )
#' @export
leo_heterogeneity_p <- function(hrs, p_values, subgroup_names = NULL) {

  # --- 1. Input Validation ---
  if (length(hrs) != length(p_values)) stop("Error: The length of 'hrs' and 'p_values' must be the same.")
  if (length(hrs) < 2) stop("Error: At least two subgroups are required to test for heterogeneity.")
  if (any(hrs <= 0) || any(p_values <= 0) || any(p_values > 1)) stop("Error: Invalid input. HRs must be > 0 and P-values must be between 0 and 1.")

  # --- 2. Calculate log(HR) and Standard Errors (SE) ---
  # The test statistic for a single HR is z = log(HR) / SE(log(HR))
  # The p-value is derived from z. We can reverse this to find SE.
  # p = 2 * pnorm(-abs(z)) => z = abs(qnorm(p / 2))

  log_hrs <- log(hrs)
  z_scores <- abs(qnorm(p_values / 2))

  # Handle cases where HR is exactly 1 (log_hr is 0)
  # In this case, SE can't be calculated by division.
  # If HR=1, p should be 1, and z=0. SE would be infinite, weight=0.
  # We add a small epsilon to avoid division by zero if p is close to 1.
  z_scores[z_scores == 0] <- 1e-10

  ses <- abs(log_hrs / z_scores)

  # If a log_hr is 0 (HR=1), its SE should result in a weight of 0
  # unless p-value is also 1. If p-value is not 1, z is not 0, so SE is 0.
  # A SE of 0 gives infinite weight, which is problematic.
  # A robust way is to assign a very small weight (or handle as a special case),
  # but in practice, if log(HR) is 0, its contribution to the Q-stat is 0
  # if weighted properly. Let's use weights directly.

  # --- 3. Calculate Weights (Inverse Variance) ---
  weights <- 1 / (ses^2)

  # Handle cases where SE was 0 (due to log_hr being 0 but p not 1)
  weights[is.infinite(weights)] <- 1e10 # Assign a very large weight

  # --- 4. Perform Cochran's Q Test ---
  # Calculate the weighted average of the log(HR)s
  pooled_log_hr <- sum(weights * log_hrs) / sum(weights)

  # Calculate the Q statistic
  # Q = sum(w_i * (log_hr_i - pooled_log_hr)^2)
  Q <- sum(weights * ((log_hrs - pooled_log_hr)^2))

  # Degrees of freedom = number of subgroups - 1
  df <- length(hrs) - 1

  # --- 5. Calculate Heterogeneity P-value and I^2 ---
  # The Q statistic follows a Chi-squared distribution with df degrees of freedom
  p_heterogeneity <- pchisq(Q, df, lower.tail = FALSE)

  # Calculate I^2 statistic
  # I^2 = ((Q - df) / Q) * 100%
  # It represents the percentage of variation across subgroups due to heterogeneity
  i_squared <- max(0, (Q - df) / Q) * 100

  # --- 6. Format and Print Output ---
  if (is.null(subgroup_names)) {
    subgroup_names <- paste("Subgroup", 1:length(hrs))
  }

  cat("--- Test for Heterogeneity in Subgroup Analysis ---\n\n")
  cat("Input Data:\n")
  print(data.frame(
    Subgroup = subgroup_names,
    HR = hrs,
    P_value = p_values,
    log_HR = round(log_hrs, 4),
    SE_log_HR = round(ses, 4)
  ))

  cat("\nResults of Cochran's Q Test:\n")
  cat("Q-statistic:", round(Q, 4), "\n")
  cat("Degrees of Freedom (df):", df, "\n")
  cat("P-value for Heterogeneity:", format.pval(p_heterogeneity, digits = 3, eps = 0.001), "\n\n")

  cat("Heterogeneity Metrics:\n")
  cat("I^2 statistic:", paste0(round(i_squared, 2), "%"), "\n")
  cat("Interpretation of I^2: ~0-40% (might not be important); 30-60% (may represent moderate heterogeneity); 50-90% (may represent substantial heterogeneity); 75-100% (considerable heterogeneity).\n")
  cat("\n---------------------------------------------------\n")

  # Return results as a list
  invisible(list(
    Q_statistic = Q,
    df = df,
    p_value_heterogeneity = p_heterogeneity,
    I2_statistic_percent = i_squared
  ))
}

# Helper for regression ----

#' Check whether a variable should be treated as continuous or categorical
#'
#' @param x Input vector.
#' @param var_name Character scalar used in messages.
#' @param var_type One of `"auto"`, `"continuous"`, or `"categorical"`.
#' @param verbose Logical; whether to print warning messages.
#' @return Character scalar: `"continuous"` or `"categorical"`.
#' @examples
#' leo.ukb:::.check_var_type(c(0.2, 1.1, -0.3), var_name = "prs")
#' leo.ukb:::.check_var_type(factor(c("Low", "High")), var_name = "group")
#' leo.ukb:::.check_var_type(c(0, 1, 2, 0, 1), var_name = "code", var_type = "categorical")
#' @noRd
.check_var_type <- function(x, var_name = "variable", var_type = c("auto", "continuous", "categorical"), verbose = TRUE) {
  var_type <- match.arg(var_type)
  x_non_missing <- stats::na.omit(x)
  n_unique <- length(unique(x_non_missing))
  if (var_type == "categorical") {
    if (n_unique < 2) stop(var_name, " must contain at least 2 non-missing values when var_type = 'categorical'.", call. = FALSE)
    if (is.factor(x) || is.character(x) || is.logical(x)) return("categorical")
    if (!is.numeric(x)) stop(var_name, " has an unsupported type for categorical handling.", call. = FALSE)
    unique_values <- sort(unique(x_non_missing))
    integer_like <- all(abs(unique_values - round(unique_values)) < 1e-8)
    if (!integer_like || length(unique_values) > 10) stop(var_name, " looks continuous and cannot be forced to categorical. Convert it to factor() first if this is intentional.", call. = FALSE)
    return("categorical")
  }
  if (var_type == "continuous") {
    if (is.character(x) || is.factor(x) || is.logical(x)) stop(var_name, " requires a numeric or integer variable when var_type = 'continuous'.", call. = FALSE)
    return("continuous")
  }
  if (is.factor(x) || is.character(x) || is.logical(x)) return("categorical")
  if (!is.numeric(x)) stop(var_name, " has an unsupported type for auto detection.", call. = FALSE)
  unique_values <- sort(unique(x_non_missing))
  integer_like <- length(unique_values) > 1 && all(abs(unique_values - round(unique_values)) < 1e-8)
  if (integer_like && length(unique_values) <= 5) {
    leo.basic::leo_log("Exposure '{var_name}' looks like a coded categorical variable ({paste(unique_values, collapse = ', ')}). Convert it to factor() or set x_exp_type = 'categorical' if that is your intent.", level = "warning", verbose = verbose)
  }
  "continuous"
}

#' Helper function to format p-values according to common publication standards
#'
#' @param p_value Numeric p-value
#' @param fmt Format style: `"threshold"` (default, uses `<0.001`) or
#'   `"scientific"` (uses scientific notation for p < 0.001, e.g. `2.300e-04`).
#' @return Formatted p-value string
#' @examples
#' leo.ukb:::.format_p_value(0.2)
#' leo.ukb:::.format_p_value(0.0314)
#' leo.ukb:::.format_p_value(0.0002)
#' leo.ukb:::.format_p_value(0.0002, fmt = "scientific")
#' @noRd
.format_p_value <- function(p_value, fmt = "threshold") {
  if (is.na(p_value)) return("NA")
  if (fmt == "scientific" && p_value < 0.001) return(sprintf("%.3e", p_value))
  if (p_value < 0.001) return("<0.001")
  sprintf("%.3f", round(p_value, 3))
}

#' Normalize model input into a named covariate list
#'
#' @param x_cov `NULL`, a character vector, or a list of character vectors.
#' @param df_colnames Optional vector of available column names for validation.
#' @return A named list of covariate vectors for downstream regression.
#' @examples
#' leo.ukb:::.normalize_model_list(NULL)
#' leo.ukb:::.normalize_model_list(c("age", "sex"), df_colnames = c("age", "sex", "bmi"))
#' leo.ukb:::.normalize_model_list(list(NULL, c("age"), c("age", "sex")), df_colnames = c("age", "sex"))
#' @noRd
.normalize_model_list <- function(x_cov, df_colnames = NULL) {
  if (is.null(x_cov)) {
    model_list <- list(NULL)
  } else if (is.character(x_cov)) {
    model_list <- list(x_cov)
  } else if (is.list(x_cov)) {
    model_list <- x_cov
  } else {
    stop("x_cov must be NULL, a character vector, or a list of character vectors.", call. = FALSE)
  }

  model_list <- lapply(model_list, function(x) {
    if (is.null(x)) return(NULL)
    if (!is.character(x)) stop("Each model in x_cov must be NULL or a character vector.", call. = FALSE)
    unique(stats::na.omit(x))
  })
  default_names <- paste0("model_", seq_along(model_list))
  model_names <- names(model_list)
  if (is.null(model_names)) model_names <- rep("", length(model_list))
  model_names[is.na(model_names) | model_names == ""] <- default_names[is.na(model_names) | model_names == ""]
  names(model_list) <- make.unique(model_names)

  if (!is.null(df_colnames)) {
    missing_vars <- lapply(model_list, function(x) setdiff(if (is.null(x)) character(0) else x, df_colnames))
    bad_models <- which(vapply(missing_vars, length, integer(1)) > 0)
    if (length(bad_models) > 0) {
      msg <- paste(vapply(bad_models, function(i) {
        paste0(names(model_list)[i], " (", paste(missing_vars[[i]], collapse = ", "), ")")
      }, character(1)), collapse = "; ")
      stop("Missing covariates in df: ", msg, call. = FALSE)
    }
  }

  model_list
}

#' Prepare a shared regression dataset with one exposure and one or more outcomes
#'
#' @param df Analysis data frame.
#' @param y_cols Character vector of outcome column names to keep.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov Covariate specification passed to `.normalize_model_list()`.
#' @return A list containing complete-case data, model definitions, and metadata.
#' @examples
#' df <- data.frame(
#'   y = c(1, 0, 1),
#'   prs = c(0.2, -0.1, 0.8),
#'   age = c(50, 60, NA)
#' )
#' leo.ukb:::.prepare_regression_data(df = df, y_cols = "y", x_exp = "prs", x_cov = "age")
#' @noRd
.prepare_regression_data <- function(df, y_cols, x_exp, x_cov = NULL) {
  if (!is.data.frame(df)) stop("df must be a data.frame.", call. = FALSE)
  if (!is.character(y_cols) || length(y_cols) < 1) stop("y_cols must be a character vector of column names.", call. = FALSE)
  if (!all(y_cols %in% names(df))) stop("All y_cols must exist in df.", call. = FALSE)
  if (!is.character(x_exp) || length(x_exp) != 1) stop("x_exp must be a single column name.", call. = FALSE)
  if (!x_exp %in% names(df)) stop("x_exp must exist in df.", call. = FALSE)

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  covariates <- unique(unlist(model_list, use.names = FALSE))
  base_df <- df[, y_cols, drop = FALSE]
  base_df$exposure <- df[[x_exp]]
  if (length(covariates) > 0) base_df[covariates] <- df[covariates]
  complete_df <- stats::na.omit(base_df)
  missing_vars <- names(colSums(is.na(base_df))[colSums(is.na(base_df)) > 0])

  list(
    data = complete_df,
    models = model_list,
    y_cols = y_cols,
    exposure_name = x_exp,
    covariates = covariates,
    n_total = nrow(base_df),
    n_after_complete_case = nrow(complete_df),
    n_removed_complete_case = nrow(base_df) - nrow(complete_df),
    missing_vars = missing_vars
  )
}

#' Prepare a single-exposure Cox regression dataset
#'
#' @param df Analysis data frame.
#' @param y_out Character vector of length 2 giving the event and follow-up time
#'   column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov Covariate specification passed to `.normalize_model_list()`;
#'   each entry should be a covariate column name. All models are fitted on the
#'   same complete-case cohort defined by `y_out`, `x_exp`, and all covariates
#'   that appear in `x_cov`.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#'   Default is `0`, which excludes pre-baseline or baseline events when defining
#'   an incident outcome.
#' @param event_value Value in the event column that indicates incident events.
#' @param verbose Logical; whether to print progress messages.
#' @return A list containing filtered data, model definitions, and analysis metadata.
#' @examples
#' df <- data.frame(
#'   outcome = c(1, 0, 1),
#'   outcome_time = c(2, 0, 5),
#'   prs = c(0.2, -0.1, 0.8),
#'   age = c(50, 60, 55)
#' )
#' leo.ukb:::.prepare_cox_regression_data(
#'   df = df,
#'   y_out = c("outcome", "outcome_time"),
#'   x_exp = "prs",
#'   x_cov = "age",
#'   min_followup_time = 0,
#'   verbose = FALSE
#' )
#' @noRd
.prepare_cox_regression_data <- function(df, y_out, x_exp, x_cov = NULL, min_followup_time = 0, event_value = 1, verbose = TRUE) {
  if (!is.data.frame(df)) stop("df must be a data.frame.", call. = FALSE)
  if (!is.character(y_out) || length(y_out) != 2) stop("y_out must be a character vector of length 2: c(event, time).", call. = FALSE)
  if (!all(y_out %in% names(df))) stop("All y_out columns must exist in df.", call. = FALSE)
  if (!is.character(x_exp) || length(x_exp) != 1) stop("x_exp must be a single column name.", call. = FALSE)
  if (!x_exp %in% names(df)) stop("x_exp must exist in df.", call. = FALSE)
  if (!is.numeric(min_followup_time) || length(min_followup_time) != 1) stop("min_followup_time must be a numeric scalar.", call. = FALSE)
  if (length(event_value) != 1 || is.na(event_value)) stop("event_value must be a single non-missing value.", call. = FALSE)

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  coerce_followup_time <- function(x, column_name) {
    x_cmp <- if (is.factor(x)) as.character(x) else x
    if (is.logical(x_cmp)) stop("Follow-up time column '", column_name, "' must be numeric or numeric-like, not logical.", call. = FALSE)
    time_num <- suppressWarnings(as.numeric(x_cmp))
    if (any(!is.na(x_cmp) & is.na(time_num))) {
      stop("Follow-up time column '", column_name, "' must be numeric or numeric-like.", call. = FALSE)
    }
    time_num
  }
  coerce_event_indicator <- function(x, event_value, column_name) {
    x_cmp <- if (is.factor(x)) as.character(x) else x
    event_value_chr <- as.character(event_value)
    if (is.logical(x_cmp)) {
      if (event_value_chr %in% c("1", "TRUE", "T", "true")) return(ifelse(is.na(x_cmp), NA_integer_, as.integer(x_cmp)))
      if (event_value_chr %in% c("0", "FALSE", "F", "false")) return(ifelse(is.na(x_cmp), NA_integer_, as.integer(!x_cmp)))
      stop("For logical event column '", column_name, "', event_value must be TRUE/FALSE (or 1/0).", call. = FALSE)
    }
    observed_chr <- unique(as.character(stats::na.omit(x_cmp)))
    if (!(event_value_chr %in% observed_chr)) {
      stop(
        "event_value='", event_value_chr, "' does not match the observed values in '", column_name, "': ",
        paste(observed_chr, collapse = ", "), ". Please set event_value explicitly.",
        call. = FALSE
      )
    }
    ifelse(is.na(x_cmp), NA_integer_, as.integer(as.character(x_cmp) == event_value_chr))
  }

  event_raw <- df[[y_out[1]]]
  time_raw <- coerce_followup_time(df[[y_out[2]]], y_out[2])
  event <- coerce_event_indicator(event_raw, event_value, y_out[1])

  base_df <- data.frame(event = event, time = time_raw, exposure = df[[x_exp]], stringsAsFactors = FALSE)
  covariates <- unique(unlist(model_list, use.names = FALSE))
  if (length(covariates) > 0) {
    base_df[covariates] <- df[covariates]
  }

  n_total <- nrow(base_df)
  filtered_df <- base_df[!is.na(base_df$time) & base_df$time > min_followup_time, , drop = FALSE]
  n_after_followup <- nrow(filtered_df)
  prep <- .prepare_regression_data(df = filtered_df, y_cols = c("event", "time"), x_exp = "exposure", x_cov = model_list)

  leo.basic::leo_log("Cox data prepared: total = {n_total}, after follow-up filter = {n_after_followup}, after complete-case filter = {prep$n_after_complete_case}.", verbose = verbose)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("Complete-case filtering removed {prep$n_removed_complete_case} row(s) with missing values across exposure/outcome/model covariates{if (length(prep$missing_vars) > 0) paste0(': ', paste(prep$missing_vars, collapse = ', ')) else ''}.", level = "warning", verbose = verbose)
  }

  list(
    data = prep$data,
    models = prep$models,
    outcome_name = y_out[1],
    outcome_time = y_out[2],
    exposure_name = x_exp,
    n_total = n_total,
    n_after_followup = n_after_followup,
    n_after_complete_case = prep$n_after_complete_case,
    n_removed_complete_case = prep$n_removed_complete_case,
    missing_vars = prep$missing_vars,
    min_followup_time = min_followup_time,
    event_value = event_value
  )
}

# Analyses ----

# Cox ----
#' Cox regression fitting and formatting helpers
#'
#' `leo_cox()` fits one or more Cox proportional hazards models for a
#' single exposure and incident outcome. `leo_cox_format()` converts
#' the returned object into a wide summary table, tidy result table, or
#' `gtsummary` output.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time
#'   column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list
#'   of covariate column-name vectors. All models are fitted on the same
#'   complete-case cohort defined by `y_out`, `x_exp`, and all covariates that
#'   appear in `x_cov`. If `x_cov` is a named list, those model names are
#'   preserved in the output columns.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#'   Default is `0`, which excludes pre-baseline or baseline events when defining
#'   an incident outcome.
#' @param event_value Value in the event column that indicates incident events.
#' @param verbose Logical; print progress messages.
#' @param x_exp_type Exposure type handling for `x_exp`. Use `"auto"` to infer from the input
#'   type, `"continuous"` to force a numeric Cox term, or `"categorical"` to
#'   force factor coding. In `"auto"` mode, small integer-coded exposures will
#'   trigger a warning because they may represent categorical groups.
#' @param simplify Output mode: `"wide"` (default) returns the wide result
#'   table directly; `"tidy"` returns the tidy result table; `FALSE` returns
#'   the full `leo_cox` object.
#' @param p_fmt P-value display format: `"threshold"` (default) uses `<0.001`;
#'   `"scientific"` uses scientific notation for p < 0.001 (e.g. `2.300e-04`).
#' @return When `simplify = FALSE`, a `leo_cox` object containing the default
#'   wide table in `$result`, the underlying tidy rows in `$result_tidy`, model
#'   metadata, and fitted `coxph` objects. When `simplify = "wide"` or
#'   `"tidy"`, the corresponding formatted table directly.
#' @export
#' @examples
#' lung_df <- stats::na.omit(
#'   dplyr::transmute(
#'     survival::lung,
#'     outcome = as.integer(status == 2),
#'     outcome_censor = time / 365.25,
#'     age = age,
#'     sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'     ecog_group = factor(ph.ecog, levels = 0:3, labels = c("ECOG0", "ECOG1", "ECOG2", "ECOG3"))
#'   )
#' ); head(lung_df)
#'
#' model_cont <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("sex"),
#'   "Model B" = c("sex", "ecog_group")
#' )
#'
#' res_age <- leo_cox(
#'   df = lung_df, y_out = c("outcome", "outcome_censor"),
#'   x_exp = "age", x_cov = model_cont, simplify = FALSE, verbose = FALSE
#' )
#' res_age$result
#' leo_cox_format(res_age, style = "tidy")
#'
#' model_cat <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("sex")
#' )
#'
#' res_ecog <- leo_cox(
#'   df = lung_df, y_out = c("outcome", "outcome_censor"),
#'   x_exp = "ecog_group", x_cov = model_cat,
#'   x_exp_type = "categorical", simplify = FALSE, verbose = FALSE
#' )
#' res_ecog$result
#'
#' if (requireNamespace("gtsummary", quietly = TRUE) && requireNamespace("broom.helpers", quietly = TRUE)) {
#'   leo_cox_format(res_age, style = "gtsummary")
#' }
leo_cox <- function(df, y_out, x_exp, x_cov = NULL, event_value = 1, min_followup_time = 0, x_exp_type = "auto", simplify = "wide", p_fmt = "threshold", verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Regression", col = "blue")
  if (!x_exp_type %in% c("auto", "continuous", "categorical")) stop("x_exp_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)

  # Normalize and validate inputs
  prep <- .prepare_cox_regression_data(
    df = df,
    y_out = y_out,
    x_exp = x_exp,
    x_cov = x_cov,
    min_followup_time = min_followup_time,
    event_value = event_value,
    verbose = verbose
  )
  if (prep$n_after_followup == 0) stop("No rows remain after follow-up filtering.", call. = FALSE)
  if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering.", call. = FALSE)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before Cox analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_cox(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }
  exposure_type <- .check_var_type(prep$data$exposure, var_name = prep$exposure_name, var_type = x_exp_type, verbose = verbose)

  # Loop over models
  model_results <- list()
  fit_results <- list()
  for (model_name in names(prep$models)) {
    covariates <- if (is.null(prep$models[[model_name]])) character(0) else prep$models[[model_name]]
    leo.basic::leo_log("Fitting {model_name} with {if (length(covariates) == 0) 'no covariates (crude model)' else paste(covariates, collapse = ', ')}.", verbose = verbose)

    # Prepare analysis data
    model_df <- prep$data[, c("event", "time", "exposure", covariates), drop = FALSE]

    # Fit Cox model
    if (exposure_type == "categorical") {
      model_df$exposure <- droplevels(factor(model_df$exposure))
      if (nlevels(model_df$exposure) < 2) stop("x_exp must contain at least 2 exposure levels after filtering.", call. = FALSE)
    } else {
      if (is.character(model_df$exposure) || is.factor(model_df$exposure) || is.logical(model_df$exposure)) {
        stop("x_exp_type = 'continuous' requires a numeric or integer exposure.", call. = FALSE)
      }
      model_df$exposure <- suppressWarnings(as.numeric(model_df$exposure))
    }
    rhs <- "exposure"
    if (length(covariates) > 0) rhs <- c(rhs, paste0("`", covariates, "`"))
    formula <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs, collapse = " + ")))
    fit <- survival::coxph(formula = formula, data = model_df)
    tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    total_case <- sum(model_df$event == 1, na.rm = TRUE)
    total_control <- sum(model_df$event == 0, na.rm = TRUE)
    total_years <- sum(model_df$time, na.rm = TRUE)
    formula_txt <- Reduce(paste, deparse(formula))

    # Extract HR / CI / p-value and build tidy result rows
    if (!is.factor(model_df$exposure)) {
      exposure_row <- tidy_fit[tidy_fit$term == "exposure", , drop = FALSE]
      model_results[[model_name]] <- data.frame(
        model = model_name,
        row_id = 1L,
        exposure = prep$exposure_name,
        outcome = prep$outcome_name,
        level = NA_character_,
        case_n = total_case,
        control_n = total_control,
        person_year = total_years,
        hr = exposure_row$estimate[1],
        hr_ci_l = exposure_row$conf.low[1],
        hr_ci_u = exposure_row$conf.high[1],
        p_value = exposure_row$p.value[1],
        exposure_class = "Continuous",
        formula = formula_txt,
        check.names = FALSE
      )
      fit_results[[model_name]] <- fit
      next
    }

    levels_x <- levels(model_df$exposure)
    level_counts <- do.call(rbind, lapply(levels_x, function(level) {
      level_df <- model_df[model_df$exposure == level, , drop = FALSE]
      data.frame(
        level = level,
        Case_N = sum(level_df$event == 1, na.rm = TRUE),
        Control_N = sum(level_df$event == 0, na.rm = TRUE),
        person_year = sum(level_df$time, na.rm = TRUE)
      )
    }))
    result_df <- data.frame(
      model = model_name,
      row_id = seq_along(levels_x),
      exposure = c(paste0(levels_x[1], " (Ref)"), levels_x[-1]),
      outcome = prep$outcome_name,
      level = levels_x,
      case_n = level_counts$Case_N[match(levels_x, level_counts$level)],
      control_n = level_counts$Control_N[match(levels_x, level_counts$level)],
      person_year = level_counts$person_year[match(levels_x, level_counts$level)],
      hr = c(1, rep(NA_real_, length(levels_x) - 1)),
      hr_ci_l = c(1, rep(NA_real_, length(levels_x) - 1)),
      hr_ci_u = c(1, rep(NA_real_, length(levels_x) - 1)),
      p_value = c(1, rep(NA_real_, length(levels_x) - 1)),
      exposure_class = if (length(levels_x) == 2) "Binary" else paste0("Categorical (", length(levels_x), " levels)"),
      formula = formula_txt,
      check.names = FALSE
    )
    if (length(levels_x) > 1) {
      exposure_rows <- tidy_fit[seq_len(length(levels_x) - 1L), , drop = FALSE]
      if (nrow(exposure_rows) != length(levels_x) - 1L) stop("Failed to recover categorical exposure coefficients from the Cox model.", call. = FALSE)
      result_df$hr[-1] <- exposure_rows$estimate
      result_df$hr_ci_l[-1] <- exposure_rows$conf.low
      result_df$hr_ci_u[-1] <- exposure_rows$conf.high
      result_df$p_value[-1] <- exposure_rows$p.value
    }
    model_results[[model_name]] <- result_df
    fit_results[[model_name]] <- fit
  }

  # Build tidy result rows
  result_tidy <- do.call(rbind, model_results)
  rownames(result_tidy) <- NULL

  # Return result object
  out <- structure(list(
    result = NULL,
    result_tidy = result_tidy,
    data_info = list(
      n_total = prep$n_total,
      n_after_followup = prep$n_after_followup,
      n_after_complete_case = prep$n_after_complete_case,
      min_followup_time = prep$min_followup_time,
      event_value = prep$event_value,
      p_fmt = p_fmt
    ),
    fit = fit_results
  ), class = "leo_cox")
  out$result <- leo_cox_format(out, style = "wide")
  leo.basic::leo_log("Cox regression completed for {prep$exposure_name} -> {prep$outcome_name} with {length(prep$models)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  if (is.character(simplify)) return(leo_cox_format(out, style = simplify))
  return(out)
}

#' Format a `leo_cox` result
#'
#' @rdname leo_cox
#' @param x Result returned by `leo_cox()`.
#' @param style One of `"wide"`, `"tidy"`, or `"gtsummary"`.
#' @section Formatting output:
#' `leo_cox_format()` returns:
#' - `"wide"`: a wide summary data frame with one set of HR / CI / p-value
#'   columns per model.
#' - `"tidy"`: the row-level tidy result table stored in `x$result_tidy`, with
#'   formatted display columns added.
#' - `"gtsummary"`: a `gtsummary` regression table, or a merged `gtsummary`
#'   table when multiple models are present.
#'
#' @export
leo_cox_format <- function(x, style = "wide") {
  if (!inherits(x, "leo_cox")) stop("x must be a leo_cox result.", call. = FALSE)
  if (!style %in% c("wide", "tidy", "gtsummary")) stop("style must be one of 'wide', 'tidy', or 'gtsummary'.", call. = FALSE)
  switch(style,
    wide = {
      if (!is.null(x$result)) return(x$result)
      result_tidy <- x$result_tidy
      model_ids <- unique(result_tidy$model)
      result_wide <- result_tidy[result_tidy$model == model_ids[1], c("row_id", "exposure", "outcome", "case_n", "control_n", "person_year", "exposure_class"), drop = FALSE]
      names(result_wide) <- c("row_id", "Exposure", "Outcome", "Case N", "Control N", "Person-years", "Class")
      p_fmt <- if (!is.null(x$data_info$p_fmt)) x$data_info$p_fmt else "threshold"
      for (model_id in model_ids) {
        model_df <- result_tidy[result_tidy$model == model_id, c("row_id", "hr", "hr_ci_l", "hr_ci_u", "p_value"), drop = FALSE]
        model_df$hr <- round(model_df$hr, 3)
        model_df$ci_95 <- ifelse(is.na(model_df$hr_ci_l) | is.na(model_df$hr_ci_u), "NA", paste0(sprintf("%.3f", round(model_df$hr_ci_l, 3)), ", ", sprintf("%.3f", round(model_df$hr_ci_u, 3))))
        model_df$p_value <- vapply(model_df$p_value, function(p) .format_p_value(p, fmt = p_fmt), character(1))
        model_df <- model_df[, c("row_id", "hr", "ci_95", "p_value"), drop = FALSE]
        names(model_df) <- c("row_id", paste(model_id, "HR"), paste(model_id, "95% CI"), paste(model_id, "P value"))
        result_wide <- merge(result_wide, model_df, by = "row_id", all.x = TRUE, sort = FALSE)
        result_wide <- result_wide[order(result_wide$row_id), , drop = FALSE]
      }
      result_wide$row_id <- NULL
      rownames(result_wide) <- NULL
      result_wide
    },
    tidy = {
      result_tidy <- x$result_tidy
      result_out <- result_tidy[, c("model", "exposure", "outcome", "level", "case_n", "control_n", "person_year"), drop = FALSE]
      p_fmt <- if (!is.null(x$data_info$p_fmt)) x$data_info$p_fmt else "threshold"
      result_out$HR <- round(result_tidy$hr, 3)
      result_out$`95% CI` <- ifelse(is.na(result_tidy$hr_ci_l) | is.na(result_tidy$hr_ci_u), "NA", paste0(sprintf("%.3f", round(result_tidy$hr_ci_l, 3)), ", ", sprintf("%.3f", round(result_tidy$hr_ci_u, 3))))
      result_out$`P value` <- vapply(result_tidy$p_value, function(p) .format_p_value(p, fmt = p_fmt), character(1))
      result_out$Class <- result_tidy$exposure_class
      names(result_out)[1:7] <- c("Model", "Exposure", "Outcome", "Level", "Case N", "Control N", "Person-years")
      rownames(result_out) <- NULL
      result_out
    },
    gtsummary = {
      if (is.null(x$fit) || length(x$fit) == 0) {
        stop("x does not contain fitted models required for style = 'gtsummary'.", call. = FALSE)
      }
      if (!requireNamespace("gtsummary", quietly = TRUE) || !requireNamespace("broom.helpers", quietly = TRUE)) {
        stop("style = 'gtsummary' requires packages 'gtsummary' and 'broom.helpers'.", call. = FALSE)
      }
      tbls <- lapply(x$fit, function(fit) gtsummary::tbl_regression(fit, exponentiate = TRUE))
      if (length(tbls) == 1) return(tbls[[1]])
      gtsummary::tbl_merge(tbls = tbls, tab_spanner = paste0("**", names(tbls), "**"))
    }
  )
}

#' Cox interaction analyses
#'
#' `r lifecycle::badge('experimental')`
#'
#' `leo_cox_interaction()` tests whether the association between an exposure and
#' an incident outcome differs by a candidate interaction variable. It compares
#' nested Cox models with and without the interaction term using
#' `stats::anova(..., test = "Chisq")`.
#'
#' `leo_cox_add_interaction()` focuses on the binary-binary setting and reports
#' additive interaction summaries from a Cox model, including `RERI`, `AP`, and
#' `S`, together with the multiplicative interaction estimate. When the original
#' binary coding does not follow the reference-group convention expected by
#' `interactionR`, `leo_cox_add_interaction()` can either keep the original
#' double-negative group as the reference or manually recode the two binary
#' variables before handing the fitted interaction model to `interactionR`.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, interaction variable, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_inter Character scalar giving the interaction variable column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list of covariate column-name vectors.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#' @param event_value Value in the event column that indicates incident events.
#' @param reference_group Reference-group strategy for additive interaction.
#'   Use `"auto"` to follow the preventive-exposure recode convention described
#'   by `interactionR`, or `"double_negative"` to keep the original
#'   double-negative group as the reference.
#' @param x_exp_type Exposure type handling for `x_exp`.
#' @param x_inter_type Interaction-variable type handling for `x_inter`.
#' @param verbose Logical; print progress messages.
#' @importFrom cli cat_rule
#' @importFrom dplyr group_by summarise
#' @importFrom survival Surv coxph
#'
#' @return A `leo_cox_interaction` object containing a display table in `$result`,
#'   the fitted no-interaction Cox models in `$fit_main`, the fitted
#'   interaction Cox models in `$fit_inter`, and the corresponding model
#'   formulas in `$formula_main` and `$formula_inter`.
#' @export
#' @examples
#' lung_df <- stats::na.omit(
#'   dplyr::transmute(
#'     survival::lung,
#'     outcome = as.integer(status == 2),
#'     outcome_censor = time / 365.25,
#'     age = age,
#'     sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'     ecog_group = factor(ph.ecog, levels = 0:3, labels = c("ECOG0", "ECOG1", "ECOG2", "ECOG3"))
#'   )
#' )
#'
#' model <- list("Crude" = NULL, "Model A" = c("ecog_group"))
#' leo_cox_interaction(
#'   df = lung_df, y_out = c("outcome", "outcome_censor"),
#'   x_exp = "age", x_inter = "sex", x_cov = model, verbose = FALSE
#' )$result
leo_cox_interaction <- function(df, y_out, x_exp, x_inter, x_cov = NULL,
                                event_value = 1, min_followup_time = 0,
                                x_exp_type = "auto", x_inter_type = "auto", verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Interaction", col = "blue")
  if (!is.character(x_inter) || length(x_inter) != 1) stop("x_inter must be a single column name.", call. = FALSE)
  if (!x_inter %in% names(df)) stop("x_inter must exist in df.", call. = FALSE)
  if (!x_exp_type %in% c("auto", "continuous", "categorical")) stop("x_exp_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)
  if (!x_inter_type %in% c("auto", "continuous", "categorical")) stop("x_inter_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  prep_models <- lapply(model_list, function(covariates) unique(c(if (is.null(covariates)) character(0) else covariates, x_inter)))
  prep <- .prepare_cox_regression_data(
    df = df,
    y_out = y_out,
    x_exp = x_exp,
    x_cov = prep_models,
    min_followup_time = min_followup_time,
    event_value = event_value,
    verbose = verbose
  )
  if (prep$n_after_followup == 0) stop("No rows remain after follow-up filtering.", call. = FALSE)
  if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering.", call. = FALSE)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before Cox interaction analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_cox_interaction(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }

  exposure_type <- .check_var_type(prep$data$exposure, var_name = x_exp, var_type = x_exp_type, verbose = verbose)
  interaction_type <- .check_var_type(prep$data[[x_inter]], var_name = x_inter, var_type = x_inter_type, verbose = verbose)
  result_rows <- list()
  fit_main <- list()
  fit_inter <- list()
  formula_main <- stats::setNames(character(length(model_list)), names(model_list))
  formula_inter <- stats::setNames(character(length(model_list)), names(model_list))
  for (model_name in names(model_list)) {
    covariates <- if (is.null(model_list[[model_name]])) character(0) else setdiff(model_list[[model_name]], x_inter)
    leo.basic::leo_log("Testing interaction in {model_name} with {x_inter}{if (length(covariates) > 0) paste0(' adjusted for ', paste(covariates, collapse = ', ')) else ''}.", verbose = verbose)
    model_df <- prep$data[, c("event", "time", "exposure", x_inter, covariates), drop = FALSE]

    if (exposure_type == "categorical") {
      model_df$exposure <- droplevels(factor(model_df$exposure))
      if (nlevels(model_df$exposure) < 2) stop("x_exp must contain at least 2 exposure levels after filtering.", call. = FALSE)
    } else {
      model_df$exposure <- suppressWarnings(as.numeric(model_df$exposure))
      if (anyNA(model_df$exposure)) stop("x_exp_type = 'continuous' requires a numeric or integer exposure.", call. = FALSE)
    }
    if (interaction_type == "categorical") {
      model_df[[x_inter]] <- droplevels(factor(model_df[[x_inter]]))
      if (nlevels(model_df[[x_inter]]) < 2) stop("x_inter must contain at least 2 levels after filtering.", call. = FALSE)
    } else {
      model_df[[x_inter]] <- suppressWarnings(as.numeric(model_df[[x_inter]]))
      if (anyNA(model_df[[x_inter]])) stop("x_inter_type = 'continuous' requires a numeric or integer interaction variable.", call. = FALSE)
    }

    rhs_main <- unique(c("exposure", paste0("`", x_inter, "`"), if (length(covariates) > 0) paste0("`", covariates, "`") else character(0)))
    rhs_inter <- unique(c(paste0("exposure * `", x_inter, "`"), if (length(covariates) > 0) paste0("`", covariates, "`") else character(0)))
    formula_main_fit <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs_main, collapse = " + ")))
    formula_inter_fit <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs_inter, collapse = " + ")))
    fit_main[[model_name]] <- survival::coxph(formula = formula_main_fit, data = model_df)
    fit_inter[[model_name]] <- survival::coxph(formula = formula_inter_fit, data = model_df)
    formula_main[[model_name]] <- Reduce(paste, deparse(formula_main_fit))
    formula_inter[[model_name]] <- Reduce(paste, deparse(formula_inter_fit))
    anova_res <- stats::anova(fit_main[[model_name]], fit_inter[[model_name]], test = "Chisq")
    p_col <- grep("^P", names(anova_res), value = TRUE)[1]
    p_interaction <- if (is.na(p_col)) NA_real_ else unname(anova_res[2, p_col])
    result_rows[[model_name]] <- data.frame(
      model = model_name,
      exposure = x_exp,
      interaction = x_inter,
      n = nrow(model_df),
      case_n = sum(model_df$event == 1, na.rm = TRUE),
      control_n = sum(model_df$event == 0, na.rm = TRUE),
      person_year = sum(model_df$time, na.rm = TRUE),
      interaction_df = unname(anova_res[2, "Df"]),
      p_interaction = p_interaction,
      exposure_class = if (is.factor(model_df$exposure)) if (nlevels(model_df$exposure) == 2) "Binary" else paste0("Categorical (", nlevels(model_df$exposure), " levels)") else "Continuous",
      interaction_class = if (is.factor(model_df[[x_inter]])) if (nlevels(model_df[[x_inter]]) == 2) "Binary" else paste0("Categorical (", nlevels(model_df[[x_inter]]), " levels)") else "Continuous",
      check.names = FALSE
    )
  }
  result_df <- do.call(rbind, result_rows)
  rownames(result_df) <- NULL
  result <- result_df[, c("model", "exposure", "interaction", "n", "case_n", "control_n", "person_year", "interaction_df", "p_interaction", "exposure_class", "interaction_class"), drop = FALSE]
  names(result) <- c("Model", "Exposure", "Interaction", "N", "Case N", "Control N", "Person-years", "Interaction DF", "P for interaction", "Exposure class", "Interaction class")
  result$`P for interaction` <- vapply(result$`P for interaction`, .format_p_value, character(1))
  out <- structure(list(result = result, fit_main = fit_main, fit_inter = fit_inter, formula_main = formula_main, formula_inter = formula_inter), class = "leo_cox_interaction")
  leo.basic::leo_log("Cox interaction analysis completed for {x_exp} x {x_inter} with {length(model_list)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' @rdname leo_cox_interaction
#' @return A `leo_cox_add_interaction` object containing a display table in
#'   `$result`, the fitted no-interaction Cox models in `$fit_main`, the fitted
#'   interaction Cox models in `$fit_inter`, the corresponding model formulas in
#'   `$formula_main` and `$formula_inter`, and the raw `interactionR` backend
#'   objects in `$backend`.
#' @export
#' @examples
#' if (requireNamespace("interactionR", quietly = TRUE)) {
#'   lung_df <- stats::na.omit(
#'     dplyr::transmute(
#'       survival::lung,
#'       outcome = as.integer(status == 2),
#'       outcome_censor = time / 365.25,
#'       smoking = factor(age > median(age), levels = c(FALSE, TRUE), labels = c("Low", "High")),
#'       sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'       ecog_group = factor(ph.ecog, levels = 0:3, labels = c("ECOG0", "ECOG1", "ECOG2", "ECOG3"))
#'     )
#'   )
#'
#'   model <- list("Crude" = NULL, "Model A" = c("ecog_group"))
#'   leo_cox_add_interaction(
#'     df = lung_df, y_out = c("outcome", "outcome_censor"),
#'     x_exp = "smoking", x_inter = "sex", x_cov = model, verbose = FALSE
#'   )$result
#' }
leo_cox_add_interaction <- function(df, y_out, x_exp, x_inter, x_cov = NULL, event_value = 1, min_followup_time = 0, reference_group = c("auto", "double_negative"), verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Additive Interaction", col = "blue")
  # interactionR stays in Suggests, so keep runtime checks and explicit namespace calls.
  if (!requireNamespace("interactionR", quietly = TRUE)) stop("leo_cox_add_interaction() requires package 'interactionR'. Please install it first.", call. = FALSE)
  if (!is.character(x_inter) || length(x_inter) != 1) stop("x_inter must be a single column name.", call. = FALSE)
  if (!x_inter %in% names(df)) stop("x_inter must exist in df.", call. = FALSE)
  reference_strategy <- match.arg(reference_group)

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  prep_models <- lapply(model_list, function(covariates) unique(c(if (is.null(covariates)) character(0) else covariates, x_inter)))
  prep <- .prepare_cox_regression_data(
    df = df,
    y_out = y_out,
    x_exp = x_exp,
    x_cov = prep_models,
    min_followup_time = min_followup_time,
    event_value = event_value,
    verbose = verbose
  )
  if (prep$n_after_followup == 0) stop("No rows remain after follow-up filtering.", call. = FALSE)
  if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering.", call. = FALSE)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before additive interaction analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_cox_add_interaction(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }

  result_rows <- list()
  fit_main <- list()
  fit_inter <- list()
  formula_main <- stats::setNames(character(length(model_list)), names(model_list))
  formula_inter <- stats::setNames(character(length(model_list)), names(model_list))
  backend <- list()
  for (model_name in names(model_list)) {
    covariates <- if (is.null(model_list[[model_name]])) character(0) else setdiff(model_list[[model_name]], x_inter)
    leo.basic::leo_log("Fitting additive interaction in {model_name} with binary {x_exp} and {x_inter}{if (length(covariates) > 0) paste0(' adjusted for ', paste(covariates, collapse = ', ')) else ''}.", verbose = verbose)
    model_df <- prep$data[, c("event", "time", "exposure", x_inter, covariates), drop = FALSE]

    exposure_factor <- if (is.factor(model_df$exposure)) droplevels(model_df$exposure) else droplevels(factor(model_df$exposure))
    if (nlevels(exposure_factor) != 2) stop("x_exp must contain exactly 2 levels after filtering for leo_cox_add_interaction().", call. = FALSE)
    interaction_factor <- if (is.factor(model_df[[x_inter]])) droplevels(model_df[[x_inter]]) else droplevels(factor(model_df[[x_inter]]))
    if (nlevels(interaction_factor) != 2) stop("x_inter must contain exactly 2 levels after filtering for leo_cox_add_interaction().", call. = FALSE)

    exposure_bin <- as.integer(exposure_factor == levels(exposure_factor)[2])
    interaction_bin <- as.integer(interaction_factor == levels(interaction_factor)[2])
    joint_counts <- table(exposure_bin, interaction_bin)
    if (!all(dim(joint_counts) == c(2, 2)) || any(joint_counts == 0)) stop("All four joint exposure groups must be present after filtering for leo_cox_add_interaction().", call. = FALSE)

    analysis_df <- data.frame(
      event = model_df$event,
      time = model_df$time,
      check.names = FALSE
    )
    analysis_df[[x_exp]] <- exposure_bin
    analysis_df[[x_inter]] <- interaction_bin
    if (length(covariates) > 0) analysis_df[covariates] <- model_df[covariates]

    rhs_main <- unique(c(paste0("`", x_exp, "`"), paste0("`", x_inter, "`"), if (length(covariates) > 0) paste0("`", covariates, "`") else character(0)))
    rhs_inter <- unique(c(paste0("`", x_exp, "` * `", x_inter, "`"), if (length(covariates) > 0) paste0("`", covariates, "`") else character(0)))
    formula_main_fit <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs_main, collapse = " + ")))
    formula_inter_fit <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs_inter, collapse = " + ")))
    fit_main_current <- survival::coxph(formula = formula_main_fit, data = analysis_df)
    fit_inter_current <- survival::coxph(formula = formula_inter_fit, data = analysis_df)

    recode_exposure <- FALSE
    recode_interaction <- FALSE
    if (reference_strategy == "auto") {
      coef_names <- names(stats::coef(fit_inter_current))
      beta1 <- coef_names[coef_names == x_exp]
      beta2 <- coef_names[coef_names == x_inter]
      beta3 <- coef_names[coef_names %in% c(paste0(x_exp, ":", x_inter), paste0(x_inter, ":", x_exp))]
      if (length(beta1) != 1 || length(beta2) != 1 || length(beta3) != 1) {
        stop("Failed to recover additive interaction coefficients from the fitted Cox model. Available coefficients: ", paste(coef_names, collapse = ", "), call. = FALSE)
      }
      b1 <- stats::coef(fit_inter_current)[beta1]
      b2 <- stats::coef(fit_inter_current)[beta2]
      b3 <- stats::coef(fit_inter_current)[beta3]
      if ((exp(b1) < 1) || (exp(b2) < 1)) {
        ref_cat <- stats::setNames(
          c(unname(exp(b1)), unname(exp(b2)), unname(exp(b1 + b2 + b3))),
          c("OR10", "OR01", "OR11")
        )
        ref_cat <- names(ref_cat)[which.min(ref_cat)]
        recode_exposure <- ref_cat %in% c("OR10", "OR11")
        recode_interaction <- ref_cat %in% c("OR01", "OR11")
      }
    }

    if (recode_exposure) analysis_df[[x_exp]] <- 1L - analysis_df[[x_exp]]
    if (recode_interaction) analysis_df[[x_inter]] <- 1L - analysis_df[[x_inter]]
    if (recode_exposure || recode_interaction) {
      fit_main_current <- survival::coxph(formula = formula_main_fit, data = analysis_df)
      fit_inter_current <- survival::coxph(formula = formula_inter_fit, data = analysis_df)
    }
    fit_main[[model_name]] <- fit_main_current
    fit_inter[[model_name]] <- fit_inter_current
    formula_main[[model_name]] <- Reduce(paste, deparse(formula_main_fit))
    formula_inter[[model_name]] <- Reduce(paste, deparse(formula_inter_fit))

    reference_group_label <- paste0(
      x_exp, "=", levels(exposure_factor)[if (recode_exposure) 2 else 1],
      ", ",
      x_inter, "=", levels(interaction_factor)[if (recode_interaction) 2 else 1]
    )
    if (reference_strategy == "double_negative") {
      reference_group_label <- paste0(
        x_exp, "=", levels(exposure_factor)[1],
        ", ",
        x_inter, "=", levels(interaction_factor)[1]
      )
    }

    anova_res <- stats::anova(fit_main[[model_name]], fit_inter[[model_name]], test = "Chisq")
    p_col <- grep("^P", names(anova_res), value = TRUE)[1]
    p_interaction <- if (is.na(p_col)) NA_real_ else unname(anova_res[2, p_col])
    backend_df <- analysis_df
    backend_exp_name <- "leo_exposure_internal"
    backend_inter_name <- "leo_interaction_internal"
    backend_df[[backend_exp_name]] <- backend_df[[x_exp]]
    backend_df[[backend_inter_name]] <- backend_df[[x_inter]]
    rhs_backend <- unique(c(
      paste0("`", backend_exp_name, "` * `", backend_inter_name, "`"),
      if (length(covariates) > 0) paste0("`", covariates, "`") else character(0)
    ))
    backend_formula_inter <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(rhs_backend, collapse = " + ")))
    backend_fit <- survival::coxph(formula = backend_formula_inter, data = backend_df)
    backend_obj <- interactionR::interactionR(
      model = backend_fit,
      exposure_names = c(backend_exp_name, backend_inter_name),
      ci.type = "delta",
      ci.level = 0.95,
      em = FALSE,
      recode = FALSE
    )
    backend[[model_name]] <- backend_obj
    backend_df <- backend_obj$dframe
    mult_row <- backend_df[backend_df$Measures == "Multiplicative scale", , drop = FALSE]
    reri_row <- backend_df[backend_df$Measures == "RERI", , drop = FALSE]
    ap_row <- backend_df[backend_df$Measures == "AP", , drop = FALSE]
    s_row <- backend_df[backend_df$Measures %in% c("SI", "S"), , drop = FALSE]

    result_rows[[model_name]] <- data.frame(
      model = model_name,
      exposure = x_exp,
      interaction = x_inter,
      n = nrow(analysis_df),
      case_n = sum(analysis_df$event == 1, na.rm = TRUE),
      control_n = sum(analysis_df$event == 0, na.rm = TRUE),
      person_year = sum(analysis_df$time, na.rm = TRUE),
      reference_group = reference_group_label,
      recode_applied = recode_exposure || recode_interaction,
      mult = mult_row$Estimates[1],
      mult_l = mult_row$CI.ll[1],
      mult_u = mult_row$CI.ul[1],
      p_interaction = p_interaction,
      reri = reri_row$Estimates[1],
      reri_l = reri_row$CI.ll[1],
      reri_u = reri_row$CI.ul[1],
      ap = ap_row$Estimates[1],
      ap_l = ap_row$CI.ll[1],
      ap_u = ap_row$CI.ul[1],
      s = s_row$Estimates[1],
      s_l = s_row$CI.ll[1],
      s_u = s_row$CI.ul[1],
      check.names = FALSE
    )
  }

  result_df <- do.call(rbind, result_rows)
  rownames(result_df) <- NULL
  format_ci <- function(lower, upper) ifelse(is.na(lower) | is.na(upper), "NA", paste0(sprintf("%.3f", round(lower, 3)), ", ", sprintf("%.3f", round(upper, 3))))
  result <- result_df[, c("model", "exposure", "interaction", "n", "case_n", "control_n", "person_year", "reference_group", "recode_applied"), drop = FALSE]
  result$`Multiplicative interaction` <- round(result_df$mult, 3)
  result$`Multiplicative 95% CI` <- format_ci(result_df$mult_l, result_df$mult_u)
  result$`P for interaction` <- vapply(result_df$p_interaction, .format_p_value, character(1))
  result$RERI <- round(result_df$reri, 3)
  result$`RERI 95% CI` <- format_ci(result_df$reri_l, result_df$reri_u)
  result$AP <- round(result_df$ap, 3)
  result$`AP 95% CI` <- format_ci(result_df$ap_l, result_df$ap_u)
  result$S <- round(result_df$s, 3)
  result$`S 95% CI` <- format_ci(result_df$s_l, result_df$s_u)
  names(result)[1:9] <- c("Model", "Exposure", "Interaction", "N", "Case N", "Control N", "Person-years", "Reference group", "Recode applied")
  out <- structure(list(
    result = result,
    fit_main = fit_main,
    fit_inter = fit_inter,
    formula_main = formula_main,
    formula_inter = formula_inter,
    backend = backend
  ), class = "leo_cox_add_interaction")
  leo.basic::leo_log("Cox additive interaction analysis completed for {x_exp} x {x_inter} with {length(model_list)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' Cox subgroup analysis
#'
#' `r lifecycle::badge('experimental')`
#'
#' `leo_cox_subgroup()` runs `leo_cox()` within each level of one or more
#' subgroup variables and appends a `P for interaction` from
#' `leo_cox_interaction()`. Optionally, it also calculates a summary-statistic
#' heterogeneity P value using `leo_heterogeneity_p()`.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, subgroup variables, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_subgroup Character vector giving one or more subgroup column names.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list of covariate column-name vectors.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#' @param event_value Value in the event column that indicates incident events.
#' @param x_exp_type Exposure type handling for `x_exp`.
#' @param x_subgroup_type Type handling for `x_subgroup`. Subgroup analyses are expected to use categorical grouping variables.
#' @param add_interaction Logical; whether to append `P for interaction` from `leo_cox_interaction()`.
#' @param add_heterogeneity Logical; whether to calculate a heterogeneity P value from subgroup summary statistics.
#' @param verbose Logical; print progress messages.
#'
#' @return A `leo_cox_subgroup` object containing a subgroup result table in
#'   `$result` and, when requested, interaction test objects in `$interaction`.
#'   When `x_subgroup` has length 1, `$interaction` is a single
#'   `leo_cox_interaction` object; otherwise it is a named list of
#'   `leo_cox_interaction` objects keyed by subgroup variable.
#' @export
#' @examples
#' lung_df <- stats::na.omit(
#'   dplyr::transmute(
#'     survival::lung,
#'     outcome = as.integer(status == 2),
#'     outcome_censor = time / 365.25,
#'     age = age,
#'     sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'     ecog_group = factor(ph.ecog, levels = 0:3, labels = c("ECOG0", "ECOG1", "ECOG2", "ECOG3"))
#'   )
#' )
#'
#' model <- list("Crude" = NULL, "Model A" = c("ecog_group"))
#' res_sub <- leo_cox_subgroup(
#'   df = lung_df, y_out = c("outcome", "outcome_censor"),
#'   x_exp = "age", x_subgroup = "sex", x_cov = model, verbose = FALSE
#' )
#' res_sub$result
#' leo_cox_subgroup_format(res_sub, style = "wide")
leo_cox_subgroup <- function(df, y_out, x_exp, x_subgroup, x_cov = NULL,
                             event_value = 1, min_followup_time = 0,
                             x_exp_type = "auto", x_subgroup_type = "categorical",
                             add_interaction = TRUE, add_heterogeneity = FALSE, verbose = TRUE) {
  t0 <- Sys.time()
  if (verbose) cli::cat_rule("Cox Subgroup", col = "blue")
  if (!is.character(x_subgroup) || length(x_subgroup) < 1) stop("x_subgroup must be a character vector of subgroup column names.", call. = FALSE)
  if (!all(x_subgroup %in% names(df))) stop("All x_subgroup columns must exist in df.", call. = FALSE)
  if (!x_exp_type %in% c("auto", "continuous", "categorical")) stop("x_exp_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)
  if (!x_subgroup_type %in% c("auto", "continuous", "categorical")) stop("x_subgroup_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)

  base_models <- .normalize_model_list(x_cov, df_colnames = names(df))
  subgroup_rows <- list()
  interaction_fits <- list()
  subgroup_index <- 1L
  for (subgroup_var in x_subgroup) {
    subgroup_models <- lapply(base_models, function(covariates) setdiff(if (is.null(covariates)) character(0) else covariates, subgroup_var))
    prep_models <- lapply(subgroup_models, function(covariates) unique(c(covariates, subgroup_var)))
    prep <- .prepare_cox_regression_data(
      df = df,
      y_out = y_out,
      x_exp = x_exp,
      x_cov = prep_models,
      min_followup_time = min_followup_time,
      event_value = event_value,
      verbose = verbose
    )
    if (prep$n_after_followup == 0 || nrow(prep$data) == 0) next
    subgroup_type <- .check_var_type(prep$data[[subgroup_var]], var_name = subgroup_var, var_type = x_subgroup_type, verbose = verbose)
    if (subgroup_type != "categorical") stop("Subgroup analysis requires a categorical grouping variable. Convert ", subgroup_var, " to factor() first if needed.", call. = FALSE)
    prep$data[[subgroup_var]] <- droplevels(factor(prep$data[[subgroup_var]]))
    interaction_fit <- NULL
    if (add_interaction) {
      interaction_fit <- tryCatch(
        leo_cox_interaction(
          df = df,
          y_out = y_out,
          x_exp = x_exp,
          x_inter = subgroup_var,
          x_cov = subgroup_models,
          event_value = event_value,
          min_followup_time = min_followup_time,
          x_exp_type = x_exp_type,
          x_inter_type = "categorical",
          verbose = FALSE
        ),
        error = function(e) {
          leo.basic::leo_log("Skipping interaction test for subgroup '{subgroup_var}': {e$message}", level = "warning", verbose = verbose)
          return(NULL)
        }
      )
    }
    interaction_fits[[subgroup_var]] <- interaction_fit
    for (subgroup_level in levels(prep$data[[subgroup_var]])) {
      sub_df <- prep$data[prep$data[[subgroup_var]] == subgroup_level, c("event", "time", "exposure", unique(unlist(subgroup_models, use.names = FALSE))), drop = FALSE]
      if (nrow(sub_df) == 0) next
      if (length(unique(stats::na.omit(sub_df$event))) < 2) {
        leo.basic::leo_log("Skipping subgroup '{subgroup_var} = {subgroup_level}' because the outcome has fewer than 2 observed states.", level = "warning", verbose = verbose)
        next
      }
      if (length(unique(stats::na.omit(sub_df$exposure))) < 2) {
        leo.basic::leo_log("Skipping subgroup '{subgroup_var} = {subgroup_level}' because exposure has fewer than 2 observed values.", level = "warning", verbose = verbose)
        next
      }
      fit_sub <- tryCatch(
        leo_cox(
          df = sub_df,
          y_out = c("event", "time"),
          x_exp = "exposure",
          x_cov = subgroup_models,
          event_value = 1,
          min_followup_time = min_followup_time,
          x_exp_type = x_exp_type,
          simplify = FALSE,
          verbose = FALSE
        ),
        error = function(e) {
          leo.basic::leo_log("Skipping subgroup '{subgroup_var} = {subgroup_level}' because Cox fitting failed: {e$message}", level = "warning", verbose = verbose)
          return(NULL)
        }
      )
      if (is.null(fit_sub)) next
      sub_rows <- fit_sub$result_tidy
      sub_rows$exposure <- x_exp
      sub_rows$outcome <- y_out[1]
      sub_rows$level[is.na(sub_rows$level) & sub_rows$exposure_class == "Continuous"] <- "Per unit increase"
      sub_rows$subgroup <- subgroup_var
      sub_rows$subgroup_level <- as.character(subgroup_level)
      sub_rows$subgroup_n <- nrow(sub_df)
      sub_rows$row_key <- paste(subgroup_var, subgroup_level, sub_rows$row_id, sep = "___")
      sub_rows$display_id <- seq.int(from = subgroup_index, length.out = nrow(sub_rows))
      subgroup_index <- subgroup_index + nrow(sub_rows)
      sub_rows$p_interaction_display <- if (is.null(interaction_fit)) NA_character_ else interaction_fit$result$`P for interaction`[match(sub_rows$model, interaction_fit$result$Model)]
      sub_rows$p_heterogeneity <- NA_real_
      subgroup_rows[[paste0(subgroup_var, "___", subgroup_level)]] <- sub_rows
    }
  }
  if (length(subgroup_rows) == 0) stop("No subgroup results were generated.", call. = FALSE)
  result_df <- do.call(rbind, subgroup_rows)
  rownames(result_df) <- NULL

  if (add_heterogeneity) {
    for (model_name in unique(result_df$model)) {
      for (subgroup_var in unique(result_df$subgroup)) {
        for (row_id in unique(result_df$row_id[result_df$subgroup == subgroup_var])) {
          idx <- which(result_df$model == model_name & result_df$subgroup == subgroup_var & result_df$row_id == row_id)
          idx <- idx[is.finite(result_df$hr[idx]) & is.finite(result_df$p_value[idx]) & !is.na(result_df$p_value[idx]) & !grepl("(Ref)", result_df$exposure[idx], fixed = TRUE)]
          if (length(idx) < 2) next
          hetero_res <- utils::capture.output(
            hetero_obj <- leo_heterogeneity_p(
              hrs = result_df$hr[idx],
              p_values = result_df$p_value[idx],
              subgroup_names = result_df$subgroup_level[idx]
            )
          )
          result_df$p_heterogeneity[idx] <- hetero_obj$p_value_heterogeneity
        }
      }
    }
  }

  result <- result_df[, c("subgroup", "subgroup_level", "subgroup_n", "model", "exposure", "outcome", "level", "case_n", "control_n", "person_year"), drop = FALSE]
  result$HR <- round(result_df$hr, 3)
  result$`95% CI` <- ifelse(is.na(result_df$hr_ci_l) | is.na(result_df$hr_ci_u), "NA", paste0(sprintf("%.3f", round(result_df$hr_ci_l, 3)), ", ", sprintf("%.3f", round(result_df$hr_ci_u, 3))))
  result$`P value` <- vapply(result_df$p_value, .format_p_value, character(1))
  if (add_interaction) result$`P for interaction` <- ifelse(is.na(result_df$p_interaction_display), "NA", result_df$p_interaction_display)
  if (add_heterogeneity) result$`P for heterogeneity` <- vapply(result_df$p_heterogeneity, .format_p_value, character(1))
  names(result)[1:10] <- c("Subgroup", "Level", "N", "Model", "Exposure", "Outcome", "Exposure level", "Case N", "Control N", "Person-years")
  rownames(result) <- NULL

  interaction_out <- interaction_fits
  if (length(interaction_out) == 1) interaction_out <- interaction_out[[1]]
  out <- structure(list(
    result = result,
    interaction = interaction_out
  ), class = "leo_cox_subgroup")
  leo.basic::leo_log("Cox subgroup analysis completed for {x_exp} across {length(x_subgroup)} subgroup variable(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' Format a `leo_cox_subgroup` result
#'
#' @rdname leo_cox_subgroup
#' @param x Result returned by `leo_cox_subgroup()`.
#' @param style One of `"wide"` or `"tidy"`.
#' @export
leo_cox_subgroup_format <- function(x, style = "wide") {
  if (!inherits(x, "leo_cox_subgroup")) stop("x must be a leo_cox_subgroup result.", call. = FALSE)
  if (!style %in% c("wide", "tidy")) stop("style must be one of 'wide' or 'tidy'.", call. = FALSE)
  switch(style,
    wide = {
      result_long <- x$result
      base_cols <- c("Subgroup", "Level", "N", "Exposure", "Outcome", "Exposure level", "Case N", "Control N", "Person-years")
      model_ids <- unique(result_long$Model)
      result_wide <- result_long[result_long$Model == model_ids[1], base_cols, drop = FALSE]
      for (model_id in model_ids) {
        model_df <- result_long[result_long$Model == model_id, c("Subgroup", "Level", "HR", "95% CI", "P value"), drop = FALSE]
        rename_cols <- c("Subgroup", "Level", paste(model_id, "HR"), paste(model_id, "95% CI"), paste(model_id, "P value"))
        if ("P for interaction" %in% names(result_long)) {
          model_df$`P for interaction` <- result_long$`P for interaction`[result_long$Model == model_id]
          rename_cols <- c(rename_cols, paste(model_id, "P for interaction"))
        }
        if ("P for heterogeneity" %in% names(result_long)) {
          model_df$`P for heterogeneity` <- result_long$`P for heterogeneity`[result_long$Model == model_id]
          rename_cols <- c(rename_cols, paste(model_id, "P for heterogeneity"))
        }
        names(model_df) <- rename_cols
        result_wide <- merge(result_wide, model_df, by = c("Subgroup", "Level"), all.x = TRUE, sort = FALSE)
      }
      rownames(result_wide) <- NULL
      result_wide
    },
    tidy = {
      x$result
    }
  )
}

#' Cox mediation analysis
#'
#' `r lifecycle::badge('experimental')`
#'
#' Performs causal mediation analysis for survival outcomes via `regmedint::regmedint()`.
#' Fits mediator model `x_med ~ x_exp + x_cov` and outcome model
#' `time ~ x_exp + x_med + x_cov` (with optional `x_exp:x_med` interaction).
#'
#' `yreg = "auto"` uses the observed event proportion in the post-filter,
#' complete-case analysis set: Cox for rare outcomes (<=10% events), AFT otherwise.
#' This 10% threshold is an empirical rare-event approximation rule.
#'
#' @param df Data frame with outcome, time, exposure, mediator, and covariates.
#' @param y_out `c(event, time)` column names.
#' @param x_exp Exposure column name.
#' @param x_med Mediator column name.
#' @param x_cov Covariate column names, or `NULL`.
#' @param min_followup_time Minimum follow-up time filter.
#' @param event_value Value indicating an event.
#' @param x_exp_a0 Reference value of `x_exp`. Defaults to first factor level;
#'   **required if x_exp is character (not factor)**. For a binary factor
#'   exposure, either use the default contrast or supply `x_exp_a0` and
#'   `x_exp_a1` together. For numeric `x_exp`, if exactly two unique observed
#'   values remain after filtering, they are used as the default contrast;
#'   otherwise, supply `x_exp_a0` and `x_exp_a1` explicitly.
#' @param x_exp_a1 Contrasted value of `x_exp`. Defaults to second factor level;
#'   **required if x_exp is character (not factor)**. For a binary factor
#'   exposure, either use the default contrast or supply `x_exp_a0` and
#'   `x_exp_a1` together. For numeric `x_exp`, if exactly two unique observed
#'   values remain after filtering, they are used as the default contrast;
#'   otherwise, supply `x_exp_a0` and `x_exp_a1` explicitly.
#' @param x_med_cde Value of `x_med` for CDE evaluation. Defaults to second factor
#'   level (binary) or median (continuous); **required if x_med is character (not factor)**.
#'   For logistic mediators, use the original observed mediator value rather than
#'   an internal 0/1 recode unless the observed values are themselves 0/1.
#' @param x_cov_cond Covariate values for effect evaluation. Named list/vector for
#'   non-numeric covariates; defaults to medians. Required for factor/binary covariates.
#' @param mreg Mediator model: `"auto"`, `"linear"`, or `"logistic"`.
#' @param yreg Outcome model: `"auto"`, `"survCox"`, `"survAFT_weibull"`, or `"survAFT_exp"`.
#' @param interaction Include `x_exp:x_med` interaction? If `TRUE`, the outcome model adds `x_exp:x_med`,
#'   i.e. `time ~ x_exp + x_med + x_exp:x_med + x_cov`.
#' @param verbose Print progress messages.
#'
#' @return `leo_cox_mediation` object with `$result`, `$result_detail`, `$evaluation`, `$fit`.
#' @export
#' @examples
#' if (requireNamespace("regmedint", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 200
#'   age <- rnorm(n, 60, 8)
#'   exposure <- rbinom(n, 1, 0.5)
#'   mediator <- rnorm(n, 0.6 * exposure + 0.02 * (age - 60), 1)
#'   time_event <- rexp(n, rate = exp(-5.8 + 0.55 * exposure + 0.20 * mediator + 0.02 * (age - 60)))
#'   time_censor <- rexp(n, rate = 0.08)
#'   med_df <- data.frame(
#'     outcome = as.integer(time_event <= time_censor),
#'     outcome_censor = pmax(pmin(time_event, time_censor), 0.1),
#'     exposure = exposure, mediator = mediator, age = age
#'   )
#'   res_med <- leo_cox_mediation(
#'     df = med_df, y_out = c("outcome", "outcome_censor"),
#'     x_exp = "exposure", x_med = "mediator", x_cov = "age",
#'     verbose = FALSE
#'   )
#'   res_med$result
#' }
#'
#' # Continuous exposure requires x_exp_a0 and x_exp_a1
#' if (requireNamespace("regmedint", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 200
#'   cont_df <- data.frame(
#'     outcome = rbinom(n, 1, 0.06),
#'     outcome_censor = rexp(n, 0.1),
#'     exposure = rnorm(n), mediator = rnorm(n), age = rnorm(n, 60, 8)
#'   )
#'   res_cont <- leo_cox_mediation(
#'     df = cont_df, y_out = c("outcome", "outcome_censor"),
#'     x_exp = "exposure", x_med = "mediator", x_cov = "age",
#'     x_exp_a0 = -1, x_exp_a1 = 1, verbose = FALSE
#'   )
#'   res_cont$result
#' }
leo_cox_mediation <- function(df, y_out, x_exp, x_med, x_cov = NULL,
                              event_value = 1, min_followup_time = 0,
                              x_exp_a0 = NULL, x_exp_a1 = NULL,
                              x_med_cde = NULL, x_cov_cond = NULL,
                              mreg = c("auto", "linear", "logistic"), 
                              yreg = c("auto", "survCox", "survAFT_weibull", "survAFT_exp"), 
                              interaction = TRUE, verbose = TRUE) {

  # initialization and validation ----
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Mediation", col = "blue")
  if (!requireNamespace("regmedint", quietly = TRUE)) stop("leo_cox_mediation() requires package 'regmedint'. Please install it first.", call. = FALSE)
  if (!is.character(x_med) || length(x_med) != 1) stop("x_med must be a single column name.", call. = FALSE)
  if (!x_med %in% names(df)) stop("x_med must exist in df.", call. = FALSE)
  mreg <- match.arg(mreg, c("auto", "linear", "logistic"))
  yreg <- match.arg(yreg)

  # helper functions ----
  format_num <- function(x) sprintf("%.3f", round(x, 3))
  normalize_factor_like <- function(x) {
    if (is.factor(x)) return(droplevels(x))
    x_chr <- as.character(x)
    if (is.logical(x)) {
      level_order <- c("FALSE", "TRUE")
      level_order <- level_order[level_order %in% x_chr]
    } else {
      level_order <- unique(x_chr)
    }
    factor(x_chr, levels = level_order)
  }
  build_cvar_matrix <- function(cov_df, covariates, x_cov_cond) {
    if (length(covariates) == 0) return(list(covariate_df = NULL, cvar_use = NULL, c_cond_use = NULL, c_cond_label = NA_character_))
    cov_model_df <- cov_df
    binary_numeric_covariates <- covariates[vapply(cov_df, function(x) {
      if (!(is.numeric(x) || is.integer(x))) return(FALSE)
      x_num <- suppressWarnings(as.numeric(x))
      x_num <- sort(unique(stats::na.omit(x_num[is.finite(x_num)])))
      if (length(x_num) != 2) return(FALSE)
      all(abs(x_num - round(x_num)) < 1e-8)
    }, logical(1))]
    for (covariate in covariates) {
      if (is.character(cov_model_df[[covariate]]) || is.logical(cov_model_df[[covariate]])) cov_model_df[[covariate]] <- normalize_factor_like(cov_model_df[[covariate]])
    }
    cov_formula <- stats::reformulate(covariates)
    cov_matrix <- stats::model.matrix(cov_formula, data = cov_model_df)
    cov_matrix <- cov_matrix[, colnames(cov_matrix) != "(Intercept)", drop = FALSE]
    if (is.null(x_cov_cond)) {
      if (!all(vapply(cov_df, function(x) is.numeric(x) || is.integer(x), logical(1)))) {
        stop("Please supply x_cov_cond explicitly when x_cov contains non-numeric covariates for leo_cox_mediation().", call. = FALSE)
      }
      if (length(binary_numeric_covariates) > 0) {
        stop(
          "Please supply x_cov_cond explicitly when x_cov contains binary numeric covariates (",
          paste(binary_numeric_covariates, collapse = ", "),
          "), because the default median can create impossible evaluation values such as 0.5.",
          call. = FALSE
        )
      }
      c_cond_model_df <- as.data.frame(as.list(vapply(cov_df, stats::median, numeric(1), na.rm = TRUE)), stringsAsFactors = FALSE)
      c_cond_label <- paste(paste0(covariates, "=", format_num(unlist(c_cond_model_df[1, covariates, drop = FALSE]))), collapse = "; ")
    } else if (is.numeric(x_cov_cond) && is.null(names(x_cov_cond))) {
      if (!all(vapply(cov_df, function(x) is.numeric(x) || is.integer(x), logical(1)))) {
        stop("For non-numeric covariates, x_cov_cond must be a named list or named vector using original covariate values.", call. = FALSE)
      }
      if (length(x_cov_cond) != length(covariates)) stop("x_cov_cond must have the same length as the selected covariates.", call. = FALSE)
      c_cond_model_df <- as.data.frame(as.list(as.numeric(x_cov_cond)), stringsAsFactors = FALSE)
      names(c_cond_model_df) <- covariates
      c_cond_label <- paste(paste0(covariates, "=", format_num(as.numeric(x_cov_cond))), collapse = "; ")
    } else {
      c_cond_list <- if (is.list(x_cov_cond)) x_cov_cond else as.list(x_cov_cond)
      if (is.null(names(c_cond_list)) || any(names(c_cond_list) == "")) {
        stop("Named x_cov_cond values are required when supplying covariate-specific evaluation settings.", call. = FALSE)
      }
      if (!all(covariates %in% names(c_cond_list))) {
        stop("Named x_cov_cond must include every selected covariate.", call. = FALSE)
      }
      c_cond_model_df <- cov_model_df[1, covariates, drop = FALSE]
      c_cond_label_parts <- character(0)
      for (covariate in covariates) {
        value <- c_cond_list[[covariate]]
        if (length(value) != 1 || is.na(value[1])) stop("x_cov_cond value for ", covariate, " must be a single non-missing value.", call. = FALSE)
        value <- value[1]
        if (is.factor(cov_model_df[[covariate]])) {
          value_chr <- as.character(value)
          if (!value_chr %in% levels(cov_model_df[[covariate]])) {
            stop("x_cov_cond value for ", covariate, " must match one of the observed factor levels.", call. = FALSE)
          }
          c_cond_model_df[[covariate]] <- factor(value_chr, levels = levels(cov_model_df[[covariate]]))
          c_cond_label_parts <- c(c_cond_label_parts, paste0(covariate, "=", value_chr))
        } else {
          value_num <- suppressWarnings(as.numeric(value))
          if (length(value_num) != 1 || is.na(value_num)) stop("x_cov_cond value for ", covariate, " must be numeric.", call. = FALSE)
          c_cond_model_df[[covariate]] <- value_num
          c_cond_label_parts <- c(c_cond_label_parts, paste0(covariate, "=", format_num(value_num)))
        }
      }
      c_cond_label <- paste(c_cond_label_parts, collapse = "; ")
    }
    c_cond_matrix <- stats::model.matrix(cov_formula, data = c_cond_model_df)
    c_cond_matrix <- c_cond_matrix[, colnames(c_cond_matrix) != "(Intercept)", drop = FALSE]
    cvar_use <- colnames(cov_matrix)
    c_cond_use <- as.numeric(c_cond_matrix[1, cvar_use, drop = TRUE])
    names(c_cond_use) <- cvar_use
    return(list(
      covariate_df = as.data.frame(cov_matrix, check.names = FALSE),
      cvar_use = cvar_use,
      c_cond_use = c_cond_use,
      c_cond_label = c_cond_label
    ))
  }

  # data preparation and yreg selection ----
  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))

  # main loop ----
  result_rows <- list()
  fit_results <- list()
  evaluation_rows <- list()
  for (model_name in names(model_list)) {
    covariates <- if (is.null(model_list[[model_name]])) character(0) else model_list[[model_name]]
    prep_covariates <- unique(c(covariates, x_med))
    prep <- .prepare_cox_regression_data(
      df = df,
      y_out = y_out,
      x_exp = x_exp,
      x_cov = prep_covariates,
      min_followup_time = min_followup_time,
      event_value = event_value,
      verbose = verbose
    )
    if (prep$n_after_followup == 0) stop("No rows remain after follow-up filtering for ", model_name, ".", call. = FALSE)
    if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering for ", model_name, ".", call. = FALSE)
    if (prep$n_removed_complete_case > 0) {
      leo.basic::leo_log(
        "For {model_name}, if you want to keep more rows, consider imputing missing values before Cox mediation analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_cox_mediation(df = {df_name}_imputed, ...)`.",
        level = "warning",
        verbose = verbose
      )
    }
    model_df <- prep$data[, c("event", "time", "exposure", x_med, covariates), drop = FALSE]
    event_rate <- sum(model_df$event == 1, na.rm = TRUE) / nrow(model_df)
    yreg_use <- if (yreg == "auto") {
      if (is.finite(event_rate) && event_rate <= 0.1) "survCox" else "survAFT_weibull"
    } else {
      yreg
    }
    if (yreg == "auto") {
      leo.basic::leo_log(
        "In {model_name}, observed event proportion is {format_num(100 * event_rate)}%; auto-selecting `yreg = '{yreg_use}'` for leo_cox_mediation().",
        verbose = verbose
      )
    }
    if (identical(yreg_use, "survCox") && is.finite(event_rate) && event_rate > 0.1) {
      warning(
        "Observed event proportion is ", format_num(100 * event_rate),
        "% in ", model_name,
        ". `regmedint` with `yreg = 'survCox'` relies on the rare-event approximation, so interpret the mediation estimates with caution when the outcome is not rare.",
        call. = FALSE
      )
    }
    a0_use <- x_exp_a0
    a1_use <- x_exp_a1
    m_cde_use <- x_med_cde
    m_cde_internal_use <- x_med_cde
    m_cde_output_use <- x_med_cde
    exposure_contrast_display <- NULL
    mediator_reference_display <- NULL

    ## exposure processing ----
    exposure_raw <- model_df$exposure
    if (is.factor(exposure_raw) || is.character(exposure_raw) || is.logical(exposure_raw)) {
      exposure_factor <- normalize_factor_like(exposure_raw)
      if (nlevels(exposure_factor) != 2) stop("leo_cox_mediation() currently supports binary factor exposures or numeric exposures.", call. = FALSE)
      if (is.character(exposure_raw) && (is.null(x_exp_a0) || is.null(x_exp_a1))) {
        stop(
          "x_exp is stored as character. Please explicitly specify x_exp_a0 (reference level) and x_exp_a1 (contrasted level). ",
          "Observed levels: ", paste(levels(exposure_factor), collapse = ", "),
          ". Alternatively, convert x_exp to a factor with `factor(x, levels = c(reference, contrasted))`.",
          call. = FALSE
        )
      }
      model_df$exposure <- as.integer(exposure_factor == levels(exposure_factor)[2])
      exposure_levels <- levels(exposure_factor)
      map_exposure_level <- function(value, arg_name) {
        if (length(value) != 1 || is.na(value[1])) {
          stop(arg_name, " must be a single non-missing value coded as 0/1 or matching one of the observed exposure levels.", call. = FALSE)
        }
        value <- value[1]
        if (is.numeric(value) && length(value) == 1 && !is.na(value) && value %in% c(0, 1)) return(as.numeric(value))
        value_chr <- as.character(value)
        if (value_chr %in% exposure_levels) return(match(value_chr, exposure_levels) - 1L)
        stop(
          "For binary factor exposures, ", arg_name,
          " must be coded as 0/1 or match one of the observed exposure levels: ",
          paste(exposure_levels, collapse = ", "), ".",
          call. = FALSE
        )
      }
      if (xor(is.null(a0_use), is.null(a1_use))) {
        stop("Please provide both x_exp_a0 and x_exp_a1 together when overriding the default exposure contrast for a binary factor exposure.", call. = FALSE)
      }
      if (is.null(a0_use)) a0_use <- 0
      if (is.null(a1_use)) a1_use <- 1
      a0_use <- map_exposure_level(a0_use, "x_exp_a0")
      a1_use <- map_exposure_level(a1_use, "x_exp_a1")
      exposure_contrast_display <- paste0(exposure_levels[a0_use + 1], " -> ", exposure_levels[a1_use + 1])
    } else {
      model_df$exposure <- suppressWarnings(as.numeric(exposure_raw))
      if (anyNA(model_df$exposure)) stop("Exposure must be numeric or a binary factor.", call. = FALSE)
      unique_exposure <- sort(unique(stats::na.omit(model_df$exposure)))
      if (length(unique_exposure) == 2 && is.null(a0_use) && is.null(a1_use)) {
        a0_use <- unique_exposure[1]
        a1_use <- unique_exposure[2]
      }
      if (is.null(a0_use) || is.null(a1_use)) stop("Please provide both x_exp_a0 and x_exp_a1 for Cox mediation when x_exp is not a binary factor.", call. = FALSE)
      exposure_contrast_display <- paste0(format_num(a0_use), " -> ", format_num(a1_use))
    }
    if (!is.numeric(a0_use) || length(a0_use) != 1 || is.na(a0_use)) stop("x_exp_a0 must be a single numeric value.", call. = FALSE)
    if (!is.numeric(a1_use) || length(a1_use) != 1 || is.na(a1_use)) stop("x_exp_a1 must be a single numeric value.", call. = FALSE)
    if (isTRUE(all.equal(a0_use, a1_use))) stop("x_exp_a0 and x_exp_a1 must be different for Cox mediation analysis.", call. = FALSE)

    ## mediator processing ----
    mediator_raw <- model_df[[x_med]]
    mediator_mode_use <- mreg
    if (mediator_mode_use == "auto") {
      if (is.factor(mediator_raw) || is.character(mediator_raw) || is.logical(mediator_raw)) {
        mediator_levels <- levels(normalize_factor_like(mediator_raw))
        if (length(mediator_levels) != 2) {
          stop(
            "Automatic mediator detection currently supports binary factors or numeric mediators. ",
            "For multi-level categorical mediators, please recode x_med to a scientifically justified binary mediator, ",
            "or keep an ordered score as numeric and set mreg = 'linear'.",
            call. = FALSE
          )
        }
        mediator_mode_use <- "logistic"
      } else {
        mediator_num <- suppressWarnings(as.numeric(mediator_raw))
        unique_mediator <- sort(unique(stats::na.omit(mediator_num)))
        mediator_is_integer_like <- length(unique_mediator) > 0 && all(abs(unique_mediator - round(unique_mediator)) < 1e-8)
        if (length(unique_mediator) == 2 && mediator_is_integer_like) {
          mediator_mode_use <- "logistic"
        } else {
          if (length(unique_mediator) >= 3 && length(unique_mediator) <= 10) {
            leo.basic::leo_log(
              "Mediator '{x_med}' has only {length(unique_mediator)} unique numeric values ({paste(unique_mediator, collapse = ', ')}). ",
              "In `mreg = 'auto'`, it will be treated as a continuous mediator and fit with `linear`. ",
              "Please confirm this is scientifically intended and ignore this warning if expected.",
              level = "warning",
              verbose = verbose
            )
          }
          mediator_mode_use <- "linear"
        }
      }
    }
    if (mediator_mode_use == "logistic") {
      if (is.numeric(mediator_raw) || is.integer(mediator_raw)) {
        mediator_levels <- as.character(sort(unique(stats::na.omit(as.numeric(mediator_raw)))))
        mediator_factor <- factor(as.character(mediator_raw), levels = mediator_levels)
      } else {
        mediator_factor <- normalize_factor_like(mediator_raw)
      }
      if (nlevels(mediator_factor) != 2) {
        stop(
          "Mediator must have exactly 2 levels when mreg = 'logistic'. ",
          "For ordered score-like mediators such as integer severity indices, keep x_med numeric and use mreg = 'linear'.",
          call. = FALSE
        )
      }
      if (is.character(mediator_raw) && is.null(x_med_cde)) {
        stop(
          "x_med is stored as character. Please explicitly specify x_med_cde (the mediator value for CDE evaluation). ",
          "Observed levels: ", paste(levels(mediator_factor), collapse = ", "),
          ". Alternatively, convert x_med to a factor with `factor(x, levels = ...)`.",
          call. = FALSE
        )
      }
      model_df[[x_med]] <- as.integer(mediator_factor == levels(mediator_factor)[2])
      mediator_levels <- levels(mediator_factor)
      mediator_level_output <- function(level_chr) {
        if (is.numeric(mediator_raw) || is.integer(mediator_raw)) return(as.numeric(level_chr))
        if (is.logical(mediator_raw)) return(level_chr %in% c("TRUE", "T", "true", "1"))
        level_chr
      }
      map_mediator_level <- function(value) {
        if (length(value) != 1 || is.na(value[1])) {
          stop("x_med_cde must be a single non-missing value matching one of the observed mediator levels.", call. = FALSE)
        }
        value <- value[1]
        value_chr <- as.character(value)
        if (value_chr %in% mediator_levels) {
          level_index <- match(value_chr, mediator_levels)
          return(list(
            internal = unname(level_index - 1L),
            output = mediator_level_output(mediator_levels[level_index])
          ))
        }
        if (identical(mediator_levels, c("0", "1")) && is.numeric(value) && !is.na(value) && value %in% c(0, 1)) {
          return(list(
            internal = unname(as.numeric(value)),
            output = unname(as.numeric(value))
          ))
        }
        stop(
          "x_med_cde must match one of the observed mediator levels when mreg = 'logistic': ",
          paste(mediator_levels, collapse = ", "), ".",
          call. = FALSE
        )
      }
      if (is.null(m_cde_use)) m_cde_use <- mediator_level_output(mediator_levels[2])
      mediator_mapping <- map_mediator_level(m_cde_use)
      m_cde_internal_use <- mediator_mapping$internal
      m_cde_output_use <- mediator_mapping$output
      mediator_reference_display <- mediator_levels[m_cde_internal_use + 1]
    } else {
      if (is.factor(mediator_raw) || is.character(mediator_raw) || is.logical(mediator_raw)) {
        stop(
          "Mediator must be truly numeric when mreg = 'linear'. ",
          "If x_med is binary, use mreg = 'auto' or 'logistic'; if it is a multi-level categorical mediator, recode it to a justified binary mediator before using the current function.",
          call. = FALSE
        )
      }
      model_df[[x_med]] <- suppressWarnings(as.numeric(mediator_raw))
      if (anyNA(model_df[[x_med]])) stop("Mediator must be numeric when mreg = 'linear'.", call. = FALSE)
      if (is.null(m_cde_use)) m_cde_use <- stats::median(model_df[[x_med]], na.rm = TRUE)
      if (!is.numeric(m_cde_use) || length(m_cde_use) != 1 || is.na(m_cde_use)) stop("x_med_cde must be a single numeric value.", call. = FALSE)
      m_cde_internal_use <- m_cde_use
      m_cde_output_use <- m_cde_use
      mediator_reference_display <- format_num(m_cde_output_use)
    }

    ## covariate matrix ----
    cov_build <- build_cvar_matrix(model_df[, covariates, drop = FALSE], covariates, x_cov_cond)
    reg_df <- model_df[, c("event", "time", "exposure", x_med), drop = FALSE]
    if (!is.null(cov_build$covariate_df)) reg_df <- cbind(reg_df, cov_build$covariate_df)

    ## regmedint fitting ----
    leo.basic::leo_log("Fitting mediation {model_name} with mediator {x_med}{if (length(covariates) > 0) paste0(' adjusted for ', paste(covariates, collapse = ', ')) else ''}.", verbose = verbose)
    fit_results[[model_name]] <- regmedint::regmedint(
      data = reg_df,
      yvar = "time",
      avar = "exposure",
      mvar = x_med,
      cvar = cov_build$cvar_use,
      eventvar = "event",
      a0 = a0_use,
      a1 = a1_use,
      m_cde = m_cde_internal_use,
      c_cond = cov_build$c_cond_use,
      mreg = mediator_mode_use,
      yreg = yreg_use,
      interaction = interaction,
      casecontrol = FALSE,
      na_omit = FALSE
    )

    ## result extraction and pm stability ----
    med_summary <- summary(fit_results[[model_name]], exponentiate = TRUE)
    med_coef <- as.data.frame(coef(med_summary))
    names(med_coef) <- gsub("\\(", ".", gsub("\\)", "", names(med_coef)))
    med_coef$effect <- rownames(med_coef)
    rownames(med_coef) <- NULL
    med_coef$model <- model_name
    med_coef$exposure <- x_exp
    med_coef$mediator <- x_med
    med_coef$outcome <- y_out[1]
    med_coef$n <- nrow(model_df)
    med_coef$case_n <- sum(model_df$event == 1, na.rm = TRUE)
    med_coef$control_n <- sum(model_df$event == 0, na.rm = TRUE)
    med_coef$person_year <- sum(model_df$time, na.rm = TRUE)
    med_coef$non_event_n <- med_coef$control_n
    med_coef$followup_total <- med_coef$person_year
    med_coef$mreg <- mediator_mode_use
    med_coef$yreg <- yreg_use
    med_coef$x_exp_a0 <- a0_use
    med_coef$x_exp_a1 <- a1_use
    med_coef$x_med_cde <- m_cde_output_use
    med_coef$x_med_cde_internal <- m_cde_internal_use
    med_coef$exposure_contrast_label <- exposure_contrast_display
    med_coef$mediator_reference_label <- mediator_reference_display
    med_coef$x_cov_cond <- cov_build$c_cond_label
    te_row <- med_coef$effect == "te"
    te_ci_includes_null <- FALSE
    pm_note <- NA_character_
    if (sum(te_row, na.rm = TRUE) == 1L) {
      te_ci_includes_null <- is.finite(med_coef$lower[te_row]) &&
        is.finite(med_coef$upper[te_row]) &&
        med_coef$lower[te_row] <= 0 &&
        med_coef$upper[te_row] >= 0
      if (te_ci_includes_null) {
        pm_note <- "Total effect 95% CI includes the null value; proportion mediated (pm) may be unstable and should be interpreted cautiously."
        leo.basic::leo_log(
          "Total effect 95% CI includes the null value in {model_name}; proportion mediated (`pm`) may be unstable and should be interpreted cautiously.",
          level = "warning",
          verbose = verbose
        )
      }
    }
    med_coef$pm_unstable <- te_ci_includes_null
    med_coef$pm_note <- pm_note
    result_rows[[model_name]] <- med_coef

    ## evaluation summary ----
    evaluation_rows[[model_name]] <- data.frame(
      model = model_name,
      exposure = x_exp,
      mediator = x_med,
      outcome = y_out[1],
      n = nrow(model_df),
      case_n = sum(model_df$event == 1, na.rm = TRUE),
      non_event_n = sum(model_df$event == 0, na.rm = TRUE),
      followup_total = sum(model_df$time, na.rm = TRUE),
      event_rate = sum(model_df$event == 1, na.rm = TRUE) / nrow(model_df),
      exposure_contrast = exposure_contrast_display,
      exposure_contrast_value = paste0(format_num(a0_use), " -> ", format_num(a1_use)),
      x_exp_a0 = a0_use,
      x_exp_a1 = a1_use,
      mediator_reference = mediator_reference_display,
      x_med_cde = m_cde_output_use,
      x_med_cde_internal = m_cde_internal_use,
      mreg = mediator_mode_use,
      yreg = yreg_use,
      pm_unstable = te_ci_includes_null,
      pm_note = pm_note,
      interaction = interaction,
      covariates = if (length(covariates) > 0) paste(covariates, collapse = ", ") else NA_character_,
      x_cov_cond = cov_build$c_cond_label,
      stringsAsFactors = FALSE
    )
  }

  # result assembly ----
  result_detail <- do.call(rbind, result_rows)
  rownames(result_detail) <- NULL
  effect_labels <- c(
    cde = "Controlled direct effect",
    pnde = "Pure natural direct effect",
    tnie = "Total natural indirect effect",
    tnde = "Total natural direct effect",
    pnie = "Pure natural indirect effect",
    te = "Total effect",
    pm = "Proportion mediated"
  )
  result_detail$effect_label <- unname(effect_labels[result_detail$effect])
  result_detail$effect_scale <- ifelse(
    result_detail$effect == "pm",
    "Proportion",
    ifelse(result_detail$yreg == "survCox", "Hazard ratio", "Time ratio")
  )
  result_detail$exposure_contrast <- ifelse(is.na(result_detail$exposure_contrast_label), paste0(result_detail$x_exp_a0, " -> ", result_detail$x_exp_a1), result_detail$exposure_contrast_label)
  result_detail$mediator_reference <- ifelse(
    result_detail$effect == "cde",
    ifelse(is.na(result_detail$mediator_reference_label), format_num(result_detail$x_med_cde), result_detail$mediator_reference_label),
    NA_character_
  )
  result_detail$exposure_contrast_label <- NULL
  result_detail$mediator_reference_label <- NULL
  result_detail$control_n <- NULL
  result_detail$person_year <- NULL
  result <- result_detail[, c("model", "effect", "effect_label", "effect_scale", "exposure", "mediator", "outcome", "exposure_contrast", "mediator_reference", "n", "case_n", "non_event_n", "followup_total", "est", "lower", "upper", "p", "exp.est", "exp.lower", "exp.upper", "mreg", "yreg"), drop = FALSE]
  result$Estimate <- ifelse(result$effect == "pm", round(result$est, 3), round(result$exp.est, 3))
  result$`95% CI` <- ifelse(
    result$effect == "pm",
    paste0(sprintf("%.3f", round(result$lower, 3)), ", ", sprintf("%.3f", round(result$upper, 3))),
    paste0(sprintf("%.3f", round(result$exp.lower, 3)), ", ", sprintf("%.3f", round(result$exp.upper, 3)))
  )
  result$`P value` <- vapply(result$p, .format_p_value, character(1))
  result <- result[, c("model", "effect", "effect_label", "effect_scale", "exposure", "mediator", "outcome", "exposure_contrast", "mediator_reference", "n", "case_n", "non_event_n", "followup_total", "Estimate", "95% CI", "P value", "mreg", "yreg"), drop = FALSE]
  names(result) <- c("Model", "Effect code", "Effect", "Scale", "Exposure", "Mediator", "Outcome", "Exposure contrast", "Mediator reference", "N", "Case N", "Non-event N", "Total follow-up", "Estimate", "95% CI", "P value", "Mediator model", "Outcome model")
  evaluation <- do.call(rbind, evaluation_rows)
  rownames(evaluation) <- NULL
  out <- structure(list(result = result, result_detail = result_detail, result_tidy = result_detail, evaluation = evaluation, fit = fit_results), class = "leo_cox_mediation")
  leo.basic::leo_log("Cox mediation analysis completed for {x_exp} -> {x_med} -> {y_out[1]} with {length(model_list)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' Paper-style Cox mediation plot
#'
#' `leo_cox_mediation_plot()` turns a `leo_cox_mediation` result into a compact
#' pathway figure for a selected model. The plot keeps all annotation text on a
#' shared vertical center axis, uses equal-sized boxes, and draws only straight
#' arrows between the exposure, mediator, and outcome nodes.
#'
#' @param x A `leo_cox_mediation` object returned by [leo_cox_mediation()].
#' @param model Character scalar giving the model to display. Defaults to the last model in `x$result`.
#' @param exposure_label Character label shown in the left box. Manual `\n` line breaks are respected.
#' @param mediator_label Character label shown in the middle box. Manual `\n` line breaks are respected.
#' @param outcome_label Character label shown in the right box. Manual `\n` line breaks are respected.
#' @param language One of `"en"` or `"zh"`.
#' @param palette One of `"jama"`, `"jco"`, `"lancet"`, or `"nejm"`.
#' @param add_note Logical; whether to show tutorial-style note text such as the
#'   exposure contrast, CDE mediator reference, and a TE-null caution when
#'   applicable.
#' @param font_family Optional graphics font family. Use generic families such as
#'   `"sans"`, `"serif"`, or `"mono"` for the most robust cross-device output.
#'
#' @return For `leo_cox_mediation()`, a `leo_cox_mediation` object containing a display table in `$result`, a detailed mediation table in `$result_detail`, an evaluation summary in `$evaluation`, and fitted `regmedint` objects in `$fit`. For `leo_cox_mediation_plot()`, an invisible `recordedplot` object after drawing the pathway figure on the active graphics device.
#' @export
#' @rdname leo_cox_mediation
#' @examples
#' if (requireNamespace("regmedint", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 200
#'   age <- rnorm(n, 60, 8)
#'   exposure_num <- rbinom(n, 1, 0.5)
#'   mediator_num <- rbinom(n, 1, plogis(-0.4 + 0.9 * exposure_num + 0.02 * (age - 60)))
#'   time_event <- rexp(n, rate = exp(-5.8 + 0.55 * exposure_num + 0.45 * mediator_num + 0.02 * (age - 60)))
#'   time_censor <- rexp(n, rate = 0.08)
#'   med_df <- data.frame(
#'     outcome = as.integer(time_event <= time_censor),
#'     outcome_censor = pmax(pmin(time_event, time_censor), 0.1),
#'     exposure = factor(exposure_num, levels = c(0, 1), labels = c("Low risk", "High risk")),
#'     mediator = factor(mediator_num, levels = c(0, 1), labels = c("Low inflammation", "High inflammation")),
#'     age = age
#'   )
#'   res_med <- leo_cox_mediation(
#'     df = med_df, y_out = c("outcome", "outcome_censor"),
#'     x_exp = "exposure", x_med = "mediator", x_cov = "age",
#'     verbose = FALSE
#'   )
#'   # with default labels and note
#'   leo_cox_mediation_plot(
#'     res_med, model = "model_1",
#'     exposure_label = "Metabolic\nrisk",
#'     mediator_label = "Inflammation",
#'     outcome_label = "Incident\noutcome"
#'   )
#'   # with no note
#'   leo_cox_mediation_plot(
#'     res_med, model = "model_1",
#'     exposure_label = "Metabolic\nrisk",
#'     mediator_label = "Inflammation",
#'     outcome_label = "Incident\noutcome",
#'     add_note = FALSE
#'   )
#' }
leo_cox_mediation_plot <- function(x, model = NULL, exposure_label = "Exposure", mediator_label = "Mediator", outcome_label = "Outcome",
                                   language = c("en", "zh"), palette = c("jama", "jco", "lancet", "nejm"), add_note = TRUE, font_family = NULL) {
  if (!inherits(x, "leo_cox_mediation")) stop("x must be a `leo_cox_mediation` object returned by leo_cox_mediation().", call. = FALSE)
  if (!is.data.frame(x$result)) stop("x$result is missing or malformed.", call. = FALSE)
  language <- match.arg(language)
  palette <- match.arg(palette)
  if (!is.logical(add_note) || length(add_note) != 1L || is.na(add_note)) stop("add_note must be TRUE or FALSE.", call. = FALSE)

  model_names <- unique(x$result$Model)
  if (length(model_names) == 0) stop("No mediation results are available for plotting.", call. = FALSE)
  if (is.null(model)) model <- tail(model_names, 1)
  if (!model %in% model_names) stop("model must match one of the fitted models in x$result.", call. = FALSE)

  plot_df <- x$result[x$result$Model == model, , drop = FALSE]
  needed_effects <- c("te", "tnie", "pm", "pnde")
  if (!all(needed_effects %in% plot_df$`Effect code`)) stop("The selected model must contain te, tnie, pm, and pnde rows.", call. = FALSE)
  te_row <- plot_df[plot_df$`Effect code` == "te", , drop = FALSE]
  tnie_row <- plot_df[plot_df$`Effect code` == "tnie", , drop = FALSE]
  pm_row <- plot_df[plot_df$`Effect code` == "pm", , drop = FALSE]
  pnde_row <- plot_df[plot_df$`Effect code` == "pnde", , drop = FALSE]
  pm_ci <- suppressWarnings(as.numeric(trimws(strsplit(pm_row$`95% CI`, ",", fixed = TRUE)[[1]])))
  if (length(pm_ci) != 2 || anyNA(pm_ci)) pm_ci <- c(NA_real_, NA_real_)
  evaluation_df <- x$evaluation
  evaluation_row <- if (is.data.frame(evaluation_df)) evaluation_df[evaluation_df$model == model, , drop = FALSE] else NULL
  pm_unstable_flag <- is.data.frame(evaluation_row) && nrow(evaluation_row) == 1L && isTRUE(evaluation_row$pm_unstable[1])
  scale_tag <- if (identical(te_row$Scale, "Hazard ratio")) "HR" else if (identical(te_row$Scale, "Time ratio")) "TR" else "Ratio"
  exposure_contrast_label <- as.character(te_row$`Exposure contrast`[1])
  mediator_reference_label <- if (is.data.frame(evaluation_row) && nrow(evaluation_row) == 1L) {
    as.character(evaluation_row$mediator_reference[1])
  } else {
    as.character(stats::na.omit(plot_df$`Mediator reference`)[1])
  }

  if (is.null(font_family) || !nzchar(font_family)) font_family <- "sans"
  device_font_family <- if (font_family %in% c("sans", "serif", "mono")) font_family else ""

  palette_def <- switch(
    palette,
    jama = list(fill = c(exposure = "#374E55FF", mediator = "#79AF97FF", outcome = "#B24745FF"), text = c(exposure = "white", mediator = "white", outcome = "white")),
    jco = list(fill = c(exposure = "#0073C2FF", mediator = "#EFC000FF", outcome = "#CD534CFF"), text = c(exposure = "white", mediator = "#1F1F1F", outcome = "white")),
    lancet = list(fill = c(exposure = "#00468BFF", mediator = "#42B540FF", outcome = "#ED0000FF"), text = c(exposure = "white", mediator = "white", outcome = "white")),
    nejm = list(fill = c(exposure = "#0072B5FF", mediator = "#E18727FF", outcome = "#BC3C29FF"), text = c(exposure = "white", mediator = "#1F1F1F", outcome = "white"))
  )

  if (language == "zh") {
    top_label <- sprintf("间接效应 (TNIE, %s): %.3f", scale_tag, tnie_row$Estimate)
    pm_label <- if (all(is.finite(pm_ci))) sprintf("中介比例: %.1f%%\n(95%% CI %.1f%%, %.1f%%)", 100 * pm_row$Estimate, 100 * pm_ci[1], 100 * pm_ci[2]) else sprintf("中介比例: %.1f%%", 100 * pm_row$Estimate)
    if (add_note && pm_unstable_flag) pm_label <- paste0(pm_label, "\n总效应 CI 跨空值，需谨慎解释")
    direct_label <- sprintf("直接效应 (PNDE, %s): %.3f", scale_tag, pnde_row$Estimate)
    total_label <- sprintf("总效应 (TE, %s): %.3f", scale_tag, te_row$Estimate)
    context_label <- if (add_note) sprintf("暴露对比: %s\nCDE 中介参考值: %s", exposure_contrast_label, mediator_reference_label) else ""
  } else {
    top_label <- sprintf("Indirect effect (TNIE, %s): %.3f", scale_tag, tnie_row$Estimate)
    pm_label <- if (all(is.finite(pm_ci))) sprintf("Proportion mediated: %.1f%%\n(95%% CI %.1f%%, %.1f%%)", 100 * pm_row$Estimate, 100 * pm_ci[1], 100 * pm_ci[2]) else sprintf("Proportion mediated: %.1f%%", 100 * pm_row$Estimate)
    if (add_note && pm_unstable_flag) pm_label <- paste0(pm_label, "\nTE 95% CI includes the null; interpret cautiously")
    direct_label <- sprintf("Direct effect (PNDE, %s): %.3f", scale_tag, pnde_row$Estimate)
    total_label <- sprintf("Total effect (TE, %s): %.3f", scale_tag, te_row$Estimate)
    context_label <- if (add_note) sprintf("Exposure contrast: %s\nCDE mediator reference: %s", exposure_contrast_label, mediator_reference_label) else ""
  }

  label_cex <- if (language == "zh") 1.12 else 1.05
  top_cex <- if (language == "zh") 1.18 else 1.16
  body_cex <- if (language == "zh") 0.98 else 0.96
  pad_x <- 0.42
  pad_y <- 0.24
  edge_col <- "#44546A"
  xlim <- c(0.2, 11.8)
  ylim <- c(1.0, 8.1)
  plot_width <- diff(xlim)
  plot_height <- diff(ylim)
  plot_bottom <- ylim[1]
  plot_top <- ylim[2]

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(bg = "white", family = device_font_family, xpd = NA, mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim)

  label_lines <- strsplit(c(exposure_label, mediator_label, outcome_label), "\n", fixed = TRUE)
  safe_strwidth <- function(text, cex) {
    tryCatch(
      graphics::strwidth(text, cex = cex, units = "user"),
      error = function(e) plot_width * 0.018 * max(nchar(enc2utf8(text), type = "width"), 1) * cex
    )
  }
  safe_strheight <- function(text, cex) {
    tryCatch(
      graphics::strheight(text, cex = cex, units = "user"),
      error = function(e) plot_height * 0.045 * cex
    )
  }
  draw_text_safe <- function(x, y, labels, cex, col = "black", font = 1, family = device_font_family, adj = NULL) {
    text_args <- list(x = x, y = y, labels = labels, cex = cex, col = col, font = font, family = family)
    if (!is.null(adj)) text_args$adj <- adj
    tryCatch(
      suppressWarnings(do.call(graphics::text, text_args)),
      error = function(e) {
        text_args$labels <- iconv(labels, from = "", to = "ASCII//TRANSLIT", sub = "?")
        suppressWarnings(do.call(graphics::text, text_args))
      }
    )
  }
  max_line_width <- max(vapply(label_lines, function(lines) max(vapply(lines, safe_strwidth, numeric(1), cex = label_cex)), numeric(1)))
  max_nline <- max(vapply(label_lines, length, integer(1)))
  box_width <- max_line_width + 2 * pad_x
  box_height <- max_nline * safe_strheight("M", cex = label_cex) * 1.15 + 2 * pad_y
  x_center <- mean(xlim)
  side_x_pad <- max(box_width * 0.35, plot_width * 0.085)
  left_x <- xlim[1] + box_width / 2 + side_x_pad
  mid_x <- x_center
  right_x <- xlim[2] - box_width / 2 - side_x_pad
  side_y <- plot_bottom + plot_height * 0.415
  mid_y <- side_y + box_height * 1.55 + plot_height * 0.03
  lower_gap <- side_y - plot_bottom
  direct_y <- plot_bottom + lower_gap * 0.38
  total_y <- plot_bottom + lower_gap * 0.24
  context_y <- plot_bottom + lower_gap * 0.09

  draw_box <- function(x, y, label, fill, text_col) {
    xleft <- x - box_width / 2
    xright <- x + box_width / 2
    ybottom <- y - box_height / 2
    ytop <- y + box_height / 2
    graphics::rect(xleft, ybottom, xright, ytop, col = fill, border = fill, lwd = 1.8)
    draw_text_safe(x, y, labels = label, cex = label_cex, col = text_col, font = 2, family = device_font_family, adj = c(0.5, 0.5))
    list(x = x, y = y, xleft = xleft, xright = xright, ybottom = ybottom, ytop = ytop)
  }

  graphics::rect(xlim[1], ylim[1], xlim[2], ylim[2], col = "white", border = NA)

  left_box <- draw_box(left_x, side_y, exposure_label, palette_def$fill["exposure"], palette_def$text["exposure"])
  mid_box <- draw_box(mid_x, mid_y, mediator_label, palette_def$fill["mediator"], palette_def$text["mediator"])
  right_box <- draw_box(right_x, side_y, outcome_label, palette_def$fill["outcome"], palette_def$text["outcome"])
  pm_y <- side_y + (mid_box$ybottom - side_y) * 0.55
  top_y <- min(plot_top - plot_height * 0.08, mid_box$ytop + plot_height * 0.14)

  graphics::arrows(left_box$xright, side_y + 0.05, mid_box$xleft, mid_box$ybottom + 0.06, length = 0.085, lwd = 2.1, col = edge_col)
  graphics::arrows(mid_box$xright, mid_box$ybottom + 0.06, right_box$xleft, side_y + 0.05, length = 0.085, lwd = 2.1, col = edge_col)
  graphics::arrows(left_box$xright, side_y, right_box$xleft, side_y, length = 0.085, lwd = 2.1, col = edge_col)

  draw_text_safe(x_center, top_y, top_label, cex = top_cex, font = 2, family = device_font_family)
  draw_text_safe(x_center, pm_y, pm_label, cex = body_cex, family = device_font_family)
  draw_text_safe(x_center, direct_y, direct_label, cex = top_cex, font = 2, family = device_font_family)
  draw_text_safe(x_center, total_y, total_label, cex = body_cex, family = device_font_family)
  if (add_note && nzchar(context_label)) draw_text_safe(x_center, context_y, context_label, cex = body_cex * 0.86, family = device_font_family)

  out <- grDevices::recordPlot()
  return(invisible(out))
}

# Linear ----

#' Linear regression fitting and formatting helpers
#'
#' `r lifecycle::badge('experimental')`
#'
#' `leo_linear()` fits one or more linear regression models for a
#' single exposure and continuous outcome. `leo_linear_format()`
#' converts the returned object into a wide summary table, tidy result table, or
#' `gtsummary` output.
#'
#' @param df Data frame containing the continuous outcome, exposure, and covariates.
#' @param y_out Character scalar giving the continuous outcome column name.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list
#'   of covariate column-name vectors. All models are fitted on the same
#'   complete-case cohort defined by `y_out`, `x_exp`, and all covariates that
#'   appear in `x_cov`. If `x_cov` is a named list, those model names are
#'   preserved in the output columns.
#' @param x_exp_type Exposure type handling for `x_exp`. Use `"auto"` to infer
#'   from the input type, `"continuous"` to force a numeric linear term, or
#'   `"categorical"` to force factor coding.
#' @param verbose Logical; print progress messages.
#'
#' @return A `leo_linear` object containing the default wide table in
#'   `$result`, the underlying tidy rows in `$result_tidy`, model metadata, and
#'   fitted `lm` objects.
#' @export
#' @examples
#' linear_df <- dplyr::transmute(
#'   mtcars,
#'   outcome = mpg,
#'   wt = wt,
#'   am = factor(am, levels = c(0, 1), labels = c("Auto", "Manual")),
#'   cyl_group = factor(cyl, levels = c(4, 6, 8), labels = c("Cyl4", "Cyl6", "Cyl8"))
#' )
#'
#' model_cont <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("am"),
#'   "Model B" = c("am", "cyl_group")
#' )
#'
#' res_wt <- leo_linear(
#'   df = linear_df, y_out = "outcome", x_exp = "wt",
#'   x_cov = model_cont, verbose = FALSE
#' )
#' res_wt$result
#' leo_linear_format(res_wt, style = "tidy")
#'
#' model_cat <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("am")
#' )
#'
#' res_cyl <- leo_linear(
#'   df = linear_df, y_out = "outcome", x_exp = "cyl_group",
#'   x_cov = model_cat, x_exp_type = "categorical", verbose = FALSE
#' )
#' res_cyl$result
#'
#' if (requireNamespace("gtsummary", quietly = TRUE) && requireNamespace("broom.helpers", quietly = TRUE)) {
#'   leo_linear_format(res_wt, style = "gtsummary")
#' }
leo_linear <- function(df, y_out, x_exp, x_cov = NULL, x_exp_type = "auto", verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Linear Regression", col = "blue")
  if (!is.character(y_out) || length(y_out) != 1) stop("y_out must be a single column name.", call. = FALSE)
  if (!x_exp_type %in% c("auto", "continuous", "categorical")) stop("x_exp_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)

  prep <- .prepare_regression_data(df = df, y_cols = y_out, x_exp = x_exp, x_cov = x_cov)
  if (prep$n_total == 0) stop("df contains no rows.", call. = FALSE)
  if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering.", call. = FALSE)
  if (!is.numeric(prep$data[[y_out]]) && !is.integer(prep$data[[y_out]])) stop("y_out must be numeric for linear regression.", call. = FALSE)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("Complete-case filtering removed {prep$n_removed_complete_case} row(s) with missing values across exposure/outcome/model covariates{if (length(prep$missing_vars) > 0) paste0(': ', paste(prep$missing_vars, collapse = ', ')) else ''}.", level = "warning", verbose = verbose)
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before linear analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_linear(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }
  exposure_type <- .check_var_type(prep$data$exposure, var_name = prep$exposure_name, var_type = x_exp_type, verbose = verbose)
  leo.basic::leo_log("Linear data prepared: total = {prep$n_total}, after complete-case filter = {prep$n_after_complete_case}.", verbose = verbose)

  model_results <- list()
  fit_results <- list()
  for (model_name in names(prep$models)) {
    covariates <- if (is.null(prep$models[[model_name]])) character(0) else prep$models[[model_name]]
    leo.basic::leo_log("Fitting {model_name} with {if (length(covariates) == 0) 'no covariates (crude model)' else paste(covariates, collapse = ', ')}.", verbose = verbose)

    model_df <- prep$data[, c(y_out, "exposure", covariates), drop = FALSE]
    if (exposure_type == "categorical") {
      model_df$exposure <- droplevels(factor(model_df$exposure))
      if (nlevels(model_df$exposure) < 2) stop("x_exp must contain at least 2 exposure levels after filtering.", call. = FALSE)
    } else {
      if (is.character(model_df$exposure) || is.factor(model_df$exposure) || is.logical(model_df$exposure)) {
        stop("x_exp_type = 'continuous' requires a numeric or integer exposure.", call. = FALSE)
      }
      model_df$exposure <- suppressWarnings(as.numeric(model_df$exposure))
    }

    rhs <- "exposure"
    if (length(covariates) > 0) rhs <- c(rhs, paste0("`", covariates, "`"))
    formula <- stats::as.formula(paste("`", y_out, "` ~ ", paste(rhs, collapse = " + "), sep = ""))
    fit <- stats::lm(formula = formula, data = model_df)
    tidy_fit <- broom::tidy(fit, conf.int = TRUE)
    total_n <- nrow(model_df)
    formula_txt <- Reduce(paste, deparse(formula))

    if (!is.factor(model_df$exposure)) {
      exposure_row <- tidy_fit[tidy_fit$term == "exposure", , drop = FALSE]
      model_results[[model_name]] <- data.frame(
        model = model_name,
        row_id = 1L,
        exposure = prep$exposure_name,
        outcome = y_out,
        level = NA_character_,
        n = total_n,
        beta = exposure_row$estimate[1],
        beta_ci_l = exposure_row$conf.low[1],
        beta_ci_u = exposure_row$conf.high[1],
        p_value = exposure_row$p.value[1],
        exposure_class = "Continuous",
        formula = formula_txt,
        check.names = FALSE
      )
      fit_results[[model_name]] <- fit
      next
    }

    levels_x <- levels(model_df$exposure)
    level_counts <- data.frame(level = levels_x, N = vapply(levels_x, function(level) sum(model_df$exposure == level, na.rm = TRUE), integer(1)))
    result_df <- data.frame(
      model = model_name,
      row_id = seq_along(levels_x),
      exposure = c(paste0(levels_x[1], " (Ref)"), levels_x[-1]),
      outcome = y_out,
      level = levels_x,
      n = level_counts$N[match(levels_x, level_counts$level)],
      beta = c(0, rep(NA_real_, length(levels_x) - 1)),
      beta_ci_l = c(0, rep(NA_real_, length(levels_x) - 1)),
      beta_ci_u = c(0, rep(NA_real_, length(levels_x) - 1)),
      p_value = c(NA_real_, rep(NA_real_, length(levels_x) - 1)),
      exposure_class = if (length(levels_x) == 2) "Binary" else paste0("Categorical (", length(levels_x), " levels)"),
      formula = formula_txt,
      check.names = FALSE
    )
    if (length(levels_x) > 1) {
      exposure_rows <- tidy_fit[tidy_fit$term != "(Intercept)", , drop = FALSE]
      exposure_rows <- exposure_rows[seq_len(length(levels_x) - 1L), , drop = FALSE]
      if (nrow(exposure_rows) != length(levels_x) - 1L) stop("Failed to recover categorical exposure coefficients from the linear model.", call. = FALSE)
      result_df$beta[-1] <- exposure_rows$estimate
      result_df$beta_ci_l[-1] <- exposure_rows$conf.low
      result_df$beta_ci_u[-1] <- exposure_rows$conf.high
      result_df$p_value[-1] <- exposure_rows$p.value
    }
    model_results[[model_name]] <- result_df
    fit_results[[model_name]] <- fit
  }

  result_tidy <- do.call(rbind, model_results)
  rownames(result_tidy) <- NULL
  out <- structure(list(
    result = NULL,
    result_tidy = result_tidy,
    data_info = list(
      n_total = prep$n_total,
      n_after_complete_case = prep$n_after_complete_case
    ),
    fit = fit_results
  ), class = "leo_linear")
  out$result <- leo_linear_format(out, style = "wide")
  leo.basic::leo_log("Linear regression completed for {prep$exposure_name} -> {y_out} with {length(prep$models)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' Format a `leo_linear` result
#'
#' @rdname leo_linear
#' @param x Result returned by `leo_linear()`.
#' @param style One of `"wide"`, `"tidy"`, or `"gtsummary"`.
#' @export
leo_linear_format <- function(x, style = "wide") {
  if (!inherits(x, "leo_linear")) stop("x must be a leo_linear result.", call. = FALSE)
  if (!style %in% c("wide", "tidy", "gtsummary")) stop("style must be one of 'wide', 'tidy', or 'gtsummary'.", call. = FALSE)
  switch(style,
    wide = {
      if (!is.null(x$result)) return(x$result)
      result_tidy <- x$result_tidy
      model_ids <- unique(result_tidy$model)
      result_wide <- result_tidy[result_tidy$model == model_ids[1], c("row_id", "exposure", "outcome", "n", "exposure_class"), drop = FALSE]
      names(result_wide) <- c("row_id", "Exposure", "Outcome", "N", "Class")
      for (model_id in model_ids) {
        model_df <- result_tidy[result_tidy$model == model_id, c("row_id", "beta", "beta_ci_l", "beta_ci_u", "p_value"), drop = FALSE]
        model_df$beta <- round(model_df$beta, 3)
        model_df$ci_95 <- ifelse(is.na(model_df$beta_ci_l) | is.na(model_df$beta_ci_u), "NA", paste0(sprintf("%.3f", round(model_df$beta_ci_l, 3)), ", ", sprintf("%.3f", round(model_df$beta_ci_u, 3))))
        model_df$p_value <- vapply(model_df$p_value, .format_p_value, character(1))
        model_df <- model_df[, c("row_id", "beta", "ci_95", "p_value"), drop = FALSE]
        names(model_df) <- c("row_id", paste(model_id, "Beta"), paste(model_id, "95% CI"), paste(model_id, "P value"))
        result_wide <- merge(result_wide, model_df, by = "row_id", all.x = TRUE, sort = FALSE)
        result_wide <- result_wide[order(result_wide$row_id), , drop = FALSE]
      }
      result_wide$row_id <- NULL
      rownames(result_wide) <- NULL
      result_wide
    },
    tidy = {
      result_tidy <- x$result_tidy
      result_out <- result_tidy[, c("model", "exposure", "outcome", "level", "n"), drop = FALSE]
      result_out$Beta <- round(result_tidy$beta, 3)
      result_out$`95% CI` <- ifelse(is.na(result_tidy$beta_ci_l) | is.na(result_tidy$beta_ci_u), "NA", paste0(sprintf("%.3f", round(result_tidy$beta_ci_l, 3)), ", ", sprintf("%.3f", round(result_tidy$beta_ci_u, 3))))
      result_out$`P value` <- vapply(result_tidy$p_value, .format_p_value, character(1))
      result_out$Class <- result_tidy$exposure_class
      result_out$Formula <- result_tidy$formula
      names(result_out)[1:5] <- c("Model", "Exposure", "Outcome", "Level", "N")
      rownames(result_out) <- NULL
      result_out
    },
    gtsummary = {
      if (is.null(x$fit) || length(x$fit) == 0) stop("x does not contain fitted models required for style = 'gtsummary'.", call. = FALSE)
      if (!requireNamespace("gtsummary", quietly = TRUE) || !requireNamespace("broom.helpers", quietly = TRUE)) {
        stop("style = 'gtsummary' requires packages 'gtsummary' and 'broom.helpers'.", call. = FALSE)
      }
      tbls <- lapply(x$fit, function(fit) gtsummary::tbl_regression(fit))
      if (length(tbls) == 1) return(tbls[[1]])
      gtsummary::tbl_merge(tbls = tbls, tab_spanner = paste0("**", names(tbls), "**"))
    }
  )
}

# Logistic ----

#' Logistic regression fitting and formatting helpers
#'
#' `r lifecycle::badge('experimental')`
#'
#' `leo_logistic()` fits one or more logistic regression models for a
#' single exposure and binary outcome. `leo_logistic_format()`
#' converts the returned object into a wide summary table, tidy result table, or
#' `gtsummary` output.
#'
#' @param df Data frame containing the binary outcome, exposure, and covariates.
#' @param y_out Character scalar giving the binary outcome column name.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list
#'   of covariate column-name vectors. All models are fitted on the same
#'   complete-case cohort defined by `y_out`, `x_exp`, and all covariates that
#'   appear in `x_cov`. If `x_cov` is a named list, those model names are
#'   preserved in the output columns.
#' @param case_value Value in `y_out` that indicates case status.
#' @param x_exp_type Exposure type handling for `x_exp`. Use `"auto"` to infer
#'   from the input type, `"continuous"` to force a numeric logistic term, or
#'   `"categorical"` to force factor coding.
#' @param verbose Logical; print progress messages.
#'
#' @return A `leo_logistic` object containing the default wide table
#'   in `$result`, the underlying tidy rows in `$result_tidy`, model metadata,
#'   and fitted `glm` objects.
#' @export
#' @examples
#' set.seed(123)
#' n <- 300
#' logi_df <- data.frame(
#'   outcome = rbinom(n, 1, 0.35),
#'   prs = rnorm(n),
#'   sex = factor(sample(c("Male", "Female"), n, TRUE)),
#'   bmi_group = factor(sample(c("Low", "Mid", "High"), n, TRUE), levels = c("Low", "Mid", "High"))
#' )
#'
#' model_cont <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("sex"),
#'   "Model B" = c("sex", "bmi_group")
#' )
#'
#' res_prs <- leo_logistic(
#'   df = logi_df, y_out = "outcome", x_exp = "prs",
#'   x_cov = model_cont, verbose = FALSE
#' )
#' res_prs$result
#' leo_logistic_format(res_prs, style = "tidy")
#'
#' model_cat <- list(
#'   "Crude" = NULL,
#'   "Model A" = c("sex")
#' )
#'
#' res_bmi <- leo_logistic(
#'   df = logi_df, y_out = "outcome", x_exp = "bmi_group",
#'   x_cov = model_cat, x_exp_type = "categorical", verbose = FALSE
#' )
#' res_bmi$result
#'
#' if (requireNamespace("gtsummary", quietly = TRUE) && requireNamespace("broom.helpers", quietly = TRUE)) {
#'   leo_logistic_format(res_prs, style = "gtsummary")
#' }
leo_logistic <- function(df, y_out, x_exp, x_cov = NULL, case_value = 1, x_exp_type = "auto", verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Logistic Regression", col = "blue")
  if (!is.character(y_out) || length(y_out) != 1) stop("y_out must be a single column name.", call. = FALSE)
  if (!x_exp_type %in% c("auto", "continuous", "categorical")) stop("x_exp_type must be one of 'auto', 'continuous', or 'categorical'.", call. = FALSE)

  prep <- .prepare_regression_data(df = df, y_cols = y_out, x_exp = x_exp, x_cov = x_cov)
  if (prep$n_total == 0) stop("df contains no rows.", call. = FALSE)
  if (nrow(prep$data) == 0) stop("No rows remain after complete-case filtering.", call. = FALSE)
  if (prep$n_removed_complete_case > 0) {
    leo.basic::leo_log("Complete-case filtering removed {prep$n_removed_complete_case} row(s) with missing values across exposure/outcome/model covariates{if (length(prep$missing_vars) > 0) paste0(': ', paste(prep$missing_vars, collapse = ', ')) else ''}.", level = "warning", verbose = verbose)
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before logistic analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_logistic(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }

  outcome_raw <- prep$data[[y_out]]
  outcome_values <- unique(stats::na.omit(as.character(outcome_raw)))
  if (length(outcome_values) != 2) stop("y_out must contain exactly 2 non-missing values for logistic regression.", call. = FALSE)
  outcome <- ifelse(is.na(outcome_raw), NA_integer_, as.integer(as.character(outcome_raw) == as.character(case_value)))
  if (length(unique(stats::na.omit(outcome))) != 2) stop("case_value must match one of the two observed values in y_out.", call. = FALSE)

  logistic_df <- prep$data
  logistic_df$outcome <- outcome
  exposure_type <- .check_var_type(logistic_df$exposure, var_name = prep$exposure_name, var_type = x_exp_type, verbose = verbose)
  leo.basic::leo_log("Logistic data prepared: total = {prep$n_total}, after complete-case filter = {prep$n_after_complete_case}.", verbose = verbose)

  model_results <- list()
  fit_results <- list()
  for (model_name in names(prep$models)) {
    covariates <- if (is.null(prep$models[[model_name]])) character(0) else prep$models[[model_name]]
    leo.basic::leo_log("Fitting {model_name} with {if (length(covariates) == 0) 'no covariates (crude model)' else paste(covariates, collapse = ', ')}.", verbose = verbose)

    model_df <- logistic_df[, c("outcome", "exposure", covariates), drop = FALSE]
    if (exposure_type == "categorical") {
      model_df$exposure <- droplevels(factor(model_df$exposure))
      if (nlevels(model_df$exposure) < 2) stop("x_exp must contain at least 2 exposure levels after filtering.", call. = FALSE)
    } else {
      if (is.character(model_df$exposure) || is.factor(model_df$exposure) || is.logical(model_df$exposure)) {
        stop("x_exp_type = 'continuous' requires a numeric or integer exposure.", call. = FALSE)
      }
      model_df$exposure <- suppressWarnings(as.numeric(model_df$exposure))
    }

    rhs <- "exposure"
    if (length(covariates) > 0) rhs <- c(rhs, paste0("`", covariates, "`"))
    formula <- stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))
    fit <- stats::glm(formula = formula, data = model_df, family = stats::binomial())
    tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    total_case <- sum(model_df$outcome == 1, na.rm = TRUE)
    total_control <- sum(model_df$outcome == 0, na.rm = TRUE)
    formula_txt <- Reduce(paste, deparse(formula))

    if (!is.factor(model_df$exposure)) {
      exposure_row <- tidy_fit[tidy_fit$term == "exposure", , drop = FALSE]
      model_results[[model_name]] <- data.frame(
        model = model_name,
        row_id = 1L,
        exposure = prep$exposure_name,
        outcome = y_out,
        level = NA_character_,
        case_n = total_case,
        control_n = total_control,
        or = exposure_row$estimate[1],
        or_ci_l = exposure_row$conf.low[1],
        or_ci_u = exposure_row$conf.high[1],
        p_value = exposure_row$p.value[1],
        exposure_class = "Continuous",
        formula = formula_txt,
        check.names = FALSE
      )
      fit_results[[model_name]] <- fit
      next
    }

    levels_x <- levels(model_df$exposure)
    level_counts <- do.call(rbind, lapply(levels_x, function(level) {
      level_df <- model_df[model_df$exposure == level, , drop = FALSE]
      data.frame(
        level = level,
        Case_N = sum(level_df$outcome == 1, na.rm = TRUE),
        Control_N = sum(level_df$outcome == 0, na.rm = TRUE)
      )
    }))
    result_df <- data.frame(
      model = model_name,
      row_id = seq_along(levels_x),
      exposure = c(paste0(levels_x[1], " (Ref)"), levels_x[-1]),
      outcome = y_out,
      level = levels_x,
      case_n = level_counts$Case_N[match(levels_x, level_counts$level)],
      control_n = level_counts$Control_N[match(levels_x, level_counts$level)],
      or = c(1, rep(NA_real_, length(levels_x) - 1)),
      or_ci_l = c(1, rep(NA_real_, length(levels_x) - 1)),
      or_ci_u = c(1, rep(NA_real_, length(levels_x) - 1)),
      p_value = c(NA_real_, rep(NA_real_, length(levels_x) - 1)),
      exposure_class = if (length(levels_x) == 2) "Binary" else paste0("Categorical (", length(levels_x), " levels)"),
      formula = formula_txt,
      check.names = FALSE
    )
    if (length(levels_x) > 1) {
      exposure_rows <- tidy_fit[tidy_fit$term != "(Intercept)", , drop = FALSE]
      exposure_rows <- exposure_rows[seq_len(length(levels_x) - 1L), , drop = FALSE]
      if (nrow(exposure_rows) != length(levels_x) - 1L) stop("Failed to recover categorical exposure coefficients from the logistic model.", call. = FALSE)
      result_df$or[-1] <- exposure_rows$estimate
      result_df$or_ci_l[-1] <- exposure_rows$conf.low
      result_df$or_ci_u[-1] <- exposure_rows$conf.high
      result_df$p_value[-1] <- exposure_rows$p.value
    }
    model_results[[model_name]] <- result_df
    fit_results[[model_name]] <- fit
  }

  result_tidy <- do.call(rbind, model_results)
  rownames(result_tidy) <- NULL
  out <- structure(list(
    result = NULL,
    result_tidy = result_tidy,
    data_info = list(
      n_total = prep$n_total,
      n_after_complete_case = prep$n_after_complete_case,
      case_value = case_value
    ),
    fit = fit_results
  ), class = "leo_logistic")
  out$result <- leo_logistic_format(out, style = "wide")
  leo.basic::leo_log("Logistic regression completed for {prep$exposure_name} -> {y_out} with {length(prep$models)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
}

#' Format a `leo_logistic` result
#'
#' @rdname leo_logistic
#' @param x Result returned by `leo_logistic()`.
#' @param style One of `"wide"`, `"tidy"`, or `"gtsummary"`.
#' @export
leo_logistic_format <- function(x, style = "wide") {
  if (!inherits(x, "leo_logistic")) stop("x must be a leo_logistic result.", call. = FALSE)
  if (!style %in% c("wide", "tidy", "gtsummary")) stop("style must be one of 'wide', 'tidy', or 'gtsummary'.", call. = FALSE)
  switch(style,
    wide = {
      if (!is.null(x$result)) return(x$result)
      result_tidy <- x$result_tidy
      model_ids <- unique(result_tidy$model)
      result_wide <- result_tidy[result_tidy$model == model_ids[1], c("row_id", "exposure", "outcome", "case_n", "control_n", "exposure_class"), drop = FALSE]
      names(result_wide) <- c("row_id", "Exposure", "Outcome", "Case N", "Control N", "Class")
      for (model_id in model_ids) {
        model_df <- result_tidy[result_tidy$model == model_id, c("row_id", "or", "or_ci_l", "or_ci_u", "p_value"), drop = FALSE]
        model_df$or <- round(model_df$or, 3)
        model_df$ci_95 <- ifelse(is.na(model_df$or_ci_l) | is.na(model_df$or_ci_u), "NA", paste0(sprintf("%.3f", round(model_df$or_ci_l, 3)), ", ", sprintf("%.3f", round(model_df$or_ci_u, 3))))
        model_df$p_value <- vapply(model_df$p_value, .format_p_value, character(1))
        model_df <- model_df[, c("row_id", "or", "ci_95", "p_value"), drop = FALSE]
        names(model_df) <- c("row_id", paste(model_id, "OR"), paste(model_id, "95% CI"), paste(model_id, "P value"))
        result_wide <- merge(result_wide, model_df, by = "row_id", all.x = TRUE, sort = FALSE)
        result_wide <- result_wide[order(result_wide$row_id), , drop = FALSE]
      }
      result_wide$row_id <- NULL
      rownames(result_wide) <- NULL
      result_wide
    },
    tidy = {
      result_tidy <- x$result_tidy
      result_out <- result_tidy[, c("model", "exposure", "outcome", "level", "case_n", "control_n"), drop = FALSE]
      result_out$OR <- round(result_tidy$or, 3)
      result_out$`95% CI` <- ifelse(is.na(result_tidy$or_ci_l) | is.na(result_tidy$or_ci_u), "NA", paste0(sprintf("%.3f", round(result_tidy$or_ci_l, 3)), ", ", sprintf("%.3f", round(result_tidy$or_ci_u, 3))))
      result_out$`P value` <- vapply(result_tidy$p_value, .format_p_value, character(1))
      result_out$Class <- result_tidy$exposure_class
      result_out$Formula <- result_tidy$formula
      names(result_out)[1:6] <- c("Model", "Exposure", "Outcome", "Level", "Case N", "Control N")
      rownames(result_out) <- NULL
      result_out
    },
    gtsummary = {
      if (is.null(x$fit) || length(x$fit) == 0) stop("x does not contain fitted models required for style = 'gtsummary'.", call. = FALSE)
      if (!requireNamespace("gtsummary", quietly = TRUE) || !requireNamespace("broom.helpers", quietly = TRUE)) {
        stop("style = 'gtsummary' requires packages 'gtsummary' and 'broom.helpers'.", call. = FALSE)
      }
      tbls <- lapply(x$fit, function(fit) gtsummary::tbl_regression(fit, exponentiate = TRUE))
      if (length(tbls) == 1) return(tbls[[1]])
      gtsummary::tbl_merge(tbls = tbls, tab_spanner = paste0("**", names(tbls), "**"))
    }
  )
}
