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
#' @return Formatted p-value string
#' @examples
#' leo.ukb:::.format_p_value(0.2)
#' leo.ukb:::.format_p_value(0.0314)
#' leo.ukb:::.format_p_value(0.0002)
#' @noRd
.format_p_value <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("<0.001")
  if (p_value < 0.01) return(as.character(round(p_value, 3)))
  if (p_value < 0.05) return(as.character(round(p_value, 3)))
  return(as.character(round(p_value, 3)))
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

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  event_raw <- df[[y_out[1]]]
  time_raw <- suppressWarnings(as.numeric(df[[y_out[2]]]))
  event <- ifelse(is.na(event_raw), NA_integer_, as.integer(as.character(event_raw) == as.character(event_value)))

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
#'
#' @return A `leo_cox` object containing the default wide table in
#'   `$result`, the underlying tidy rows in `$result_tidy`, model metadata, and
#'   fitted `coxph` objects.
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
#'   x_exp = "age", x_cov = model_cont, verbose = FALSE
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
#'   x_exp_type = "categorical", verbose = FALSE
#' )
#' res_ecog$result
#'
#' if (requireNamespace("gtsummary", quietly = TRUE) && requireNamespace("broom.helpers", quietly = TRUE)) {
#'   leo_cox_format(res_age, style = "gtsummary")
#' }
leo_cox <- function(df, y_out, x_exp, x_cov = NULL, event_value = 1, min_followup_time = 0, x_exp_type = "auto", verbose = TRUE) {
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
      exposure = c("Ref", levels_x[-1]),
      outcome = prep$outcome_name,
      level = levels_x,
      case_n = level_counts$Case_N[match(levels_x, level_counts$level)],
      control_n = level_counts$Control_N[match(levels_x, level_counts$level)],
      person_year = level_counts$person_year[match(levels_x, level_counts$level)],
      hr = c(1, rep(NA_real_, length(levels_x) - 1)),
      hr_ci_l = c(1, rep(NA_real_, length(levels_x) - 1)),
      hr_ci_u = c(1, rep(NA_real_, length(levels_x) - 1)),
      p_value = c(NA_real_, rep(NA_real_, length(levels_x) - 1)),
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
      event_value = prep$event_value
    ),
    fit = fit_results
  ), class = "leo_cox")
  out$result <- leo_cox_format(out, style = "wide")
  leo.basic::leo_log("Cox regression completed for {prep$exposure_name} -> {prep$outcome_name} with {length(prep$models)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
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
      for (model_id in model_ids) {
        model_df <- result_tidy[result_tidy$model == model_id, c("row_id", "hr", "hr_ci_l", "hr_ci_u", "p_value"), drop = FALSE]
        model_df$hr <- round(model_df$hr, 3)
        model_df$ci_95 <- ifelse(is.na(model_df$hr_ci_l) | is.na(model_df$hr_ci_u), "NA", paste0(sprintf("%.3f", round(model_df$hr_ci_l, 3)), ", ", sprintf("%.3f", round(model_df$hr_ci_u, 3))))
        model_df$p_value <- vapply(model_df$p_value, .format_p_value, character(1))
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
      result_out$HR <- round(result_tidy$hr, 3)
      result_out$`95% CI` <- ifelse(is.na(result_tidy$hr_ci_l) | is.na(result_tidy$hr_ci_u), "NA", paste0(sprintf("%.3f", round(result_tidy$hr_ci_l, 3)), ", ", sprintf("%.3f", round(result_tidy$hr_ci_u, 3))))
      result_out$`P value` <- vapply(result_tidy$p_value, .format_p_value, character(1))
      result_out$Class <- result_tidy$exposure_class
      result_out$Formula <- result_tidy$formula
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
#' binary coding does not use the lowest-incidence joint stratum as the
#' reference group, `leo_cox_add_interaction()` recodes the two binary variables
#' before handing the fitted interaction model to `interactionR`.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, interaction variable, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_inter Character scalar giving the interaction variable column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list of covariate column-name vectors.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#' @param event_value Value in the event column that indicates incident events.
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
leo_cox_add_interaction <- function(df, y_out, x_exp, x_inter, x_cov = NULL, event_value = 1, min_followup_time = 0, verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Additive Interaction", col = "blue")
  # interactionR stays in Suggests, so keep runtime checks and explicit namespace calls.
  if (!requireNamespace("interactionR", quietly = TRUE)) stop("leo_cox_add_interaction() requires package 'interactionR'. Please install it first.", call. = FALSE)
  if (!is.character(x_inter) || length(x_inter) != 1) stop("x_inter must be a single column name.", call. = FALSE)
  if (!x_inter %in% names(df)) stop("x_inter must exist in df.", call. = FALSE)

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

    incidence_df <- dplyr::summarise(
      dplyr::group_by(
        data.frame(event = model_df$event, time = model_df$time, exposure_bin = exposure_bin, interaction_bin = interaction_bin),
        exposure_bin,
        interaction_bin
      ),
      case_n = sum(event),
      person_year = sum(time),
      .groups = "drop"
    )
    incidence_df$incidence_rate <- incidence_df$case_n / incidence_df$person_year
    ref_idx <- which.min(incidence_df$incidence_rate)
    recode_exposure <- incidence_df$exposure_bin[ref_idx] == 1
    recode_interaction <- incidence_df$interaction_bin[ref_idx] == 1
    if (recode_exposure) exposure_bin <- 1L - exposure_bin
    if (recode_interaction) interaction_bin <- 1L - interaction_bin
    reference_group <- paste0(
      x_exp, "=", levels(exposure_factor)[if (recode_exposure) 2 else 1],
      ", ",
      x_inter, "=", levels(interaction_factor)[if (recode_interaction) 2 else 1]
    )

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
    fit_main[[model_name]] <- survival::coxph(formula = formula_main_fit, data = analysis_df)
    fit_inter[[model_name]] <- survival::coxph(formula = formula_inter_fit, data = analysis_df)
    formula_main[[model_name]] <- Reduce(paste, deparse(formula_main_fit))
    formula_inter[[model_name]] <- Reduce(paste, deparse(formula_inter_fit))

    anova_res <- stats::anova(fit_main[[model_name]], fit_inter[[model_name]], test = "Chisq")
    p_col <- grep("^P", names(anova_res), value = TRUE)[1]
    p_interaction <- if (is.na(p_col)) NA_real_ else unname(anova_res[2, p_col])
    backend_obj <- suppressWarnings(interactionR::interactionR(
      model = fit_inter[[model_name]],
      exposure_names = c(x_exp, x_inter),
      ci.type = "delta",
      ci.level = 0.95,
      em = FALSE,
      recode = FALSE
    ))
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
      reference_group = reference_group,
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
          idx <- idx[is.finite(result_df$hr[idx]) & is.finite(result_df$p_value[idx]) & !is.na(result_df$p_value[idx]) & result_df$exposure[idx] != "Ref"]
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
#' `leo_cox_mediation()` performs regression-based causal mediation analysis for
#' a survival outcome using `regmedint::regmedint()` with `yreg = "survCox"`.
#' Because the official `regmedint` interface requires explicit evaluation
#' settings, continuous exposures should usually be accompanied by `a0` and `a1`.
#' The display table in `$result` reports both the standard mediation-effect
#' codes (for example `cde`, `pnde`, and `pm`) and their full labels, together
#' with the exposure contrast and mediator reference value used for evaluation.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, mediator, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_med Character scalar giving the mediator column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list of covariate column-name vectors.
#' @param min_followup_time Numeric scalar; keep rows with `time > min_followup_time`.
#' @param event_value Value in the event column that indicates incident events.
#' @param a0 Numeric scalar; reference exposure value.
#' @param a1 Numeric scalar; contrasted exposure value.
#' @param m_cde Numeric scalar; mediator value at which the controlled direct effect is evaluated.
#' @param c_cond Optional numeric vector of covariate values at which effects are evaluated. If `NULL`, medians are used for numeric covariates.
#' @param mediator_model One of `"auto"`, `"linear"`, or `"logistic"`.
#' @param interaction Logical; whether to include exposure-mediator interaction in the outcome model.
#' @param verbose Logical; print progress messages.
#'
#' @return A `leo_cox_mediation` object containing a display table in `$result`,
#'   the raw mediation table in `$result_tidy`, and fitted `regmedint` objects.
#' @export
#' @examples
#' if (requireNamespace("regmedint", quietly = TRUE)) {
#'   set.seed(123)
#'   n <- 200
#'   med_df <- data.frame(
#'     outcome = rbinom(n, 1, 0.4),
#'     outcome_censor = rexp(n, rate = 0.1),
#'     exposure = rbinom(n, 1, 0.5),
#'     mediator = rnorm(n),
#'     age = rnorm(n, 60, 8)
#'   )
#'   res_med <- leo_cox_mediation(
#'     df = med_df, y_out = c("outcome", "outcome_censor"),
#'     x_exp = "exposure", x_med = "mediator", x_cov = "age",
#'     a0 = 0, a1 = 1, verbose = FALSE
#'   )
#'   res_med$result
#' }
leo_cox_mediation <- function(df, y_out, x_exp, x_med, x_cov = NULL, event_value = 1, min_followup_time = 0, a0 = NULL, a1 = NULL, m_cde = NULL, c_cond = NULL, mediator_model = "auto", interaction = TRUE, verbose = TRUE) {
  t0 <- Sys.time()
  df_name <- deparse(substitute(df))
  if (verbose) cli::cat_rule("Cox Mediation", col = "blue")
  if (!requireNamespace("regmedint", quietly = TRUE)) stop("leo_cox_mediation() requires package 'regmedint'. Please install it first.", call. = FALSE)
  if (!is.character(x_med) || length(x_med) != 1) stop("x_med must be a single column name.", call. = FALSE)
  if (!x_med %in% names(df)) stop("x_med must exist in df.", call. = FALSE)
  mediator_model <- match.arg(mediator_model, c("auto", "linear", "logistic"))

  model_list <- .normalize_model_list(x_cov, df_colnames = names(df))
  prep_models <- lapply(model_list, function(covariates) unique(c(if (is.null(covariates)) character(0) else covariates, x_med)))
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
    leo.basic::leo_log("If you want to keep more rows, consider imputing missing values before Cox mediation analysis, e.g. `{df_name}_imputed <- leo_impute_na({df_name}, method = \"mean\")` or `{df_name}_imputed <- leo_impute_na({df_name}, method = \"rf\")`, then rerun `leo_cox_mediation(df = {df_name}_imputed, ...)`.", level = "warning", verbose = verbose)
  }

  result_rows <- list()
  fit_results <- list()
  for (model_name in names(model_list)) {
    covariates <- if (is.null(model_list[[model_name]])) character(0) else model_list[[model_name]]
    model_df <- prep$data[, c("event", "time", "exposure", x_med, covariates), drop = FALSE]
    a0_use <- a0
    a1_use <- a1
    m_cde_use <- m_cde

    exposure_raw <- model_df$exposure
    if (is.factor(exposure_raw) || is.character(exposure_raw) || is.logical(exposure_raw)) {
      exposure_factor <- droplevels(factor(exposure_raw))
      if (nlevels(exposure_factor) != 2) stop("leo_cox_mediation() currently supports binary factor exposures or numeric exposures.", call. = FALSE)
      model_df$exposure <- as.integer(exposure_factor == levels(exposure_factor)[2])
      if (is.null(a0_use)) a0_use <- 0
      if (is.null(a1_use)) a1_use <- 1
    } else {
      model_df$exposure <- suppressWarnings(as.numeric(exposure_raw))
      unique_exposure <- sort(unique(stats::na.omit(model_df$exposure)))
      if (length(unique_exposure) == 2 && is.null(a0_use) && is.null(a1_use)) {
        a0_use <- unique_exposure[1]
        a1_use <- unique_exposure[2]
      }
      if (is.null(a0_use) || is.null(a1_use)) stop("Please provide both a0 and a1 for Cox mediation when x_exp is not a binary factor.", call. = FALSE)
    }

    mediator_raw <- model_df[[x_med]]
    mediator_mode_use <- mediator_model
    if (mediator_mode_use == "auto") {
      if (is.factor(mediator_raw) || is.character(mediator_raw) || is.logical(mediator_raw)) {
        mediator_levels <- levels(droplevels(factor(mediator_raw)))
        if (length(mediator_levels) != 2) stop("Automatic mediator detection only supports binary factors or numeric mediators.", call. = FALSE)
        mediator_mode_use <- "logistic"
      } else {
        mediator_num <- suppressWarnings(as.numeric(mediator_raw))
        unique_mediator <- sort(unique(stats::na.omit(mediator_num)))
        mediator_mode_use <- if (length(unique_mediator) == 2 && all(abs(unique_mediator - round(unique_mediator)) < 1e-8)) "logistic" else "linear"
      }
    }
    if (mediator_mode_use == "logistic") {
      mediator_factor <- droplevels(factor(mediator_raw))
      if (nlevels(mediator_factor) != 2) stop("Mediator must have exactly 2 levels when mediator_model = 'logistic'.", call. = FALSE)
      model_df[[x_med]] <- as.integer(mediator_factor == levels(mediator_factor)[2])
      if (is.null(m_cde_use)) m_cde_use <- 1
    } else {
      model_df[[x_med]] <- suppressWarnings(as.numeric(mediator_raw))
      if (anyNA(model_df[[x_med]])) stop("Mediator must be numeric when mediator_model = 'linear'.", call. = FALSE)
      if (is.null(m_cde_use)) m_cde_use <- stats::median(model_df[[x_med]], na.rm = TRUE)
    }

    if (length(covariates) > 0) {
      cov_df <- model_df[, covariates, drop = FALSE]
      if (is.null(c_cond)) {
        if (!all(vapply(cov_df, function(x) is.numeric(x) || is.integer(x), logical(1)))) {
          stop("Please supply c_cond explicitly when x_cov contains non-numeric covariates for leo_cox_mediation().", call. = FALSE)
        }
        c_cond_use <- vapply(cov_df, stats::median, numeric(1), na.rm = TRUE)
      } else {
        c_cond_use <- as.numeric(c_cond)
        if (length(c_cond_use) != length(covariates)) stop("c_cond must have the same length as the selected covariates.", call. = FALSE)
      }
    } else {
      c_cond_use <- NULL
    }

    leo.basic::leo_log("Fitting mediation {model_name} with mediator {x_med}{if (length(covariates) > 0) paste0(' adjusted for ', paste(covariates, collapse = ', ')) else ''}.", verbose = verbose)
    fit_results[[model_name]] <- regmedint::regmedint(
      data = model_df,
      yvar = "time",
      avar = "exposure",
      mvar = x_med,
      cvar = if (length(covariates) > 0) covariates else NULL,
      eventvar = "event",
      a0 = a0_use,
      a1 = a1_use,
      m_cde = m_cde_use,
      c_cond = c_cond_use,
      mreg = mediator_mode_use,
      yreg = "survCox",
      interaction = interaction,
      casecontrol = FALSE,
      na_omit = FALSE
    )
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
    med_coef$mediator_model <- mediator_mode_use
    med_coef$a0 <- a0_use
    med_coef$a1 <- a1_use
    med_coef$m_cde <- m_cde_use
    result_rows[[model_name]] <- med_coef
  }

  result_tidy <- do.call(rbind, result_rows)
  rownames(result_tidy) <- NULL
  effect_labels <- c(
    cde = "Controlled direct effect",
    pnde = "Pure natural direct effect",
    tnie = "Total natural indirect effect",
    tnde = "Total natural direct effect",
    pnie = "Pure natural indirect effect",
    te = "Total effect",
    pm = "Proportion mediated"
  )
  result_tidy$effect_label <- unname(effect_labels[result_tidy$effect])
  result_tidy$effect_scale <- ifelse(result_tidy$effect == "pm", "Proportion", "Hazard ratio")
  result_tidy$exposure_contrast <- paste0(result_tidy$a0, " -> ", result_tidy$a1)
  result_tidy$mediator_reference <- ifelse(is.na(result_tidy$m_cde), NA_character_, sprintf("%.3f", round(result_tidy$m_cde, 3)))
  result <- result_tidy[, c("model", "effect", "effect_label", "effect_scale", "exposure", "mediator", "outcome", "exposure_contrast", "mediator_reference", "n", "case_n", "control_n", "person_year", "est", "lower", "upper", "p", "exp.est", "exp.lower", "exp.upper", "mediator_model"), drop = FALSE]
  result$Estimate <- ifelse(result$effect == "pm", round(result$est, 3), round(result$exp.est, 3))
  result$`95% CI` <- ifelse(
    result$effect == "pm",
    paste0(sprintf("%.3f", round(result$lower, 3)), ", ", sprintf("%.3f", round(result$upper, 3))),
    paste0(sprintf("%.3f", round(result$exp.lower, 3)), ", ", sprintf("%.3f", round(result$exp.upper, 3)))
  )
  result$`P value` <- vapply(result$p, .format_p_value, character(1))
  result <- result[, c("model", "effect", "effect_label", "effect_scale", "exposure", "mediator", "outcome", "exposure_contrast", "mediator_reference", "n", "case_n", "control_n", "person_year", "Estimate", "95% CI", "P value", "mediator_model"), drop = FALSE]
  names(result) <- c("Model", "Effect code", "Effect", "Scale", "Exposure", "Mediator", "Outcome", "Exposure contrast", "Mediator reference", "N", "Case N", "Control N", "Person-years", "Estimate", "95% CI", "P value", "Mediator model")
  out <- structure(list(result = result, result_tidy = result_tidy, fit = fit_results), class = "leo_cox_mediation")
  leo.basic::leo_log("Cox mediation analysis completed for {x_exp} -> {x_med} -> {y_out[1]} with {length(model_list)} model(s).", level = "success", verbose = verbose)
  if (verbose) leo.basic::leo_time_elapsed(t0)
  return(out)
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
      exposure = c("Ref", levels_x[-1]),
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
      exposure = c("Ref", levels_x[-1]),
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
