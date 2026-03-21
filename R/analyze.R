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
#' `r lifecycle::badge('stable')`
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
  names(model_list) <- paste0("model_", seq_along(model_list))

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

#' Prepare a single-exposure Cox regression dataset
#'
#' @param df Analysis data frame.
#' @param y_out Character vector of length 2 giving the event and follow-up time
#'   column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov Covariate specification passed to `.normalize_model_list()`;
#'   each entry should be a covariate column name.
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

  leo.basic::leo_log("Cox data prepared: total = {n_total}, after follow-up filter = {n_after_followup}.", verbose = verbose)

  list(
    data = filtered_df,
    models = model_list,
    outcome_name = y_out[1],
    outcome_time = y_out[2],
    exposure_name = x_exp,
    n_total = n_total,
    n_after_followup = n_after_followup,
    min_followup_time = min_followup_time,
    event_value = event_value
  )
}

#' Cox regression fitting and formatting helpers
#'
#' `leo_cox_regression()` fits one or more Cox proportional hazards models for a
#' single exposure and incident outcome. `leo_cox_regression_format()` converts
#' the returned object into a wide summary table, tidy result table, or
#' `gtsummary` output.
#'
#' @param df Data frame containing the outcome, follow-up time, exposure, and covariates.
#' @param y_out Character vector of length 2 giving the event and follow-up time
#'   column names: `c(event, time)`.
#' @param x_exp Character scalar giving the exposure column name.
#' @param x_cov `NULL`, a character vector of covariate column names, or a list
#'   of covariate column-name vectors.
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
#' @return A `leo_cox_regression` object containing the default wide table in
#'   `$result`, the underlying tidy rows in `$result_tidy`, model metadata, and
#'   optional fitted `coxph` objects.
#' @export
#' @examples
#' lung_df <- stats::na.omit(
#'   dplyr::transmute(
#'     survival::lung,
#'     outcome = as.integer(status == 2),
#'     outcome_censor = time,
#'     age = age,
#'     sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
#'     ph.ecog = factor(ph.ecog)
#'   )
#' ); head(lung_df)
#'
#' model <- list(
#'   crude = NULL,
#'   model_a = c("sex"),
#'   model_b = c("sex", "ph.ecog")
#' )
#'
#' res <- leo_cox_regression(
#'   df = lung_df,
#'   y_out = c("outcome", "outcome_censor"),
#'   x_exp = "age",
#'   x_cov = model
#' )
#'
#' leo_cox_regression_format(res, style = "wide")
#' leo_cox_regression_format(res, style = "tidy")
#' if (requireNamespace("gtsummary", quietly = TRUE) && requireNamespace("broom.helpers", quietly = TRUE)) {
#'   leo_cox_regression_format(res, style = "gtsummary")
#' }
leo_cox_regression <- function(df, y_out, x_exp, x_cov = NULL, event_value = 1,
                               min_followup_time = 0, x_exp_type = "auto",
                               verbose = TRUE) {
  t0 <- Sys.time()
  if (verbose) cli::cat_rule("Cox Regression", col = "blue")
  if (! x_exp_type )
  x_exp_type <- match.arg(x_exp_type)

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
  if (nrow(prep$data) == 0) stop("No rows remain after follow-up filtering.", call. = FALSE)
  exposure_type <- .check_var_type(prep$data$exposure, var_name = prep$exposure_name, var_type = x_exp_type, verbose = verbose)

  # Loop over models
  model_results <- list()
  fit_results <- list()
  for (model_name in names(prep$models)) {
    covariates <- if (is.null(prep$models[[model_name]])) character(0) else prep$models[[model_name]]
    leo.basic::leo_log("Fitting {model_name} with {if (length(covariates) == 0) 'no covariates (crude model)' else paste(covariates, collapse = ', ')}.", verbose = verbose)

    # Prepare analysis data
    model_df <- prep$data[, c("event", "time", "exposure", covariates), drop = FALSE]
    model_df <- stats::na.omit(model_df)
    if (nrow(model_df) == 0) stop("No complete-case rows remain for ", model_name, ".", call. = FALSE)

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
        expose = prep$exposure_name,
        outcome = prep$outcome_name,
        level = NA_character_,
        Case_N = total_case,
        Control_N = total_control,
        Person_years = total_years,
        HR = exposure_row$estimate[1],
        HR_CI_l = exposure_row$conf.low[1],
        HR_CI_u = exposure_row$conf.high[1],
        p = exposure_row$p.value[1],
        class = "continue",
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
        Person_years = sum(level_df$time, na.rm = TRUE)
      )
    }))
    result_df <- data.frame(
      model = model_name,
      row_id = seq_along(levels_x),
      expose = c("Ref", levels_x[-1]),
      outcome = prep$outcome_name,
      level = levels_x,
      Case_N = level_counts$Case_N[match(levels_x, level_counts$level)],
      Control_N = level_counts$Control_N[match(levels_x, level_counts$level)],
      Person_years = level_counts$Person_years[match(levels_x, level_counts$level)],
      HR = c(1, rep(NA_real_, length(levels_x) - 1)),
      HR_CI_l = c(1, rep(NA_real_, length(levels_x) - 1)),
      HR_CI_u = c(1, rep(NA_real_, length(levels_x) - 1)),
      p = c(NA_real_, rep(NA_real_, length(levels_x) - 1)),
      class = paste0(length(levels_x), "_class"),
      formula = formula_txt,
      check.names = FALSE
    )
    if (length(levels_x) > 1) {
      for (i in seq_along(levels_x[-1])) {
        tidy_row <- tidy_fit[tidy_fit$term == paste0("exposure", levels_x[-1][i]), , drop = FALSE]
        result_df$HR[i + 1] <- tidy_row$estimate[1]
        result_df$HR_CI_l[i + 1] <- tidy_row$conf.low[1]
        result_df$HR_CI_u[i + 1] <- tidy_row$conf.high[1]
        result_df$p[i + 1] <- tidy_row$p.value[1]
      }
    }
    model_results[[model_name]] <- result_df
    fit_results[[model_name]] <- fit
  }

  # Build tidy result rows
  result_tidy <- do.call(rbind, model_results)

  # Return result object
  out <- structure(list(
    result = NULL,
    result_tidy = result_tidy,
    data_info = list(
      n_total = prep$n_total,
      n_after_followup = prep$n_after_followup,
      min_followup_time = prep$min_followup_time,
      event_value = prep$event_value
    ),
    fit = fit_results
  ), class = "leo_cox_regression")
  out$result <- leo_cox_regression_format(out, style = "wide")
  leo.basic::leo_log("Cox regression completed for {prep$exposure_name} -> {prep$outcome_name} with {length(prep$models)} model(s).", level = "success", verbose = verbose)
  out
}

#' Format a `leo_cox_regression` result
#'
#' @rdname leo_cox_regression
#' @param x Result returned by `leo_cox_regression()`.
#' @param style One of `"wide"`, `"tidy"`, or `"gtsummary"`.
#' @section Formatting output:
#' `leo_cox_regression_format()` returns:
#' - `"wide"`: a wide summary data frame with one set of HR / CI / p-value
#'   columns per model.
#' - `"tidy"`: the row-level tidy result table stored in `x$result_tidy`, with
#'   formatted display columns added.
#' - `"gtsummary"`: a `gtsummary` regression table, or a merged `gtsummary`
#'   table when multiple models are present.
#'
#' @export
leo_cox_regression_format <- function(x, style = "wide") {
  if (!inherits(x, "leo_cox_regression")) stop("x must be a leo_cox_regression result.", call. = FALSE)
  if (!style %in% c("wide", "tidy", "gtsummary")) stop("style must be one of 'wide', 'tidy', or 'gtsummary'.", call. = FALSE)
  switch(style,
    wide = {
      if (!is.null(x$result)) return(x$result)
      result_tidy <- x$result_tidy
      model_ids <- unique(result_tidy$model)
      result_wide <- result_tidy[result_tidy$model == model_ids[1], c("row_id", "expose", "outcome", "Case_N", "Control_N", "Person_years", "class"), drop = FALSE]
      for (model_id in model_ids) {
        model_df <- result_tidy[result_tidy$model == model_id, c("row_id", "HR", "HR_CI_l", "HR_CI_u", "p"), drop = FALSE]
        model_df$HR <- round(model_df$HR, 3)
        model_df$`95%_CI` <- ifelse(is.na(model_df$HR_CI_l) | is.na(model_df$HR_CI_u), "NA", paste0(sprintf("%.3f", round(model_df$HR_CI_l, 3)), ", ", sprintf("%.3f", round(model_df$HR_CI_u, 3))))
        model_df$`P value` <- vapply(model_df$p, .format_p_value, character(1))
        model_df <- model_df[, c("row_id", "HR", "95%_CI", "P value"), drop = FALSE]
        names(model_df) <- c("row_id", paste0(model_id, " HR"), paste0(model_id, " 95%_CI"), paste0(model_id, " P value"))
        if (grepl("^model_[0-9]+$", model_id)) {
          idx <- sub("^model_", "", model_id)
          names(model_df)[2:4] <- paste0("model_", idx, c(" HR", " 95%_CI", " P value"))
        }
        result_wide <- merge(result_wide, model_df, by = "row_id", all.x = TRUE, sort = FALSE)
        result_wide <- result_wide[order(result_wide$row_id), , drop = FALSE]
      }
      result_wide$row_id <- NULL
      result_wide
    },
    tidy = {
      result_tidy <- x$result_tidy
      result_tidy$HR <- round(result_tidy$HR, 3)
      result_tidy$HR_CI_l <- round(result_tidy$HR_CI_l, 3)
      result_tidy$HR_CI_u <- round(result_tidy$HR_CI_u, 3)
      result_tidy$`95%_CI` <- ifelse(is.na(result_tidy$HR_CI_l) | is.na(result_tidy$HR_CI_u), "NA", paste0(sprintf("%.3f", result_tidy$HR_CI_l), ", ", sprintf("%.3f", result_tidy$HR_CI_u)))
      result_tidy$`P value` <- vapply(result_tidy$p, .format_p_value, character(1))
      result_tidy
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


# Analyses ----













# Linear regression analysis for PRS vs clinical outcomes
#' Linear regression analysis
#' `r lifecycle::badge('experimental')`
#'
#' This function performs linear regression analysis between PRS and clinical outcomes
#' with options for handling missing values through complete-case analysis or imputation.
#'
#' @param data A data frame containing the analysis data
#' @param y_var The name of the outcome variable (response)
#' @param x_var The name of the predictor variable (PRS)
#' @param covariates Optional vector of covariate names to adjust for
#' @param impute_na Logical indicating whether to impute missing values (default: FALSE)
#' @param impute_method Imputation method to use if impute_na = TRUE
#' @param ... Additional arguments passed to the imputation function
#'
#' @return A list containing regression results including model, coefficients, and fit statistics
#' @export
#' @importFrom stats lm as.formula coef
#' @examples
#' # Complete case analysis
#' result_cc <- leo_linear_regression(prs_zh8, "va1y", "PRS",
#'                                      covariates = c("age_onset", "gender"))
#'
#' # With mean imputation
#' result_imputed <- leo_linear_regression(prs_zh8, "va1y", "PRS",
#'                                          covariates = c("age_onset", "gender"),
#'                                          impute_na = TRUE, impute_method = "mean")
leo_linear_regression <- function(data, y_var, x_var, covariates = NULL,
                                  impute_na = FALSE, impute_method = "mean", ...) {
  # Build formula
  if (is.null(covariates)) {
    formula <- stats::as.formula(paste(y_var, "~", x_var))
  } else {
    formula <- stats::as.formula(paste(y_var, "~", x_var, "+", paste(covariates, collapse = " + ")))
  }

  # Select relevant columns
  complete_data <- data %>% dplyr::select(dplyr::all_of(c(y_var, x_var, covariates)))

  # Handle missing values
  if (impute_na) {
    complete_data <- leo_impute_na(complete_data, method = impute_method, ...)
    leo_log("Using {impute_method} imputation for {y_var} ~ {x_var}, sample size: {nrow(complete_data)}")
  } else {
    complete_data <- complete_data %>% stats::na.omit()
    leo_log("Using complete-case analysis for {y_var} ~ {x_var}, sample size: {nrow(complete_data)}")
  }

  # Check sample size
  if (nrow(complete_data) < 10) {
    leo_log("Insufficient sample size (<10) for {y_var}", level = "warning")
    return(NULL)
  }

  # Fit linear model
  model <- stats::lm(formula, data = complete_data)
  model_summary <- summary(model)

  # Extract key results
  coefficients <- stats::coef(model_summary)
  x_coef <- coefficients[x_var, ]

  result <- list(
    formula = formula,
    model = model,
    summary = model_summary,
    n = nrow(complete_data),
    r_squared = model_summary$r.squared,
    adj_r_squared = model_summary$adj.r.squared,
    coefficients = coefficients,
    x_effect = data.frame(beta = x_coef["Estimate"], se = x_coef["Std. Error"],
                          t_value = x_coef["t value"], p_value = x_coef["Pr(>|t|)"]),
    imputation = if (impute_na) impute_method else "complete_case"
  )
  if (impute_na) {
    result$imputed_data <- complete_data
  }
  leo_log("Regression completed: beta = {round(x_coef['Estimate'], 3)}, p = {round(x_coef['Pr(>|t|)'], 4)}")
  return(result)
}

#' Extract and format linear regression results for publication
#' `r lifecycle::badge("experimental")`
#'
#' This function extracts key results from linear regression models and formats them
#' in a publication-ready table with standardized effect sizes and confidence intervals.
#'
#' @param regression_result List object returned by leo_linear_regression function
#' @param round_digits Number of digits to round numeric values (default: 3)
#' @param include_ci Whether to include confidence intervals (default: TRUE)
#' @param include_stats Whether to include model statistics (R², sample size) (default: TRUE)
#'
#' @return A data frame with formatted regression results suitable for publication
#' @export
#' @importFrom dplyr tibble
#' @importFrom stats confint
#' @examples
#' # Run regression
#' result <- leo_linear_regression(prs_zh8, "va1y", "PRS", covariates = c("age_onset", "gender"))
#'
#' # Format results for publication
#' pub_table <- leo_regression_result(result)
#' print(pub_table)
leo_regression_result <- function(regression_result, round_digits = 3,
                                  include_ci = TRUE, include_stats = TRUE) {

  if (is.null(regression_result)) {
    leo_log("Regression result is NULL", level = "warning")
    return(NULL)
  }

  # Extract main effect for the predictor variable (x_var)
  x_effect <- regression_result$x_effect
  model <- regression_result$model

  # Calculate confidence intervals
  ci <- tryCatch({
    stats::confint(model, level = 0.95)
  }, error = function(e) {
    matrix(NA, nrow = nrow(regression_result$coefficients), ncol = 2)
  })

  # Get the row index for the main predictor
  pred_names <- rownames(regression_result$coefficients)
  x_index <- which(pred_names == names(regression_result$model$coefficients)[2])[1]  # Second coefficient is usually the main predictor

  if (length(x_index) == 0 || is.na(x_index)) {
    x_index <- 2  # Fallback to second coefficient
  }

  # Create basic result table
  result_df <- dplyr::tibble(
    Predictor = names(regression_result$model$coefficients)[x_index],
    Outcome = as.character(regression_result$formula)[2],
    Beta = round(x_effect$beta, round_digits),
    SE = round(x_effect$se, round_digits),
    `t-value` = round(x_effect$t_value, round_digits),
    `P-value` = .format_p_value(x_effect$p_value)
  )

  # Add confidence intervals if requested
  if (include_ci && nrow(ci) >= x_index) {
    result_df$`CI Lower` <- round(ci[x_index, 1], round_digits)
    result_df$`CI Upper` <- round(ci[x_index, 2], round_digits)
  }

  # Add model statistics if requested
  if (include_stats) {
    result_df$N <- regression_result$n
    result_df$`R-squared` <- round(regression_result$r_squared, round_digits)
    result_df$`Adj. R-squared` <- round(regression_result$adj_r_squared, round_digits)
  }

  # Add imputation method if used
  if (!is.null(regression_result$imputation)) {
    result_df$`Missing Data` <- regression_result$imputation
  }

  # Add formula for reference
  result_df$Formula <- as.character(regression_result$formula)[3]

  leo_log("Regression results formatted for {result_df$Predictor} ~ {result_df$Outcome}")
  return(result_df)
}

#' Batch linear regression analysis for multiple predictors
#'
#' This function performs linear regression analysis for multiple predictor variables
#' against a single outcome, with consistent formatting of results.
#'
#' @param data Data frame containing all variables
#' @param y_var Outcome variable name
#' @param x_vars Vector of predictor variable names to test
#' @param covariates Covariates to adjust for in all models
#' @param ... Additional arguments passed to leo_linear_regression
#'
#' @return A data frame with consolidated results from all regression models
#' @export
#'
#' @examples
#' # Analyze multiple PRS methods against the same outcome
#' x_vars <- c("PRS_plink", "PRS_catboost", "PRS_lasso")
#' results <- leo_linear_regression_loop(prs_zh8, "va1y", x_vars,
#'                                      covariates = c("age_onset", "gender"))
leo_linear_regression_loop <- function(data, y_var, x_vars, covariates = NULL, ...) {
  # Placeholder function - to be implemented
  leo_log("leo_linear_regression_loop function placeholder - to be implemented")
  leo_log("Would analyze {length(x_vars)} predictors against outcome {y_var}")
  return(NULL)
}

#' General looped logistic regression: y (binary) ~ x (+ covariates)
#'
#' @description Loop over multiple exposures (x) to fit logistic regressions against a single binary outcome (y).
#' Input data frames must have the first column as individual ID. Covariates are optional.
#'
#' @param x Data frame (first col = ID) of exposures; remaining columns are candidate predictors (numeric or factor).
#' @param y Data frame (first col = ID) of binary outcome; one column is used (see y_col).
#'          If numeric, it can be 0/1 or any two distinct numeric codes; if factor/character, must have exactly two levels.
#' @param cov Optional data frame (first col = ID) of covariates; remaining columns are covariates.
#' @param x_col Character vector of exposure column names in x; default: all columns from the 2nd.
#' @param y_col Outcome column name in y; default: the 2nd column of y.
#' @param id_col Shared ID column name in x/y/cov; default "id".
#' @param y_level Two labels found in y that define reference/contrast order (first = reference); default c("0","1").
#' @param y_label Descriptive names to print in logs parallel to y_level; default c("0","1").
#' @param scale_x Logical; if TRUE, scale() numeric exposures; default FALSE.
#' @param na_action One of "na.omit" or "na.exclude"; default "na.omit".
#' @param verbose Logical; if FALSE, suppress logs; default TRUE.
#'
#' @return A tibble; each exposure may produce multiple rows when categorical:
#' \itemize{
#'   \item \code{x_name}: exposure variable name
#'   \item \code{x_type}: "numeric" or "factor"
#'   \item \code{x_level}: "(continuous)" for numeric; for factor: level name（参考水平标注为 "(ref)"）
#'   \item \code{ref_level}: reference level（factor 才有；numeric 为 NA）
#'   \item \code{y_name}: outcome variable name
#'   \item \code{n}: number of non-missing samples used
#'   \item \code{r2_mcfadden}: McFadden's pseudo-R² = 1 - (model deviance / null deviance)
#'   \item \code{beta}, \code{se}, \code{z}, \code{p}: coefficient on logit scale (per 1SD if \code{scale_x=TRUE})
#'   \item \code{ci_l}, \code{ci_u}: 95\% CI of \code{beta} (logit scale)
#'   \item \code{OR}, \code{OR_CI_l}, \code{OR_CI_u}: odds ratio and 95\% CI（参考水平 OR=1）
#'   \item \code{p_adj}: BH-adjusted p-value across all rows in this run
#' }
#'
#' @examples
#' # --- Example 1: numeric outcome 0/1, numeric exposure (scaled), no covariates ---
#' set.seed(10)
#' x <- tibble::tibble(id = 1:300, PRS = rnorm(300))
#' linpred <- 0.8 * scale(x$PRS)[,1]
#' y <- tibble::tibble(id = 1:300, outcome = rbinom(300, 1, prob = stats::plogis(linpred)))
#' res1 <- leo_logistic_regression(x = x, y = y, x_col = "PRS", y_col = "outcome",
#'                                 id_col = "id", scale_x = T,
#'                                 y_level = c("0","1"), y_label = c("control","case"))
#' head(res1)
#'
#' # --- Example 2: factor exposure with explicit reference level, plus covariates ---
#' set.seed(11)
#' x2 <- tibble::tibble(id = 1:400,
#'                      group = factor(sample(c("A","B","C"), 400, T), levels = c("A","B","C")))  # A is ref
#' cov2 <- tibble::tibble(id = 1:400, age = rnorm(400, 50, 10))
#' y2  <- tibble::tibble(id = 1:400, outcome = rbinom(400, 1, 0.4))
#' res2 <- leo_logistic_regression(x = x2, y = y2, cov = cov2,
#'                                 x_col = "group", y_col = "outcome",
#'                                 id_col = "id", scale_x = F,
#'                                 y_level = c("0","1"), y_label = c("control","case"))
#' res2
#'
#' @export
#' @importFrom dplyr select any_of left_join bind_rows
#' @importFrom stats glm binomial na.omit na.exclude nobs as.formula
#' @importFrom broom tidy glance
#' @importFrom tibble tibble
#' @importFrom glue glue
leo_logistic_regression <- function(x, y, cov = NULL, x_col = NULL, y_col = NULL, id_col = "id",
                                    y_level = c("0","1"), y_label = c("0","1"),
                                    scale_x = F, na_action = "na.omit",
                                    verbose = TRUE) {
  # checks
  if (!id_col %in% names(x)) return(leo.basic::leo_log("id_col '{id_col}' not found in x", level = "danger", verbose = verbose))
  if (!id_col %in% names(y)) return(leo.basic::leo_log("id_col '{id_col}' not found in y", level = "danger", verbose = verbose))
  if (!is.null(cov) && !id_col %in% names(cov)) return(leo.basic::leo_log("id_col '{id_col}' not found in cov", level = "danger", verbose = verbose))
  if (is.null(x_col)) x_col <- names(x)[seq_len(ncol(x))[-1L]]
  if (length(x_col) == 0) return(leo.basic::leo_log("No exposure columns to analyze in x", level = "danger", verbose = verbose))
  if (is.null(y_col)) { if (ncol(y) < 2) return(leo.basic::leo_log("y has no outcome column", level = "danger", verbose = verbose)); y_col <- names(y)[2L] }
  if (!y_col %in% names(y)) return(leo.basic::leo_log("y_col '{y_col}' not found in y", level = "danger", verbose = verbose))
  if (!all(x_col %in% names(x))) return(leo.basic::leo_log("Some x_col not found in x: {paste(setdiff(x_col, names(x)), collapse=', ')}", level = "danger", verbose = verbose))
  if (!na_action %in% c("na.omit","na.exclude")) return(leo.basic::leo_log("na_action must be 'na.omit' or 'na.exclude'", level = "danger", verbose = verbose))

  # helpers
  na_fun <- if (na_action == "na.omit") stats::na.omit else stats::na.exclude
  bt <- function(z) paste0("`", z, "`")

  # prep data
  x_df <- dplyr::select(x, dplyr::any_of(c(id_col, x_col))); names(x_df)[1] <- "eid"
  y_df <- dplyr::select(y, dplyr::any_of(c(id_col, y_col))); names(y_df)[1] <- "eid"

  # y must be binary (numeric with 2 codes, or factor/character with 2 levels)
  y_vec <- y_df[[2]]
  if (is.numeric(y_vec)) {
    uniq <- sort(unique(na.omit(y_vec)))
    if (length(uniq) != 2) return(leo.basic::leo_log("Outcome '{y_col}' must have exactly two numeric codes or be a 2-level factor.", level = "danger", verbose = verbose))
    lev_order <- if (all(as.character(y_level) %in% as.character(uniq))) as.character(y_level) else as.character(uniq)
    y_df[[2]] <- factor(as.character(y_vec), levels = lev_order)
  } else {
    levs <- unique(na.omit(as.character(y_vec)))
    if (length(levs) != 2) return(leo.basic::leo_log("Outcome '{y_col}' must have exactly two levels. Current: {paste(levs, collapse=', ')}", level = "danger", verbose = verbose))
    lev_order <- if (all(as.character(y_level) %in% levs)) as.character(y_level) else levs
    y_df[[2]] <- factor(as.character(y_vec), levels = lev_order)
  }

  if (!is.null(cov)) {
    cov_df <- dplyr::select(cov, dplyr::any_of(c(id_col, setdiff(names(cov), id_col)))); names(cov_df)[1] <- "eid"
    cov_names <- setdiff(names(cov_df), "eid")
  } else { cov_df <- NULL; cov_names <- character(0) }

  base_df <- y_df; if (!is.null(cov_df)) base_df <- dplyr::left_join(base_df, cov_df, by = "eid")

  # header log
  leo.basic::leo_log("Performing logistic regression with settings:\n y = {y_col} ({lev_order[1]}={y_label[1]}, {lev_order[2]}={y_label[2]})\n x = {paste(x_col, collapse=', ')}\n cov = {ifelse(length(cov_names)==0, '(none)', paste(cov_names, collapse=', '))}",
                     level = "info", verbose = verbose)

  # loop
  rows <- list()
  for (xn in x_col) {
    df <- dplyr::left_join(base_df, x_df[, c("eid", xn), drop = F], by = "eid")

    # remind for categorical exposure
    if (is.character(df[[xn]]) || is.factor(df[[xn]])) {
      levs_x <- levels(factor(df[[xn]]))
      leo.basic::leo_log("Exposure '{xn}' is categorical; please ensure factor() with first level as reference. Levels: {paste(levs_x, collapse=', ')} (ref = {levs_x[1]}).", level = "warning", verbose = verbose)
    }

    rhs_x <- if (is.numeric(df[[xn]]) && scale_x) glue::glue("scale({bt(xn)})") else bt(xn)
    rhs_cov <- setdiff(names(df), c("eid", y_col, xn))
    rhs <- c(rhs_x, if (length(rhs_cov) > 0) paste0("`", rhs_cov, "`"))
    form <- stats::as.formula(paste0(bt(y_col), " ~ ", paste(rhs, collapse = " + ")))
    fit <- stats::glm(form, data = df, family = stats::binomial(), na.action = na_fun)

    gl <- broom::glance(fit); r2_mcfadden <- 1 - (gl$deviance / gl$null.deviance)
    tt <- broom::tidy(fit, conf.int = T)

    if (is.numeric(df[[xn]])) {
      # numeric exposure: single row
      exp_pat <- if (scale_x) paste0("^scale\\(", gsub("([\\W])","\\\\\\1", xn), "\\)$") else paste0("^", gsub("([\\W])","\\\\\\1", xn), "$")
      exp_row <- tt[grepl(exp_pat, tt$term), , drop = F]
      if (nrow(exp_row) == 0) exp_row <- data.frame(estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
      beta <- exp_row$estimate[1]; se <- exp_row$std.error[1]; z <- exp_row$statistic[1]; p <- exp_row$p.value[1]
      ci_l <- exp_row$conf.low[1]; ci_u <- exp_row$conf.high[1]
      OR <- base::exp(beta); OR_CI_l <- base::exp(ci_l); OR_CI_u <- base::exp(ci_u)
      rows[[length(rows) + 1]] <- tibble::tibble(x_name = xn, x_type = "numeric", x_level = "(continuous)", ref_level = NA_character_,
                                                 y_name = y_col, n = stats::nobs(fit), r2_mcfadden = r2_mcfadden,
                                                 beta = beta, se = se, z = z, p = p,
                                                 ci_l = ci_l, ci_u = ci_u, OR = OR, OR_CI_l = OR_CI_l, OR_CI_u = OR_CI_u)
    } else {
      # factor exposure: include explicit reference row (OR=1), then contrasts
      levs_x <- levels(factor(df[[xn]])); ref <- levs_x[1]
      rows[[length(rows) + 1]] <- tibble::tibble(x_name = xn, x_type = "factor", x_level = "(ref)", ref_level = ref,
                                                 y_name = y_col, n = stats::nobs(fit), r2_mcfadden = r2_mcfadden,
                                                 beta = NA_real_, se = NA_real_, z = NA_real_, p = NA_real_,
                                                 ci_l = NA_real_, ci_u = NA_real_,
                                                 OR = 1, OR_CI_l = 1, OR_CI_u = 1)
      base_pat <- paste0("^`?", gsub("([\\W])","\\\\\\1", xn))
      sub_tt <- tt[grepl(base_pat, tt$term), , drop = F]
      for (i in seq_len(nrow(sub_tt))) {
        trm <- sub_tt$term[i]
        lvl <- gsub("^`|`$", "", gsub(paste0("^`?", gsub("([\\W])","\\\\\\1", xn)), "", trm))
        beta <- sub_tt$estimate[i]; se <- sub_tt$std.error[i]; z <- sub_tt$statistic[i]; p <- sub_tt$p.value[i]
        ci_l <- sub_tt$conf.low[i]; ci_u <- sub_tt$conf.high[i]
        OR <- base::exp(beta); OR_CI_l <- base::exp(ci_l); OR_CI_u <- base::exp(ci_u)
        rows[[length(rows) + 1]] <- tibble::tibble(x_name = xn, x_type = "factor", x_level = lvl, ref_level = ref,
                                                   y_name = y_col, n = stats::nobs(fit), r2_mcfadden = r2_mcfadden,
                                                   beta = beta, se = se, z = z, p = p,
                                                   ci_l = ci_l, ci_u = ci_u, OR = OR, OR_CI_l = OR_CI_l, OR_CI_u = OR_CI_u)
      }
      leo.basic::leo_log("For categorical exposure '{xn}', interpretation is vs reference '{ref}'.", level = "info", verbose = verbose)
    }
  }

  out <- dplyr::bind_rows(rows); out$p_adj <- stats::p.adjust(out$p, method = "BH")
  n_min <- min(out$n, na.rm = T); n_max <- max(out$n, na.rm = T)
  leo.basic::leo_log("Completed. Processed {length(x_col)} exposure(s). N range: {n_min}–{n_max}.", level = "success", verbose = verbose)
  return(out)
}
