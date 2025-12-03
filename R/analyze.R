# Analyses ----

#' @title Calculate P-value for Heterogeneity from Subgroup Summary Statistics
#' @description This function performs a test for heterogeneity (interaction)
#'              based on summary statistics (Hazard Ratios and P-values) from
#'              a subgroup analysis of a Cox model. It uses Cochran's Q test.
#'
#' @param hrs A numeric vector of Hazard Ratios (HR) for each subgroup.
#' @param p_values A numeric vector of the corresponding P-values for each HR.
#' @param subgroup_names An optional character vector of subgroup names for display.
#'
#' @return A list containing the Q statistic, degrees of freedom (df),
#'         the P-value for heterogeneity, and the I-squared statistic.
#'         It also prints a summary to the console.
#'
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
leo_heterogeneity_p <- function(hrs, p_values, subgroup_names = NULL) {

  # --- 1. Input Validation ---
  if (length(hrs) != length(p_values)) {
    stop("Error: The length of 'hrs' and 'p_values' must be the same.")
  }
  if (length(hrs) < 2) {
    stop("Error: At least two subgroups are required to test for heterogeneity.")
  }
  if (any(hrs <= 0) || any(p_values <= 0) || any(p_values > 1)) {
    stop("Error: Invalid input. HRs must be > 0 and P-values must be between 0 and 1.")
  }

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

# #' Deprecated: General looped linear regression: y ~ x (+ covariates)
# #'
# #' @description Loop over multiple exposures (x) to fit linear models against a single outcome (y).
# #' Input data frames must have the first column as individual ID. Covariates are optional.
# #'
# #' @param x Data frame (first col = ID) of exposures; remaining columns are candidate predictors.
# #' @param y Data frame (first col = ID) of outcome; one outcome column is used (see y_col).
# #' @param cov Optional data frame (first col = ID) of covariates; remaining columns are covariates.
# #' @param x_col Character vector of exposure column names in x; default: all columns from the 2nd.
# #' @param y_col Outcome column name in y; default: the 2nd column of y.
# #' @param id_col Shared ID column name in x/y/cov; default "id".
# #' @param scale_x Logical; if TRUE, regress on scale(x); default FALSE.
# #' @param na_action One of "na.omit" or "na.exclude"; default "na.omit".
# #' @param verbose Logical; if FALSE, suppress logs; default TRUE.
# #'
# #' @return A tibble, each row = one regression (exposure vs outcome).
# #' Columns:
# #' \itemize{
# #'   \item \code{x_name}: exposure variable name
# #'   \item \code{y_name}: outcome variable name
# #'   \item \code{n}: number of non-missing samples used
# #'   \item \code{r2}: R-squared of model
# #'   \item \code{adjr2}: adjusted R-squared
# #'   \item \code{beta}: coefficient for exposure (per 1SD if scale_x=TRUE)
# #'   \item \code{se}: standard error of beta
# #'   \item \code{t}: t-statistic for beta
# #'   \item \code{p}: raw p-value for beta
# #'   \item \code{ci_l}: lower 95% CI of beta
# #'   \item \code{ci_u}: upper 95% CI of beta
# #'   \item \code{p_adj}: BH-adjusted p-value across all exposures in this run
# #' }
# #'
# #' @examples
# #' # --- Example 1: multiple exposures, no covariates ---
# #' set.seed(1)
# #' x <- tibble::tibble(id = 1:120, PRS1 = rnorm(120), PRS2 = rnorm(120))
# #' y <- tibble::tibble(id = 1:120, va1y = 0.2 * scale(x$PRS1)[,1] + rnorm(120))
# #' res1 <- leo_linear_regression(x = x, y = y, x_col = c("PRS1","PRS2"),
# #'                               y_col = "va1y", id_col = "id", scale_x = TRUE)
# #' head(res1)
# #'
# #' # --- Example 2: single exposure with covariates ---
# #' set.seed(2)
# #' x2   <- tibble::tibble(id = 1:150, PRS = rnorm(150))
# #' y2   <- tibble::tibble(id = 1:150, va1y = 0.3 * scale(x2$PRS)[,1] + rnorm(150))
# #' cov2 <- tibble::tibble(id = 1:150, age = rnorm(150, 50, 10), sex = sample(c("M","F"), 150, T))
# #' res2 <- leo_linear_regression(x = x2, y = y2, cov = cov2,
# #'                               x_col = "PRS", y_col = "va1y", id_col = "id",
# #'                               scale_x = TRUE, na_action = "na.omit")
# #' head(res2)
# #'
# #' @export
# #' @importFrom dplyr select any_of left_join bind_rows
# #' @importFrom stats lm as.formula na.omit na.exclude nobs
# #' @importFrom broom tidy glance
# #' @importFrom tibble tibble
# #' @importFrom glue glue
# leo_linear_regression_old <- function(x, y, cov = NULL, x_col = NULL, y_col = NULL,
#                                   id_col = "id", scale_x = F, na_action = "na.omit",
#                                   verbose = TRUE) {
#   # checks
#   if (!id_col %in% names(x)) return(leo.basic::leo_log("id_col '{id_col}' not found in x", level = "danger", verbose = verbose))
#   if (!id_col %in% names(y)) return(leo.basic::leo_log("id_col '{id_col}' not found in y", level = "danger", verbose = verbose))
#   if (!is.null(cov) && !id_col %in% names(cov)) return(leo.basic::leo_log("id_col '{id_col}' not found in cov", level = "danger", verbose = verbose))
#   if (is.null(x_col)) x_col <- names(x)[seq_len(ncol(x))[-1L]]
#   if (length(x_col) == 0) return(leo.basic::leo_log("No exposure columns to analyze in x", level = "danger", verbose = verbose))
#   if (is.null(y_col)) { if (ncol(y) < 2) return(leo.basic::leo_log("y has no outcome column", level = "danger", verbose = verbose)); y_col <- names(y)[2L] }
#   if (!y_col %in% names(y)) return(leo.basic::leo_log("y_col '{y_col}' not found in y", level = "danger", verbose = verbose))
#   if (!all(x_col %in% names(x))) return(leo.basic::leo_log("Some x_col not found in x: {paste(setdiff(x_col, names(x)), collapse=', ')}", level = "danger", verbose = verbose))
#   if (!na_action %in% c("na.omit","na.exclude")) return(leo.basic::leo_log("na_action must be 'na.omit' or 'na.exclude'", level = "danger", verbose = verbose))
#
#   # helpers
#   na_fun <- if (na_action == "na.omit") stats::na.omit else stats::na.exclude
#   bt <- function(z) paste0("`", z, "`")
#
#   # prep data
#   x_df <- dplyr::select(x, dplyr::any_of(c(id_col, x_col))); names(x_df)[1] <- "eid"
#   y_df <- dplyr::select(y, dplyr::any_of(c(id_col, y_col))); names(y_df)[1] <- "eid"
#   if (!is.null(cov)) {
#     cov_df <- dplyr::select(cov, dplyr::any_of(c(id_col, setdiff(names(cov), id_col)))); names(cov_df)[1] <- "eid"
#     cov_names <- setdiff(names(cov_df), "eid")
#   } else { cov_df <- NULL; cov_names <- character(0) }
#   base_df <- y_df; if (!is.null(cov_df)) base_df <- dplyr::left_join(base_df, cov_df, by = "eid")
#
#   # loop exposures
#   leo.basic::leo_log("Performing linear regression with settings:\n y = {y_col}\n x = {paste(x_col, collapse=', ')}\n cov = {ifelse(length(cov_names)==0, '(none)', paste(cov_names, collapse=', '))}", level = "info", verbose = verbose)
#   rows <- list()
#   for (xn in x_col) {
#     df <- dplyr::left_join(base_df, x_df[, c("eid", xn), drop = F], by = "eid")
#     rhs_x <- if (scale_x) glue::glue("scale({bt(xn)})") else bt(xn)
#     rhs_cov <- setdiff(names(df), c("eid", y_col, xn)) # drop cov if it is the x
#     rhs <- c(rhs_x, if (length(rhs_cov) > 0) paste0("`", rhs_cov, "`"))
#     form <- stats::as.formula(paste0(bt(y_col), " ~ ", paste(rhs, collapse = " + ")))
#     fit <- stats::lm(form, data = df, na.action = na_fun)
#
#     tt <- broom::tidy(fit, conf.int = TRUE)
#     exp_pat <- if (scale_x) paste0("^scale\\(", gsub("([\\W])","\\\\\\1", xn), "\\)$") else paste0("^", gsub("([\\W])","\\\\\\1", xn), "$")
#     exp_row <- tt[grepl(exp_pat, tt$term), , drop = F]
#     if (nrow(exp_row) == 0) exp_row <- data.frame(estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_)
#     gl <- broom::glance(fit)
#
#     rows[[length(rows) + 1]] <- tibble::tibble(x_name = xn, y_name = y_col, n = stats::nobs(fit),
#                                                r2 = gl$r.squared, adjr2 = gl$adj.r.squared,
#                                                beta = exp_row$estimate[1], se = exp_row$std.error[1],
#                                                t = exp_row$statistic[1], p = exp_row$p.value[1],
#                                                ci_l = exp_row$conf.low[1], ci_u = exp_row$conf.high[1])
#   }
#
#   out <- dplyr::bind_rows(rows); out$p_adj <- p.adjust(out$p, method = "BH")
#   n_min <- min(out$n, na.rm = T); n_max <- max(out$n, na.rm = T)
#   leo.basic::leo_log("Completed. Processed {length(x_col)} exposure(s). N range: {n_min}–{n_max}.", level = "success", verbose = verbose)
#   return(out)
# }

#' Helper function to format p-values according to common publication standards
#'
#' @param p_value Numeric p-value
#' @return Formatted p-value string
format_p_value <- function(p_value) {
  if (is.na(p_value)) return("NA")
  if (p_value < 0.001) return("<0.001")
  if (p_value < 0.01) return(as.character(round(p_value, 3)))
  if (p_value < 0.05) return(as.character(round(p_value, 3)))
  return(as.character(round(p_value, 3)))
}

# Linear regression analysis for PRS vs clinical outcomes
#' Linear regression analysis
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
#'
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
    complete_data <- leo.impute_na(complete_data, method = impute_method, ...)
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
#'
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
    `P-value` = format_p_value(x_effect$p_value)
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


# NA process ----
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
#' leo.impute_na(sample_data, method = "mean")
#'
#' # Random forest imputation
#' leo.impute_na(sample_data, method = "rf")
leo.impute_na <- function(data, method = "mean", ...) {
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

