# Here, we collect clinical variables derived from the clinical information of UK Biobank participants.

#' Generate New Clinical Indicators
#'
#' This function applies specific functions to generate new clinical indicators.
#' It is designed to work with vectors (e.g., columns of a dataframe).
#'
#' @param df A dataframe containing the required columns for calculation.
#' @param types A vector of indicator types, options are "eGDR" or "TyG".
#' @param remove_assist Logical, whether to remove the assisting columns from the dataframe (default TRUE).
#' @param ... Additional parameters passed to the specific indicator calculation functions.
#'
#' @return A dataframe containing the new clinical indicators.
#' @export
#'
#' @importFrom cli cli_alert_info cli_alert_danger cli_alert_success
#' @importFrom dplyr select all_of %>% mutate
#' @examples
#' df <- data.frame(
#'   hba1c = c(48, 55, 60),
#'   waist = c(80, 85, 90),
#'   hypertension = c(1, 0, 1),
#'   triglycerides = c(1.5, 1.8, 2.0),
#'   glucose = c(5.0, 5.5, 6.0)
#' )
#' leo.ukb::leo_new_clinical_indicators(df, types = c("eGDR", "TyG"))
#' @seealso \code{\link{leo_eGDR}}, \code{\link{leo_TyG}} for Insulin Resistance indicators.
leo_new_clinical_indicators <- function(df, types = c("eGDR", "TyG"), type_id = FALSE, remove_assist = TRUE, ...) {
  # Dictionary to map indicator types to required columns and calculation functions
  indicator_dict <- list(
    "eGDR" = list(
      required_cols = c("hba1c", "waist", "hypertension"),
      required_cols_ukb = c("p30750_i0", "p48_i0", NULL),
      calc_function = leo_eGDR
    ),
    "TyG" = list(
      required_cols = c("triglycerides", "glucose"),
      required_cols_ukb = c("p30870_i0", "p30740_i0"),
      calc_function = leo_TyG
    )
    # Add more indicator types here if needed
  )

  # Loop through each type
  results <- df
  for (type in types) {
    cli::cat_rule(paste("Generating new clinical indicator:", type), col = "blue")
    # Retrieve the required columns and function from the dictionary
    indicator_info <- indicator_dict[[type]]
    required_cols <- switch(as.character(type_id),
                            "FALSE" = indicator_info$required_cols,    # use column name
                            "TRUE" = indicator_info$required_cols_ukb) # use ukb field id
    calc_function <- indicator_info$calc_function
    # check
    cli::cli_alert_info(paste(" - Checking if all required column is valid for {type} calculation."))
    missing_cols <- required_cols[!required_cols %in% colnames(df)]
    if (length(missing_cols) > 0) {
      cli::cli_alert_danger(" - Missing: {paste(missing_cols, collapse = ', ')}. Skip {type}"); next
    } else { cli::cli_alert_success(" - Pass for required columns checking.") }

    # Dynamically apply the function to calculate the indicator
    cli::cli_alert_info(" - Calculating the {type} using the following columns: {.val {required_cols}}")
    columns_to_pass <- lapply(required_cols, function(col) df[[col]])
    results <- results %>% mutate(!!type := do.call(calc_function, columns_to_pass))
    if (remove_assist) { results <- results %>% select(-all_of(required_cols)) }
    cli::cli_alert_success(paste(" - {type} calculation completed.")) # Remove assisting columns if requested
  }
  return(results)
}

# Insulin Resistance -----
#' Calculate Estimated Glucose Disposal Rate (eGDR)
#'
#' This function calculates the eGDR using waist circumference, hypertension status, and HbA1c values.
#' It works element-wise for vectors of HbA1c, waist, and hypertension.
#'
#' @param hba1c A numeric vector of HbA1c values (in mmol/mol).
#' @param waist A numeric vector of waist circumferences (in cm).
#' @param hypertension A numeric vector of hypertension status (1 for yes, 0 for no).
#'
#' @return A numeric vector of eGDR values.
#' @export
#'
#' @importFrom cli cli_alert_info
#' @note
#' ### eGDR: Estimated Glucose Disposal Rate:
#' ```
#' # eGDR is calculated using the following formula:
#' eGDR = 21.158 + (-0.09 * waist circumference [cm]) + (−3.407 * hypertension [yes=1/no=0]) + (−0.551 * HbA1c (%))
#' ```
#' In UKB, HbA1c is measured in **mmol/mol (IFCC)**, and needs to be converted to **% (NGSP)**.
#' - **HbA1c (mmol/mol)** (Field ID: 30750)
#'    - Conversion formula: NGSP (%) = [0.09148 * IFCC] + 2.152 [Ref](https://ngsp.org/ifcc.asp)
#' - **Waist circumference** (Field ID: 48)
#' - **hypertension** need to extract from icd10 or first occurrence records
#' @references
#' 1. \url{https://diabetesjournals.org/care/article/36/8/2280/32950/Use-of-the-Estimated-Glucose-Disposal-Rate-as-a} (Diabetes Care, 2013 July)
#' 2. \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC11439291/#Sec2} (Cardiovasc Diabetol, 2024 Sep)
leo_eGDR <- function(hba1c, waist, hypertension, ...) {
  cli::cli_alert_info(" - Converting HbA1c from mmol/mol(IFCC) to %(NGSP) >>> {.code NGSP=[0.09148*IFCC]+2.152} {.url https://ngsp.org/ifcc.asp}")
  hba1c_percent <- 0.09148 * hba1c + 2.152
  cli::cli_alert_info(" - Calculating eGDR using the formula: {.code eGDR = 21.158 + (-0.09 * waist) + (-3.407 * hypertension) + (-0.551 * hba1c_percent)}")
  eGDR <- 21.158 + (-0.09 * waist) + (-3.407 * hypertension) + (-0.551 * hba1c_percent)
  return(eGDR)
}

#' Calculate Triglyceride-Glucose Index (TyG)
#'
#' This function calculates the Triglyceride-Glucose Index (TyG) using `triglyceride` and `glucose` values.
#'
#' @param triglycerides A numeric vector of triglyceride values (in mmol/L).
#' @param glucose A numeric vector of glucose values (in mmol/L).
#'
#' @return A numeric vector of TyG index values.
#' @export
#'
#' @importFrom cli cli_alert_info
#' @note The formula for TyG is:
#' \deqn{TyG = \log \left( \frac{triglycerides \times glucose}{2} \right)}
#' Triglycerides and glucose are converted from mmol/L to mg/dL using the factors:
#' - Triglycerides: \deqn{mg/dL = mmol/L \times 88.5704}
#' - Glucose: \deqn{mg/dL = mmol/L \times 18.0168}
#' -- In UKB, triglycerides (Field ID: 30870) and glucose (Field ID: 30740) are measured in mmol/L.
#' @references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3662289/} (Cardiovasc Diabetol, 2013)
leo_TyG <- function(triglycerides, glucose, ...) {
  cli::cli_alert_info(" - Converting triglycerides and glucose from {.emph \"mmmol/L\"} (UKB default unit) to {.emph \"mg/dL\"}")
  # Convert triglycerides and glucose from mmol/L to mg/dL
  triglycerides_mg_dl <- triglycerides * 88.5704  # Exact conversion factor
  glucose_mg_dl <- glucose * 18.0168  # Exact conversion factor
  # Calculate TyG index element-wise for the vectors
  cli::cli_alert_info(" - Calculating TyG index using the formula: {.code TyG = ln[(triglycerides * glucose)/2]}")
  TyG <- log((triglycerides_mg_dl * glucose_mg_dl) / 2)
  return(TyG)
}

