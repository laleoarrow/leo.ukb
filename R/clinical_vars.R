# Here, we collect clinical variables derived from the clinical information of UK Biobank participants.

#' Generate New Clinical Indicators
#'
#' This function applies specific functions to generate new clinical indicators.
#' It is designed to work with vectors (e.g., columns of a dataframe).
#'
#' @param df A dataframe containing the required columns for calculation.
#' @param types A vector of indicator types, options are "eGDR", "TyG", "PhysicalActivity", "LGI"...
#' @param type_id Logical; if `FALSE`, the function expects semantic column names (e.g. `"age"`, `"BMI"`).
#'   If `TRUE` (default), UKB field-id style names (e.g. `"p21003_i0"`) are used.
#' @param remove_assist Logical, whether to remove the assisting columns (default `TRUE`).
#' @param remove_other Logical, whether to keep only `eid` and the newly
#'   generated indicator columns in the returned dataframe (default `TRUE`).
#'   If `FALSE`, all original columns in `df` are kept.
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
leo_new_clinical_indicators <- function(df, types = c("eGDR", "TyG"), type_id = TRUE, remove_assist = TRUE, remove_other = TRUE, ...) {
  # Dictionary to map indicator types to required columns and calculation functions
  indicator_dict <- list(
    # Basics ----
    "age" = list(required_cols     = c("age"),
                 required_cols_ukb = c("p21003_i0"),  # age at baseline (instance 0)
                 calc_function     = leo_age),
    "BMI" = list(required_cols     = c("BMI"),
                 required_cols_ukb = c("p21001_i0"),  # BMI at baseline (instance 0)
                 calc_function     = leo_BMI),
    "gender" = list(required_cols     = c("gender"),
                    required_cols_ukb = c("p31"),      # sex, single-instance
                    calc_function     = leo_gender),
    "ethnicity" = list(required_cols     = c("ethnicity"),
                       required_cols_ukb = c("p21000_i0"),  # ethnic background (instance 0)
                       calc_function     = leo_ethnicity),
    "ethnicity_finer" = list(required_cols     = c("ethnicity"),
                             required_cols_ukb = c("p21000_i0"),
                             calc_function     = leo_ethnicity_finer),
    "tdi" = list(required_cols     = c("tdi"),
                 required_cols_ukb = c("p22189"),     # Townsend index, recruitment
                 calc_function     = leo_tdi),
    "smoking_status" = list(required_cols     = c("smoking_status"),
                            required_cols_ukb = c("p20116_i0"),  # smoking status (instance 0)
                            calc_function     = leo_smoking_status),
    "drinking_status" = list(required_cols     = c("drinking_status"),
                             required_cols_ukb = c("p20117_i0"),  # alcohol drinker status (instance 0)
                             calc_function     = leo_drinking_status),
    # Inflammation ----
    "LGI" = list(required_cols      = c("crp", "wbc", "platelet", "neutrophil", "lymphocyte"),
                 required_cols_ukb  = c("p30710_i0", "p30000_i0", "p30080_i0", "p30140_i0", "p30120_i0"),
                 calc_function      = leo_LGI),
    # Insulin Resistance ----
    "eGDR" = list(required_cols = c("hba1c", "waist", "hypertension"),
                  required_cols_ukb = c("p30750_i0", "p48_i0", NULL),
                  calc_function = leo_eGDR),
    "TyG" = list(required_cols = c("triglycerides", "glucose"),
                 required_cols_ukb = c("p30870_i0", "p30740_i0"),
                 calc_function = leo_TyG),
    # Life Style Indicators ----
    "PhysicalActivity" = list(required_cols    = c("p884_i0","p894_i0","p904_i0","p914_i0"),
                              required_cols_ukb= c("p884_i0","p894_i0","p904_i0","p914_i0"),
                              calc_function    = leo_physical_activity)
    # Add more indicator types here if needed
  )

  # Loop through each type
  results <- df
  generated_types <- character(0)  # store successfully generated indicator names
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

    # remove assisting columns if requested
    generated_types <- c(generated_types, type)
    if (remove_assist) { results <- results %>% dplyr::select(-all_of(required_cols)) }
    cli::cli_alert_success(paste(" - {type} calculation completed.")) # Remove assisting columns if requested
  }
  if (remove_other) { results <- results %>% dplyr::select(eid, all_of(generated_types)) }
  return(results)
}

# Basics ----
#' Recode age (years)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003
#'
#' @param age Numeric/character vector for age (years) from p21003.
#' @return Numeric vector of age in years.
#' @export
leo_age <- function(age, ...) {
  age_num <- suppressWarnings(as.numeric(age)) # simple numeric conversion
  return(age_num)
}

#' Recode BMI into categories (1:<18.5, 2:18.5-<25, 3:25-<30, 4:>=30)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001
#'
#' @param BMI Numeric/character vector of BMI (kg/m^2) from p21001.
#'
#' @return Integer vector 1-4; NA if invalid or missing.
#' @export
leo_BMI <- function(BMI, ...) {
  # convert to numeric and cut into 4 categories
  bmi_num <- suppressWarnings(as.numeric(BMI))
  bmi_num[bmi_num <= 0] <- NA

  bmi_cat <- dplyr::case_when(
    is.na(bmi_num)                  ~ NA_integer_,
    bmi_num < 18.5                  ~ 1L,
    bmi_num >= 18.5 & bmi_num < 25  ~ 2L,
    bmi_num >= 25   & bmi_num < 30  ~ 3L,
    bmi_num >= 30                   ~ 4L
  )
  cli::cli_alert_info("BMI recoding:\n  1  <18.5\n  2  18.5-<25\n  3  25-<30\n  4  >=30")
  return(bmi_cat)
}

#' Recode gender (1 = Male, 0 = Female)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31
#'
#' @param gender Numeric/character vector, usually from p31.
#'
#' @return Integer vector (1 = Male, 0 = Female, NA for negative/special codes).
#' @export
leo_gender <- function(gender, ...) {
  # keep only 0/1, others -> NA
  gender_num <- suppressWarnings(as.numeric(gender))
  gender_rec <- ifelse(gender_num %in% c(0, 1), as.integer(gender_num), NA_integer_)
  cli::cli_alert_info("Gender recoding:\n  0 = Female\n  1 = Male\n")
  return(gender_rec)
}

#' Recode ethnicity to White vs non-White
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000
#'
#' @param ethnicity Numeric/character vector, usually from p21000.
#'
#' @return Integer vector (1 = White [1,1001,1002,1003],
#'         0 = non-White [all other positive codes],
#'         NA for -1/-3/missing).
#' @export
leo_ethnicity <- function(ethnicity, ...) {
  eth_num <- suppressWarnings(as.numeric(ethnicity))
  white_codes <- c(1, 1001, 1002, 1003)

  eth_rec <- dplyr::case_when(
    eth_num %in% white_codes                ~ 1L,           # White
    eth_num > 0                             ~ 0L,           # Non-White (other positive codes)
    eth_num %in% c(-1, -3) | is.na(eth_num) ~ NA_integer_,  # Unknown / PNSA
    TRUE                                    ~ NA_integer_
  )

  cli::cli_alert_info("Ethnicity recoding:\n  1/1001/1002/1003  -> White (1)\n  other >0 codes    -> non-White (0)\n  -1/-3/NA          -> NA")
  return(eth_rec)
}

#' Recode ethnicity into finer categories (top-level 1–6)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000
#'
#' @param ethnicity Numeric/character vector, usually from p21000.
#'
#' @return Integer vector:
#'   1 = White              (1,1001,1002,1003)
#'   2 = Mixed              (2,2001,2002,2003,2004)
#'   3 = Asian/Asian Brit.  (3,3001,3002,3003,3004)
#'   4 = Black/Black Brit.  (4,4001,4002,4003,4004)
#'   5 = Chinese            (5)
#'   6 = Other ethnic group (6)
#'   NA for -1/-3/missing/other unexpected codes.
#' @export
leo_ethnicity_finer <- function(ethnicity, ...) {
  eth_num <- suppressWarnings(as.numeric(ethnicity))

  white_codes   <- c(1, 1001, 1002, 1003)
  mixed_codes   <- c(2, 2001, 2002, 2003, 2004)
  asian_codes   <- c(3, 3001, 3002, 3003, 3004)
  black_codes   <- c(4, 4001, 4002, 4003, 4004)
  chinese_codes <- 5
  other_codes   <- 6

  eth_rec <- dplyr::case_when(
    eth_num %in% c(-1, -3) | is.na(eth_num) ~ NA_integer_,  # unknown / PNSA
    eth_num %in% white_codes                ~ 1L,           # White
    eth_num %in% mixed_codes                ~ 2L,           # Mixed
    eth_num %in% asian_codes                ~ 3L,           # Asian or Asian British
    eth_num %in% black_codes                ~ 4L,           # Black or Black British
    eth_num %in% chinese_codes              ~ 5L,           # Chinese
    eth_num %in% other_codes                ~ 6L,           # Other ethnic group
    TRUE                                    ~ NA_integer_   # any other unexpected code
  )

  cli::cli_alert_info(
    "Ethnicity (finer) recoding:\n  1  White (1,1001,1002,1003)\n  2  Mixed (2,2001,2002,2003,2004)\n  3  Asian/Asian Brit. (3,3001,3002,3003,3004)\n  4  Black/Black Brit. (4,4001,4002,4003,4004)\n  5  Chinese (5)\n  6  Other ethnic group (6)\n  -1/-3/NA -> NA"
  )

  return(eth_rec)
}

#' Recode Townsend deprivation index (TDI)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22189
#'
#' @param tdi Numeric/character vector, usually from p22189.
#'
#' @return Numeric vector, with UKB special codes -3/-1 recoded to NA.
#' @export
leo_tdi <- function(tdi, ...) {
  # convert to numeric and drop UKB special codes
  tdi_num <- suppressWarnings(as.numeric(tdi))
  tdi_num[tdi_num %in% c(-3, -1)] <- NA
  return(tdi_num)
}

#' Recode smoking status (0 = Never, 1 = Previous, 2 = Current)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116
#'
#' @param smoking_status Numeric/character vector, usually from p20116.
#'
#' @return Integer vector (0/1/2; NA for negative/special codes or others).
#' @export
leo_smoking_status <- function(smoking_status, ...) {
  # keep 0/1/2, recode -3/-1 and other values to NA
  x <- suppressWarnings(as.numeric(smoking_status))
  x[x %in% c(-3, -1)] <- NA
  x[!is.na(x) & !(x %in% c(0, 1, 2))] <- NA
  cli::cli_alert_info("Smoking status recoding:\n  -3  Prefer not to answer\n  0   Never\n  1   Previous\n  2   Current")
  return(as.integer(x))
}

#' Recode drinking status (0 = Never, 1 = Previous, 2 = Current)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20117
#'
#' @param drinking_status Numeric/character vector, usually from p20117.
#'
#' @return Integer vector (0/1/2; NA for negative/special codes or others).
#' @export
leo_drinking_status <- function(drinking_status, ...) {
  # keep 0/1/2, recode -3/-1 and other values to NA
  x <- suppressWarnings(as.numeric(drinking_status))
  x[x %in% c(-3, -1)] <- NA
  x[!is.na(x) & !(x %in% c(0, 1, 2))] <- NA
  cli::cli_alert_info("Drinking status recoding:\n  -3  Prefer not to answer\n  0   Never\n  1   Previous\n  2   Current")
  return(as.integer(x))
}

# Inflammation ----
#' Calculate Low-Grade Inflammation Score (LGI / INFLA-score)
#'
#' This function calculates the low-grade inflammation score (INFLA-score)
#' using four biomarkers:
#'   - high-sensitivity C-reactive protein (CRP),
#'   - white blood cell (WBC) count,
#'   - platelet (PLT) count,
#'   - neutrophil-to-lymphocyte ratio (NLR = neutrophil / lymphocyte).
#'
#' For each biomarker, within-sample deciles (1–10) are mapped to scores:
#'   - deciles 1–4  -> -4, -3, -2, -1
#'   - deciles 5–6  ->  0,  0
#'   - deciles 7–10 ->  1,  2,  3,  4
#'
#' The LGI / INFLA-score is the sum of the four component scores
#' and ranges from -16 to +16. Higher values indicate higher
#' chronic low-grade inflammation.
#'
#' @param crp Numeric vector of hs-CRP values (mg/L).
#' @param wbc Numeric vector of white blood cell counts (×10^9/L).
#' @param plt Numeric vector of platelet counts (×10^9/L).
#' @param neutrophil Numeric vector of neutrophil counts (×10^9/L).
#' @param lymphocyte Numeric vector of lymphocyte counts (×10^9/L).
#' @param ... Reserved for future use.
#'
#' @return Numeric vector of LGI / INFLA scores (-16 to +16).
#' @export
#'
#' @note
#' Suggested UK Biobank fields (instance 0):
#'   - CRP:        Field 30710  (e.g. p30710_i0)
#'   - WBC:        Field 30000  (e.g. p30000_i0)
#'   - Platelets:  Field 30080  (e.g. p30080_i0)
#'   - Neutrophil: Field 30140  (e.g. p30140_i0)
#'   - Lymphocyte: Field 30120  (e.g. p30120_i0)
#'
#' The fourth component uses NLR instead of GrL, taking advantage
#' of differential white cell counts in UK Biobank.
leo_LGI <- function(crp, wbc, plt, neutrophil, lymphocyte, ...) {
  leo.basic::leo_log("Calculating LGI (INFLA-score) from CRP, WBC, PLT and NLR", level = "info")

  # calculate NLR
  nlr <- neutrophil / lymphocyte

  #' Internal helper: score biomarker by deciles for LGI
  #'
  #' @param x Numeric vector of biomarker values.
  #' @return Integer vector in [-4, 4]; NA if input is NA.
  #' @keywords internal
  leo_LGI_decile_score <- function(x) {
    if (all(is.na(x))) return(rep(NA_integer_, length(x)))
    deciles <- stats::quantile(x, probs = seq(0.1, 0.9, by = 0.1), na.rm = T, type = 7) # deciles
    tile <- cut(x, # cut into 10 deciles (1~10)
                breaks = c(-Inf, deciles, Inf),
                labels = 1:10,
                include.lowest = T, right = T)
    tile_int <- as.integer(tile)

    # mapping: decile 1..10 -> score -4..4
    score_map <- c(-4L, -3L, -2L, -1L, 0L, 0L, 1L, 2L, 3L, 4L)
    out <- ifelse(is.na(tile_int), NA_integer_, score_map[tile_int])
    return(out)
  }

  # score each biomarker by deciles
  crp_score <- leo_LGI_decile_score(crp)
  wbc_score <- leo_LGI_decile_score(wbc)
  plt_score <- leo_LGI_decile_score(plt)
  nlr_score <- leo_LGI_decile_score(nlr)

  # sum component scores
  lgi_score <- crp_score + wbc_score + plt_score + nlr_score

  return(lgi_score)
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

# Life style Indicators -----
#' Calculate Physical Activity Indicator (with internal NA recoding)
#'
#' This function first recodes UKB special codes -3/-1 to NA, then determines
#' if participants meet any of the recommended activity thresholds:
#' - ≥150 min moderate/week,
#' - ≥75 min vigorous/week,
#' - ≥5 days moderate/week,
#' - ≥1 day vigorous/week.
#'
#' @note This function is not ready; only adapted it based on code from 医工科研,
#'       which needs further reference validation. (都还没找到参考文献 真是服了)
#'
#' @param p884_i0 Integer/numeric vector; days/week moderate activity (10+ min), codes -3/-1 → NA.
#' @param p894_i0 Integer/numeric vector; minutes/week moderate activity, codes -3/-1 → NA.
#' @param p904_i0 Integer/numeric vector; days/week vigorous activity (10+ min), codes -3/-1 → NA.
#' @param p914_i0 Integer/numeric vector; minutes/week vigorous activity, codes -3/-1 → NA.
#' @return Integer vector (1 = meets guideline, 0 = does not, NA = all missing).
#' @export
#' @importFrom cli cli_alert_info
leo_physical_activity <- function(p884_i0, p894_i0, p904_i0, p914_i0, ...) {
  cli::cli_alert_info(" - Recoding UKB special codes (-3/-1) to NA")

  # helper to coerce to numeric and recode
  recode_na <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    x_num[x_num %in% c(-3, -1)] <- NA
    x_num
  }

  p884_i0 <- recode_na(p884_i0)
  p894_i0 <- recode_na(p894_i0)
  p904_i0 <- recode_na(p904_i0)
  p914_i0 <- recode_na(p914_i0)

  cli::cli_alert_info(" - Evaluating against activity thresholds")
  moderate_days    <- ifelse(!is.na(p884_i0) & p884_i0 >= 5,   1L, 0L)
  moderate_minutes <- ifelse(!is.na(p894_i0) & p894_i0 >= 150, 1L, 0L)
  vigorous_days    <- ifelse(!is.na(p904_i0) & p904_i0 >= 1,   1L, 0L)
  vigorous_minutes <- ifelse(!is.na(p914_i0) & p914_i0 >= 75,  1L, 0L)

  meets_any <- (moderate_days + moderate_minutes + vigorous_days + vigorous_minutes) >= 1

  all_na <- is.na(p884_i0) & is.na(p894_i0) & is.na(p904_i0) & is.na(p914_i0)
  result <- ifelse(all_na, NA_integer_, as.integer(meets_any))

  return(result)
}


