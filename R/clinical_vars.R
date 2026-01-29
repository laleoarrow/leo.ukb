# Here, we collect clinical variables derived from the clinical information of UK Biobank participants.
# - Main entry: leo_new_clinical_indicators() acts as a dispatcher.
# - indicator_dict maps each type -> required columns (semantic vs UKB field IDs) + a vectorized leo_*() function.
# - For each type: check required columns -> do.call(calc_function, cols) -> mutate(new column) -> optionally drop assist cols / keep only eid + generated cols.
# - Logging: the main function prints task-level progress; sub-functions print only essential notes (e.g., unit conversion, special-code handling) and usually avoid success messages.
# - Helper recoders (e.g., recode_na31()) are internal utilities used by sub-functions and are not registered in indicator_dict.

#' @title Generate New Clinical Indicators
#'
#' @description
#' This function applies specific functions to generate new clinical indicators.
#' It is designed to work with vectors (e.g., columns of a dataframe).
#'
#' @details
#' Supported indicators (use values in `types`):
#' \itemize{
#'   \item Basics: age, BMI, gender, ethnicity, ethnicity_finer.
#'   \item Socioeconomic: tdi, education, household_income, household_income_2, career.
#'   \item Inflammation: LGI.
#'   \item Insulin resistance: hba1c, hba1c_percent, eGDR, TyG.
#'   \item Lifestyle: smoking_status, drinking_status, diet_us, physical_activity, sleep_score.
#' }
#'
#' @param df A dataframe containing the required columns for calculation.
#' @param types A vector of indicator types, options are "eGDR", "TyG", "sleep_score", "LGI"...
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
#' @importFrom leo.basic leo_time_elapsed
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
  t0 <- Sys.time()
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
    # Socioeconomic Indicators (SES) ----
    "tdi" = list(required_cols     = c("tdi"),
                 required_cols_ukb = c("p22189"),     # Townsend index, recruitment
                 calc_function     = leo_tdi),
    "education" = list(required_cols     = c("education"),
                       required_cols_ukb = c("p6138_i0"),   # qualifications (instance 0)
                       calc_function     = leo_education),
    "household_income" = list(required_cols     = c("income"),
                              required_cols_ukb = c("p738_i0"),
                              calc_function     = leo_household_income),
    "household_income_2" = list(required_cols     = c("income"),
                                required_cols_ukb = c("p738_i0"),
                                calc_function     = leo_household_income_2),
    "career" = list(required_cols     = c("career"),
                    required_cols_ukb = c("p6142_i0"),   # current employment status (instance 0)
                    calc_function     = leo_career),
    # Inflammation ----
    "LGI" = list(required_cols      = c("crp", "wbc", "platelet", "neutrophil", "lymphocyte"),
                 required_cols_ukb  = c("p30710_i0", "p30000_i0", "p30080_i0", "p30140_i0", "p30120_i0"),
                 calc_function      = leo_LGI),
    # Insulin Resistance ----
    "hba1c" = list(required_cols      = c("hba1c"),
                   required_cols_ukb  = c("p30750_i0"),  # HbA1c (mmol/mol, IFCC) instance 0
                   calc_function      = leo_hba1c),
    "hba1c_percent" = list(required_cols      = c("hba1c"),
                           required_cols_ukb  = c("p30750_i0"),
                           calc_function      = leo_hba1c_percent),
    "eGDR" = list(required_cols = c("hba1c", "waist", "hypertension"),
                  required_cols_ukb = c("p30750_i0", "p48_i0", "hypertension"),
                  calc_function = leo_eGDR),
    "TyG" = list(required_cols = c("triglycerides", "glucose"),
                 required_cols_ukb = c("p30870_i0", "p30740_i0"),
                 calc_function = leo_TyG),
    # Life Style Indicators ----
    "smoking_status" = list(required_cols     = c("smoking_status"),
                            required_cols_ukb = c("p20116_i0"),  # smoking status (instance 0)
                            calc_function     = leo_smoking_status),
    "drinking_status" = list(required_cols     = c("drinking_status"),
                             required_cols_ukb = c("p20117_i0"),  # alcohol drinker status (instance 0)
                             calc_function     = leo_drinking_status),
    "diet_us" = list(required_cols      = c("fresh_fruit", "dried_fruit", "cooked_veg", "salad_veg",
                                            "oily_fish", "non_oily_fish",
                                            "bread_intake", "bread_type", "cereal_intake", "cereal_type",
                                            "processed_meat", "age_last_meat", "poultry", "beef", "lamb", "pork"),
                     required_cols_ukb  = c("p1309_i0", "p1319_i0", "p1289_i0", "p1299_i0",
                                            "p1329_i0", "p1339_i0",
                                            "p1438_i0", "p1448_i0", "p1458_i0", "p1468_i0",
                                            "p1349_i0", "p3680_i0", "p1359_i0", "p1369_i0", "p1379_i0", "p1389_i0"),
                     calc_function      = leo_diet_us),
    "physical_activity" = list(required_cols     = c("p884_i0","p894_i0","p904_i0","p914_i0"),
                               required_cols_ukb = c("p884_i0","p894_i0","p904_i0","p914_i0"),
                               calc_function     = leo_physical_activity),
    "sleep_score" = list(required_cols     = c("sleep_duration","chronotype","insomnia","snoring","daytime_sleepiness"),
                         required_cols_ukb = c("p1160_i0", "p1180_i0", "p1200_i0", "p1210_i0", "p1220_i0"),
                         calc_function     = leo_sleep_score)
    # Add more indicator types here if needed
  )

  # Loop through each type
  results <- df
  generated_types <- character(0)  # store successfully generated indicator names
  assist_cols_to_drop <- character(0)  # accumulate assist columns from all successful types
  for (type in types) {
    cli::cat_rule(paste("Generating new clinical indicator:", type), col = "blue")
    # Retrieve the required columns and function from the dictionary
    indicator_info <- indicator_dict[[type]]
    if (is.null(indicator_info)) { cli::cli_alert_danger("Unknown type: {type}. Skip."); next }
    required_cols <- if (type_id) indicator_info$required_cols_ukb else indicator_info$required_cols
    calc_function <- indicator_info$calc_function
    # check
    cli::cli_alert_info(paste("Checking if all required column is valid for {type} calculation."))
    missing_cols <- required_cols[!required_cols %in% colnames(df)]
    if (length(missing_cols) > 0) {
      cli::cli_alert_danger("Missing: {paste(missing_cols, collapse = ', ')}. Skip {type}"); next
    } else { cli::cli_alert_success("Pass for required columns checking.") }

    # Dynamically apply the function to calculate the indicator
    cli::cli_alert_info("Calculating the {type} using the following columns: {.val {required_cols}}")
    columns_to_pass <- lapply(required_cols, function(col) df[[col]])
    results <- results %>% mutate(!!type := do.call(calc_function, c(columns_to_pass, list(...))))

    # Track successfully generated types and accumulate their assist columns
    generated_types <- c(generated_types, type)
    assist_cols_to_drop <- c(assist_cols_to_drop, required_cols)
    cli::cli_alert_success(paste("{type} calculation completed."))
  }

  # Finalize the results
  if (remove_assist && !remove_other) {
    assist_cols_to_drop <- unique(assist_cols_to_drop)  # remove duplicates for clarity
    results <- results %>% dplyr::select(-all_of(assist_cols_to_drop))
  }
  if (remove_other) { results <- results %>% dplyr::select(eid, all_of(generated_types)) }
  cli::cat_rule("All done!", col = "blue")
  leo.basic::leo_time_elapsed(t0)
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

#' Recode BMI into categories (Underweight/Normal weight/Overweight/Obese)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001
#'
#' @param BMI Numeric/character vector of BMI (kg/m^2) from p21001.
#'
#' @return Factor with levels "Underweight", "Normal weight", "Overweight", "Obese";
#'   NA if invalid or missing.
#' @export
leo_BMI <- function(BMI, ...) {
  # convert to numeric and cut into 4 categories
  bmi_num <- suppressWarnings(as.numeric(BMI))
  bmi_num[bmi_num <= 0] <- NA

  bmi_cat <- dplyr::case_when(
    is.na(bmi_num)                  ~ NA_character_,
    bmi_num < 18.5                  ~ "Underweight",
    bmi_num >= 18.5 & bmi_num < 25  ~ "Normal weight",
    bmi_num >= 25   & bmi_num < 30  ~ "Overweight",
    bmi_num >= 30                   ~ "Obese"
  )
  bmi_cat <- factor(bmi_cat, levels = c("Underweight", "Normal weight", "Overweight", "Obese"))
  return(bmi_cat)
}

#' Recode sex (Male/Female)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31
#'
#' @param gender Numeric/character vector, usually from p31.
#'
#' @return Factor with levels "Female", "Male"; NA for negative/special codes.
#' @export
leo_gender <- function(gender, ...) {
  # keep only 0/1, others -> NA
  gender_num <- suppressWarnings(as.numeric(gender))
  gender_rec <- dplyr::case_when(
    gender_num == 0 ~ "Female",
    gender_num == 1 ~ "Male",
    TRUE            ~ NA_character_
  )
  gender_rec <- factor(gender_rec, levels = c("Female", "Male"))
  cli::cli_alert_info("Sex recoding:\n  Female (0)\n  Male (1)\n")
  return(gender_rec)
}

#' Recode ethnicity to White vs non-White
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000
#'
#' @param ethnicity Numeric/character vector, usually from p21000.
#'
#' @return Factor with levels "White", "Non-White";
#'   NA for -1/-3/missing.
#' @export
leo_ethnicity <- function(ethnicity, ...) {
  eth_num <- suppressWarnings(as.numeric(ethnicity))
  white_codes <- c(1, 1001, 1002, 1003)

  eth_rec <- dplyr::case_when(
    eth_num %in% white_codes                ~ "White",
    eth_num > 0                             ~ "Non-White",
    eth_num %in% c(-1, -3) | is.na(eth_num) ~ NA_character_,  # Unknown / PNSA
    TRUE                                    ~ NA_character_
  )

  # eth_rec <- factor(eth_rec, levels = c("White", "Non-White"))
  cli::cli_alert_info("Ethnicity recoding:\n  White (1/1001/1002/1003)\n  Non-White (other >0 codes)\n  -1/-3/NA -> NA")
  return(eth_rec)
}

#' Recode ethnicity into finer categories (top-level 1–6)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000
#'
#' @param ethnicity Numeric/character vector, usually from p21000.
#'
#' @return Factor with levels:
#'   "White", "Mixed", "Asian or Asian British", "Black or Black British",
#'   "Chinese", "Other ethnic group"; NA for -1/-3/missing/other unexpected codes.
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
    eth_num %in% c(-1, -3) | is.na(eth_num) ~ NA_character_,  # unknown / PNSA
    eth_num %in% white_codes                ~ "White",
    eth_num %in% mixed_codes                ~ "Mixed",
    eth_num %in% asian_codes                ~ "Asian or Asian British",
    eth_num %in% black_codes                ~ "Black or Black British",
    eth_num %in% chinese_codes              ~ "Chinese",
    eth_num %in% other_codes                ~ "Other ethnic group",
    TRUE                                    ~ NA_character_   # any other unexpected code
  )

  eth_rec <- factor(
    eth_rec,
    levels = c("White", "Mixed", "Asian or Asian British", "Black or Black British", "Chinese", "Other ethnic group")
  )
  return(eth_rec)
}

# Socioeconomic Indicators (SES) ----
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

#' Recode UKB qualifications into 3 categories
#'
#' Recode UKB field 6138 (Qualifications; multiple selections separated by "|")
#' into Degree / No degree / Unknown.
#'
#' @param p6138_i0 UKB 6138 qualifications string (https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6138).
#' @param ... Reserved for future use.
#'
#' @return Character vector with values: "Degree", "No degree", or "Unknown".
#' @export
#' @importFrom dplyr case_when
#' @importFrom stringr str_detect
#' @importFrom cli cli_alert_info cli_alert_success
leo_education <- function(p6138_i0, ...) {
  cli::cli_alert_info("Recoding UKB 6138 into Degree/No degree/Unknown")

  x <- as.character(p6138_i0)
  unknown <- is.na(x) | x == "" | x == "-3"
  has_degree <- !unknown & stringr::str_detect(paste0("|", x, "|"), "\\|1\\|")

  education <- dplyr::case_when(unknown ~ NA_character_,
                                has_degree ~ "Degree",
                                T ~ "Non_degree")

  cli::cli_alert_success("Education recoding completed.")
  return(education)
}

#' Recode employment status
#'
#' Recode employment status into paid vs not-paid employment
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6142
#'
#' @param p6142_i0 Current employment status (UKB 6142, instance 0).
#'   Values:
#'   1 = In paid employment or self-employed
#'   2 = Retired
#'   3 = Looking after home and/or family
#'   4 = Unable to work because of sickness or disability
#'   5 = Unemployed
#'   6 = Doing unpaid or voluntary work
#'   7 = Full or part-time student
#'   -7 = None of the above
#'   -3 = Prefer not to answer
#' @param ... Reserved for future use.
#'
#' @return Factor with categories "Paid employment" and "Not paid employment"; NA for special codes.
#' @export
leo_career <- function(p6142_i0, ...) {
  x <- suppressWarnings(as.numeric(p6142_i0))
  x[x %in% c(-3, -7)] <- NA_real_
  career <- dplyr::case_when(
    is.na(x) ~ NA_character_,
    x == 1 ~ "Paid employment",
    x %in% 2:7 ~ "Not paid employment",
    TRUE ~ NA_character_
  )
  return(factor(career, levels = c("Paid employment", "Not paid employment")))
}

#' Recode household income (UKB 738) into readable categories
#'
#' Recode UKB field 738 "Average total household income before tax" into
#' 5 brackets plus "Unknown" (-1/-3/NA/blank).
#'
#' @param p738_i0 Household income before tax (UKB 738, instance 0).
#'   Values: 1–5; -1/-3 -> Unknown.
#' @param ... Reserved for future use.
#'
#' @return Character vector: "0-18,000", "18,000-30,999", "31,000-51,999",
#'   "52,000-100,000", "Above 100,000", or "Unknown".
#' @export
#' @importFrom cli cli_alert_info
#' @importFrom dplyr case_when
leo_household_income <- function(p738_i0, ...) {
  income_code <- suppressWarnings(as.numeric(p738_i0))
  income_code[income_code %in% c(-1, -3)] <- NA_real_

  income <- dplyr::case_when(
    is.na(income_code) ~ NA_character_,
    income_code == 1   ~ "£0–£18,000",
    income_code == 2   ~ "£18,000–£30,999",
    income_code == 3   ~ "£31,000–£51,999",
    income_code == 4   ~ "£52,000–£100,000",
    income_code == 5   ~ "£100,000 and above",
    TRUE               ~ NA_character_
  )
  income <- factor(income,
                   levels = c("£0–£18,000", "£18,000–£30,999", "£31,000–£51,999",
                              "£52,000–£100,000", "£100,000 and above")
  )
  return(income)
}

#' Recode household income (UKB 738) into readable categories (2 categories)
#'
#' The five levels of the average total before-tax household income were classified
#' as “less than” or “equal to or above” £31,000, which is closest to the UK median
#' household income in October 2009 (£27,530)
#' - Reference: *Air Pollution, Genetic Factors, and the Risk of Lung Cancer A Prospective Study in the UK Biobank*
#'
#' @param p738_i0 Household income before tax (UKB 738, instance 0).
#' @param ... Reserved for future use.
#'
#' @return Character vector
#' @export
#' @importFrom cli cli_alert_info
#' @importFrom dplyr case_when
leo_household_income_2 <- function(p738_i0, ...) {
  income_code <- suppressWarnings(as.numeric(p738_i0))
  income_code[income_code %in% c(-1, -3)] <- NA_real_

  income <- dplyr::case_when(
    is.na(income_code) ~ NA_character_,
    income_code == 1   ~ "Less than £31,000",
    income_code == 2   ~ "Less than £31,000",
    income_code == 3   ~ "£31,000 and above",
    income_code == 4   ~ "£31,000 and above",
    income_code == 5   ~ "£31,000 and above",
    TRUE               ~ NA_character_
  )
  income <- factor(income,
                   levels = c("Less than £31,000", "£31,000 and above"))
  return(income)
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
  .LGI_decile_score <- function(x) {
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
  crp_score <- .LGI_decile_score(crp)
  wbc_score <- .LGI_decile_score(wbc)
  plt_score <- .LGI_decile_score(plt)
  nlr_score <- .LGI_decile_score(nlr)

  # sum component scores
  lgi_score <- crp_score + wbc_score + plt_score + nlr_score

  return(lgi_score)
}

# Insulin Resistance -----

#' Recode HbA1c (IFCC, mmol/mol)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30750
#'
#' UK Biobank HbA1c is stored in mmol/mol (IFCC). This function performs:
#' - numeric coercion,
#' - recoding UKB special codes (-3/-1) to NA,
#' - dropping non-positive values to NA (conservative guard).
#'
#' @param hba1c Numeric/character vector of HbA1c values (mmol/mol, IFCC), usually from p30750.
#' @param ... Reserved for future use.
#'
#' @return Numeric vector of HbA1c in mmol/mol (IFCC), with special codes as NA.
#' @export
#'
#' @importFrom cli cli_alert_info
leo_hba1c <- function(hba1c, ...) {
  x <- suppressWarnings(as.numeric(hba1c))
  cli::cli_alert_info("HbA1c recoding (IFCC, mmol/mol)")
  return(x)
}

#' Convert HbA1c from IFCC (mmol/mol) to NGSP (%)
#' https://ngsp.org/ifcc.asp
#'
#' Conversion formula:
#' \deqn{NGSP(\%) = 0.09148 \times IFCC(mmol/mol) + 2.152}
#'
#' @param hba1c Numeric/character vector of HbA1c values (mmol/mol, IFCC), usually from p30750.
#' @param ... Reserved for future use.
#'
#' @return Numeric vector of HbA1c in % (NGSP), NA preserved.
#' @export
#'
#' @importFrom cli cli_alert_info
leo_hba1c_percent <- function(hba1c, ...) {
  cli::cli_alert_info("Converting HbA1c from mmol/mol(IFCC) to %(NGSP) >>> {.code NGSP=[0.09148*IFCC]+2.152} {.url https://ngsp.org/ifcc.asp}")
  hba1c_ifcc <- leo_hba1c(hba1c)
  hba1c_ngsp <- 0.09148 * hba1c_ifcc + 2.152
  return(hba1c_ngsp)
}

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
#' - **hypertension** need to extract from icd10 or first occurrence records first.
#' @references
#' 1. \url{https://diabetesjournals.org/care/article/36/8/2280/32950/Use-of-the-Estimated-Glucose-Disposal-Rate-as-a} (Diabetes Care, 2013 July)
#' 2. \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC11439291/#Sec2} (Cardiovasc Diabetol, 2024 Sep)
leo_eGDR <- function(hba1c, waist, hypertension, ...) {
  hba1c_percent <- leo_hba1c_percent(hba1c)
  cli::cli_alert_info("Calculating eGDR using the formula: {.code eGDR = 21.158 + (-0.09 * waist) + (-3.407 * hypertension) + (-0.551 * hba1c_percent)}")
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
  cli::cli_alert_info("Converting triglycerides and glucose from {.emph \"mmmol/L\"} (UKB default unit) to {.emph \"mg/dL\"}")
  # Convert triglycerides and glucose from mmol/L to mg/dL
  triglycerides_mg_dl <- triglycerides * 88.5704  # Exact conversion factor
  glucose_mg_dl <- glucose * 18.0168  # Exact conversion factor
  # Calculate TyG index element-wise for the vectors
  cli::cli_alert_info("Calculating TyG index using the formula: {.code TyG = ln[(triglycerides * glucose)/2]}")
  TyG <- log((triglycerides_mg_dl * glucose_mg_dl) / 2)
  return(TyG)
}

# Life style Indicators -----

#' Recode smoking status (Never/Previous/Current)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116
#'
#' @param smoking_status Numeric/character vector, usually from p20116.
#'
#' @return Factor with levels "Never", "Previous", "Current";
#'   NA for negative/special codes or others.
#' @export
leo_smoking_status <- function(smoking_status, ...) {
  # keep 0/1/2, recode -3/-1 and other values to NA
  x <- suppressWarnings(as.numeric(smoking_status))
  x[x %in% c(-3, -1)] <- NA
  x[!is.na(x) & !(x %in% c(0, 1, 2))] <- NA
  smoking <- dplyr::case_when(
    x == 0 ~ "Never",
    x == 1 ~ "Previous",
    x == 2 ~ "Current",
    TRUE   ~ NA_character_
  )
  smoking <- factor(smoking, levels = c("Never", "Previous", "Current"))
  return(smoking)
}

#' Recode drinking status (Never/Previous/Current)
#' https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20117
#'
#' @param drinking_status Numeric/character vector, usually from p20117.
#'
#' @return Factor with levels "Never", "Previous", "Current";
#'   NA for negative/special codes or others.
#' @export
leo_drinking_status <- function(drinking_status, ...) {
  # keep 0/1/2, recode -3/-1 and other values to NA
  x <- suppressWarnings(as.numeric(drinking_status))
  x[x %in% c(-3, -1)] <- NA
  x[!is.na(x) & !(x %in% c(0, 1, 2))] <- NA
  drinking <- dplyr::case_when(
    x == 0 ~ "Never",
    x == 1 ~ "Previous",
    x == 2 ~ "Current",
    TRUE   ~ NA_character_
  )
  drinking <- factor(drinking, levels = c("Never", "Previous", "Current"))
  return(drinking)
}

#' Recode drinking consumption
#'
#' This is more sophisticated indictor shoing if one is exccessively drink,
#' cobining thinking: red_wine white wine, beer, spirits, and fortified wine Numeric/character vector, usually from p20117.
#'
#' @return Integer vector (0/1/2; NA for negative/special codes or others).
#' @export
leo_drinking_consumption <- function(...) {
  # todo：
  # 202601: 主函数还未加
  return(x)
}

#' healthy diet score by U.S. Dietary Guidelines
#'
#' Calculate a 0–7 healthy diet score using seven components (vegetables, fruit,
#' fish, whole grains, refined grains, processed meat, unprocessed meat) following
#' the U.S. Dietary Guidelines targets. Each component meeting its intake target
#' contributes 1 point; higher scores indicate healthier patterns.
#'
#' Field expectations (UKB instance 0):
#' - Fruits: p1309_i0 (fresh pieces/day), p1319_i0 (dried pieces/day; 5 pieces = 1 serving)
#' - Vegetables: p1289_i0 (cooked tablespoons/day), p1299_i0 (salad tablespoons/day; 3 tbsp = 1 serving)
#' - Fish: p1329_i0 (oily fish frequency code), p1339_i0 (non-oily fish frequency code)
#' - Bread: p1438_i0 (slices/week, -10 = "less than once"), p1448_i0 (type: 1=White, 2=Brown, 3=Wholemeal, 4=Other)
#' - Cereal: p1458_i0 (bowls/week, -10 = "less than once"), p1468_i0 (type: 1=Bran, 2=Biscuit, 3=Oat, 4=Muesli, 5=Other)
#' - Processed meat: p1349_i0 (frequency code 0–5, per Data-Coding 100377); p3680_i0 == 0 implies never ate meat
#' - Unprocessed meat: p1359_i0 (poultry), p1369_i0 (beef), p1379_i0 (lamb), p1389_i0 (pork) — all frequency codes 0–5
#'
#' All frequency fields (100377) encode: 0=Never, 1=<1/wk, 2=1/wk, 3=2–4/wk, 4=5–6/wk, 5=Daily.
#' To meet component thresholds, frequency codes must be converted to weekly servings
#' (e.g., code 2 → 1/wk, code 3 → 3/wk midpoint, code 5 → 7/wk).
#'
#' Whole grains: bread_type==3 (wholemeal) + cereal_type in (1,3,4) (bran/oat/muesli)
#' Refined grains: bread_type in (1,2,4) (white/brown/other) + cereal_type in (2,5) (biscuit/other)
#'
#' @param fresh_fruit Fresh fruit pieces per day.
#' @param dried_fruit Dried fruit pieces per day (5 pieces = 1 serving).
#' @param cooked_veg Cooked vegetables tablespoons per day (3 tbsp = 1 serving).
#' @param salad_veg Salad/raw vegetables tablespoons per day (3 tbsp = 1 serving).
#' @param oily_fish Oily fish frequency code (p1329_i0, 0–5 per Data-Coding 100377).
#' @param non_oily_fish Non-oily fish frequency code (p1339_i0, 0–5 per Data-Coding 100377).
#' @param bread_intake Bread slices per week (p1438_i0); -10 treated as 0/week (less than once).
#' @param bread_type Bread type code (p1448_i0): 1=White, 2=Brown, 3=Wholemeal, 4=Other.
#' @param cereal_intake Cereal bowls per week (p1458_i0); -10 treated as 0/week.
#' @param cereal_type Cereal type code (p1468_i0): 1=Bran, 2=Biscuit, 3=Oat, 4=Muesli, 5=Other.
#' @param processed_meat Processed meat frequency code (p1349_i0, 0–5 per Data-Coding 100377).
#' @param age_last_meat Age when last ate meat (p3680_i0); 0 = never ate meat.
#' @param poultry Poultry frequency code (p1359_i0, 0–5 per Data-Coding 100377).
#' @param beef Beef frequency code (p1369_i0, 0–5 per Data-Coding 100377).
#' @param lamb Lamb/mutton frequency code (p1379_i0, 0–5 per Data-Coding 100377).
#' @param pork Pork frequency code (p1389_i0, 0–5 per Data-Coding 100377).
#' @param ... Reserved for future use.
#'
#' @return Integer vector in [0,7]; NA if any component cannot be evaluated for a row.
#' @export
#' @importFrom cli cli_alert_info
leo_diet_us <- function(fresh_fruit, dried_fruit, cooked_veg, salad_veg,
                        oily_fish, non_oily_fish,
                        bread_intake, bread_type, cereal_intake, cereal_type,
                        processed_meat, age_last_meat, poultry, beef, lamb, pork, ...) {
  cli::cli_alert_info("Calculating healthy diet score (U.S. Dietary Guidelines; 0-7)")

  # component 1: Fruits
  fruit_fresh_p1309 <- recode_na31neg10(fresh_fruit) # UKB 1309
  fruit_dried_p1319 <- recode_na31neg10(dried_fruit) # UKB 1319
  # component 2: Vegetables
  veg_cooked_p1289 <- recode_na31neg10(cooked_veg)   # UKB 1289
  veg_salad_p1299  <- recode_na31neg10(salad_veg)    # UKB 1299
  # component 3: Fish
  fish_oily_p1329     <- recode_na31neg10(oily_fish)     # UKB 1329
  fish_non_oily_p1339 <- recode_na31neg10(non_oily_fish) # UKB 1339
  # component 4-5: Whole grains / Refined grains
  bread_slices_p1438 <- recode_na31neg10(bread_intake)  # UKB 1438
  bread_type_p1448   <- recode_na31(bread_type)         # UKB 1448
  cereal_bowls_p1458 <- recode_na31neg10(cereal_intake) # UKB 1458
  cereal_type_p1468  <- recode_na31(cereal_type)        # UKB 1468
  # component 6: Processed meat
  processed_meat_p1349 <- recode_na31neg10(processed_meat) # UKB 1349
  meat_last_age_p3680  <- recode_na31(age_last_meat)       # UKB 3680
  never_meat <- !is.na(meat_last_age_p3680) & meat_last_age_p3680 == 0
  # component 7: Unprocessed meat
  poultry_p1359 <- recode_na31neg10(poultry) # UKB 1359
  beef_p1369    <- recode_na31neg10(beef)    # UKB 1369
  lamb_p1379    <- recode_na31neg10(lamb)    # UKB 1379
  pork_p1389    <- recode_na31neg10(pork)    # UKB 1389


  # Frequency code -> weekly servings (Data-Coding 100377: 0=Never, 1=<1/wk, 2=1/wk, 3=2-4/wk, 4=5-6/wk, 5=Once or more Daily; -1/-3 handled as NA)
  # Data-Coding 100377 mapping (https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100377)
  # Logic we used here is to treat "less than once" as zero intake, for other intake intervals we use midpoints.
  freq_map <- c(`0` = 0, `1` = 0, `2` = 1, `3` = 3, `4` = 5.5, `5` = 7)

  # Component 1: Fruits (UKB 1309, 1319) - >=3 servings/day; dried fruit 5 pieces = 1 serving
  fruit_serv_day <- rowSums(cbind(fruit_fresh_p1309, fruit_dried_p1319 / 5), na.rm = TRUE)
  fruit_serv_day[is.na(fruit_fresh_p1309) & is.na(fruit_dried_p1319)] <- NA_real_

  # Component 2: Vegetables (UKB 1289, 1299) - >=3 servings/day; 3 tbsp = 1 serving
  veg_serv_day <- rowSums(cbind(veg_cooked_p1289 / 3, veg_salad_p1299 / 3), na.rm = TRUE)
  veg_serv_day[is.na(veg_cooked_p1289) & is.na(veg_salad_p1299)] <- NA_real_

  # Component 3: Fish (UKB 1329, 1339) - Data-Coding 100377; >=2 servings/week
  fish_oily_week <- as.numeric(freq_map[as.character(fish_oily_p1329)])
  fish_non_oily_week <- as.numeric(freq_map[as.character(fish_non_oily_p1339)])
  fish_serv_week <- rowSums(cbind(fish_oily_week, fish_non_oily_week), na.rm = TRUE)
  fish_serv_week[is.na(fish_oily_week) & is.na(fish_non_oily_week)] <- NA_real_

  # Component 4: Whole grains (UKB 1438/1448, 1458/1468) - >=3 servings/day
  # Component 5: Refined grains (UKB 1438/1448, 1458/1468) - <=2 servings/day
  # Only require type if intake > 0.
  whole_bread_week <- dplyr::case_when(
    is.na(bread_slices_p1438) ~ NA_real_,
    bread_slices_p1438 == 0   ~ 0,
    bread_slices_p1438 > 0 & is.na(bread_type_p1448) ~ NA_real_,
    bread_slices_p1438 > 0 & bread_type_p1448 == 3   ~ bread_slices_p1438,
    TRUE                      ~ 0
  )

  refined_bread_week <- dplyr::case_when(
    is.na(bread_slices_p1438) ~ NA_real_,
    bread_slices_p1438 == 0   ~ 0,
    bread_slices_p1438 > 0 & is.na(bread_type_p1448) ~ NA_real_,
    bread_slices_p1438 > 0 & bread_type_p1448 %in% c(1, 2, 4) ~ bread_slices_p1438,
    TRUE                      ~ 0
  )

  whole_cereal_week <- dplyr::case_when(
    is.na(cereal_bowls_p1458) ~ NA_real_,
    cereal_bowls_p1458 == 0   ~ 0,
    cereal_bowls_p1458 > 0 & is.na(cereal_type_p1468) ~ NA_real_,
    cereal_bowls_p1458 > 0 & cereal_type_p1468 %in% c(1, 3, 4) ~ cereal_bowls_p1458,
    TRUE                      ~ 0
  )

  refined_cereal_week <- dplyr::case_when(
    is.na(cereal_bowls_p1458) ~ NA_real_,
    cereal_bowls_p1458 == 0   ~ 0,
    cereal_bowls_p1458 > 0 & is.na(cereal_type_p1468) ~ NA_real_,
    cereal_bowls_p1458 > 0 & cereal_type_p1468 %in% c(2, 5) ~ cereal_bowls_p1458,
    TRUE                      ~ 0
  )

  whole_grain_week <- rowSums(cbind(whole_bread_week, whole_cereal_week), na.rm = TRUE)
  whole_grain_week[is.na(whole_bread_week) & is.na(whole_cereal_week)] <- NA_real_
  whole_grain_day <- ifelse(is.na(whole_grain_week), NA_real_, whole_grain_week / 7)

  refined_grain_week <- rowSums(cbind(refined_bread_week, refined_cereal_week), na.rm = TRUE)
  refined_grain_week[is.na(refined_bread_week) & is.na(refined_cereal_week)] <- NA_real_
  refined_grain_day <- ifelse(is.na(refined_grain_week), NA_real_, refined_grain_week / 7)

  # Component 6: Processed meat (UKB 1349, 3680) - Data-Coding 100377; <=1 serving/week; 3680==0 => never meat
  processed_meat_week <- as.numeric(freq_map[as.character(processed_meat_p1349)])
  processed_meat_week[never_meat %in% TRUE] <- 0

  # Component 7: Unprocessed meat (UKB 1359/1369/1379/1389, 3680) - Data-Coding 100377; <=2 servings/week; 3680==0 => never meat
  poultry_week <- as.numeric(freq_map[as.character(poultry_p1359)])
  beef_week    <- as.numeric(freq_map[as.character(beef_p1369)])
  lamb_week    <- as.numeric(freq_map[as.character(lamb_p1379)])
  pork_week    <- as.numeric(freq_map[as.character(pork_p1389)])
  poultry_week[never_meat %in% TRUE] <- 0
  beef_week[never_meat %in% TRUE] <- 0
  lamb_week[never_meat %in% TRUE] <- 0
  pork_week[never_meat %in% TRUE] <- 0
  unprocessed_meat_week <- rowSums(cbind(poultry_week, beef_week, lamb_week, pork_week), na.rm = TRUE)
  unprocessed_meat_week[is.na(poultry_week) & is.na(beef_week) & is.na(lamb_week) & is.na(pork_week)] <- NA_real_

  component_scores <- list(
    vegetables       = ifelse(is.na(veg_serv_day), NA_integer_, ifelse(veg_serv_day >= 3, 1L, 0L)),
    fruit            = ifelse(is.na(fruit_serv_day), NA_integer_, ifelse(fruit_serv_day >= 3, 1L, 0L)),
    fish             = ifelse(is.na(fish_serv_week), NA_integer_, ifelse(fish_serv_week >= 2, 1L, 0L)),
    whole_grains     = ifelse(is.na(whole_grain_day), NA_integer_, ifelse(whole_grain_day >= 3, 1L, 0L)),
    refined_grains   = ifelse(is.na(refined_grain_day), NA_integer_, ifelse(refined_grain_day <= 2, 1L, 0L)),
    processed_meat   = ifelse(is.na(processed_meat_week), NA_integer_, ifelse(processed_meat_week <= 1, 1L, 0L)),
    unprocessed_meat = ifelse(is.na(unprocessed_meat_week), NA_integer_, ifelse(unprocessed_meat_week <= 2, 1L, 0L))
  )

  comp_mat <- do.call(cbind, component_scores)
  any_na <- apply(comp_mat, 1, function(r) any(is.na(r)))
  total <- rowSums(comp_mat, na.rm = TRUE)
  total[any_na] <- NA_integer_

  return(as.integer(total))
}

#' Calculate Physical Activity Indicator
#' `r lifecycle::badge('stable')`
#'
#' Here, using IPAQ minutes/week as input, classify regular activity according
#' to AHA (150/75/equivalent combination) rules:
#' \itemize{
#'   \item Moderate:  \eqn{\ge 150} min/week; or
#'   \item Vigorous:  \eqn{\ge 75} min/week; or
#'   \item Equivalent: \eqn{moderate + 2 \times vigorous \ge 150.}
#' }
#'
#' @note *Adults should pursue at least 150 minutes per week of moderate-intensity
#' physical activity, or 75 minutes per week of vigorous-intensity aerobic
#' physical activity, or an equivalent combination of moderate- and
#' vigorous-intensity aerobic activities.* (Circulation, 2010, [source](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.109.192703))
#' - **≥150 min/wk moderate intensity**
#' - **≥75 min/wk vigorous intensity**
#' - **Combination**
#'
#' @param p884_i0 Moderate activity days/week (10+ min). UKB 884; -1/-3 -> NA.
#' @param p894_i0 Moderate activity minutes/day (typical day). UKB 894; -1/-3 -> NA.
#' @param p904_i0 Vigorous activity days/week (10+ min). UKB 904; -1/-3 -> NA.
#' @param p914_i0 Vigorous activity minutes/day (typical day). UKB 914; -1/-3 -> NA.
#' @param ... Reserved for future use.
#'
#' @return Factor with levels c("Non-regular", "Regular"); NA if insufficient data.
#' @export
#' @importFrom cli cli_alert_info
#' @importFrom dplyr case_when coalesce
leo_physical_activity <- function(p884_i0, p894_i0, p904_i0, p914_i0, ...) {
  p884_i0 <- recode_na31(p884_i0); p894_i0 <- recode_na31(p894_i0)
  p904_i0 <- recode_na31(p904_i0); p914_i0 <- recode_na31(p914_i0)
  moderate_min_week <- dplyr::case_when(is.na(p884_i0) ~ NA_real_, p884_i0 == 0 ~ 0,
                                        is.na(p894_i0) ~ NA_real_, T ~ p884_i0 * p894_i0)
  vigorous_min_week <- dplyr::case_when(is.na(p904_i0) ~ NA_real_, p904_i0 == 0 ~ 0,
                                        is.na(p914_i0) ~ NA_real_, T ~ p904_i0 * p914_i0)

  moderate_ok <- !is.na(moderate_min_week) & moderate_min_week >= 150
  vigorous_ok <- !is.na(vigorous_min_week) & vigorous_min_week >= 75
  equivalent_min_week <- dplyr::coalesce(moderate_min_week, 0) + 2 * dplyr::coalesce(vigorous_min_week, 0)
  equivalent_ok <- (!is.na(moderate_min_week) | !is.na(vigorous_min_week)) & equivalent_min_week >= 150

  activity <- dplyr::case_when(is.na(moderate_min_week) & is.na(vigorous_min_week) ~ NA_character_,
                               moderate_ok | vigorous_ok | equivalent_ok ~ "Regular",
                               T ~ "Non-regular") %>%
    factor(levels = c("Non-regular", "Regular"))
  return(activity)
}

#' Calculate Sleep Score and Pattern
#'
#' Compute a 0-5 sleep health score based on 5 components:
#' - Sleep duration (7-8h)
#' - Chronotype (morning person)
#' - No insomnia
#' - No snoring
#' - No daytime sleepiness
#'
#' @param p1160_i0 Sleep duration (hours). UKB 1160; 7-8h scores 1.
#' @param p1180_i0 Chronotype. UKB 1180; 1/2 (morning) scores 1; -1/-3 -> NA.
#' @param p1200_i0 Insomnia. UKB 1200; 1 (never/rarely) scores 1; -3 -> NA.
#' @param p1210_i0 Snoring. UKB 1210; 2 (no) scores 1; -1/-3 -> NA.
#' @param p1220_i0 Daytime sleepiness. UKB 1220; 0/1 (never/rarely) scores 1; -1/-3 -> NA.
#' @param ... Reserved for future use.
#'
#' @return Factor with levels c("poor", "intermediate", "healthy");
#'   based on score: 0-1=poor, 2-3=intermediate, 4-5=healthy.
#' @export
#' @importFrom cli cli_alert_info
#' @importFrom dplyr case_when
leo_sleep_score <- function(p1160_i0, p1180_i0, p1200_i0, p1210_i0, p1220_i0, ...) {
  cli::cli_alert_info("Calculating sleep score (0-5) and pattern classification")
  p1160 <- recode_na31(p1160_i0)
  p1180 <- recode_na31(p1180_i0)
  p1200 <- recode_na31(p1200_i0)
  p1210 <- recode_na31(p1210_i0)
  p1220 <- recode_na31(p1220_i0)
  score_duration <- dplyr::case_when(is.na(p1160) ~ NA_integer_,
                                     p1160 >= 7 & p1160 <= 8 ~ 1L, # 7-8 hours sleep
                                     TRUE ~ 0L)
  score_chronotype <- dplyr::case_when(is.na(p1180) ~ NA_integer_,
                                       p1180 %in% c(1, 2) ~ 1L,  # Morning person
                                       TRUE ~ 0L)
  score_insomnia <- dplyr::case_when(is.na(p1200) ~ NA_integer_,
                                     p1200 == 1 ~ 1L,  # Never/rarely
                                     TRUE ~ 0L)
  score_snoring <- dplyr::case_when(is.na(p1210) ~ NA_integer_,
                                    p1210 == 2 ~ 1L,  # No snoring
                                    TRUE ~ 0L)
  score_sleepiness <- dplyr::case_when(is.na(p1220) ~ NA_integer_,
                                       p1220 %in% c(0, 1) ~ 1L,  # Never/rarely
                                       TRUE ~ 0L)

  score_matrix <- cbind(score_duration, score_chronotype, score_insomnia, score_snoring, score_sleepiness)
  sleep_score <- rowSums(score_matrix, na.rm = FALSE) # If one is NA, `na.rm = FALSE` make sure the final is NA

  sleep_pattern <- dplyr::case_when(is.na(sleep_score) ~ NA_character_,
                                    sleep_score >= 4 ~ "healthy",      # 4-5分
                                    sleep_score >= 2 ~ "intermediate", # 2-3分
                                    TRUE ~ "poor")                     # 0-1分

  sleep_pattern <- factor(sleep_pattern, levels = c("poor", "intermediate", "healthy"))
  cli::cli_alert_info("Sleep pattern distribution:\n  Poor: {sum(sleep_pattern == 'poor', na.rm=T)}\n  Intermediate: {sum(sleep_pattern == 'intermediate', na.rm=T)}\n  Healthy: {sum(sleep_pattern == 'healthy', na.rm=T)}")
  return(sleep_pattern)
}



# Helper functions ----
#' Helper function to recode UKB special values (-3/-1) to NA
#' @details Used for UK Biobank data-codings where -1 = "Do not know" and -3 = "Prefer not to answer".
#' @keywords internal
recode_na31 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[x_num %in% c(-3, -1)] <- NA
  x_num
}

#' Helper function to recode UKB special values (-3/-1 to NA and -10 to 0)
#' @details For UK Biobank fields where -10 encodes "Less than once" and should be treated as 0;
#'   -1 = "Do not know", -3 = "Prefer not to answer".
#' @keywords internal
recode_na31neg10 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  x_num[x_num %in% c(-3, -1)] <- NA
  x_num[x_num == -10] <- 0
  x_num
}
