# Here, we build table 1 for a clincial cohort study.

#' Table 1 Help
#'
#' This function provides help on how to use the leo.table1 functions and basic background knowledge.
#
#' @importFrom cli cli_h1 cli_h3 cat_line cat_bullet cli_code
#' @note
#' If you encounter a warning about an unknown RStudio theme when using `cli` functions, you can suppress it by setting:
#'
#' ```r
#' options(cli.ignore_unknown_rstudio_theme = TRUE)
#' ```
#' This option should be set in your R session to avoid the warning.
#' @export
leo.table1.help <- function() {
  cli::cli_h1("Table 1 Functions Overview & Background")
  # step 1
  cli::cli_h3("1. `leo.table1.step1`: Normality Check")
  cli::cat_line("- Sample > 5000: Use Kolmogorov-Smirnov test for normality check.")
  cli::cat_line("- Sample <= 5000: Use Shapiro-Wilk test for normality check.")
  cli::cat_line("- Sample < 3: No normality test performed.")
  # step 2
  cli::cli_h3("\n2. `leo.table1.step2`: Descriptive Stats & Group Comparisons")
  cli::cat_bullet("Descriptive stats:")
  cli::cat_line("- Mean ± SD for continuous variables (normal distribution).")
  cli::cat_line("- Median (IQR) for continuous variables (non-normal distribution).")
  cli::cat_line("- Count (proportion) for categorical variables.")
  cli::cat_bullet("Group comparisons:")
  cli::cat_line("- Independent T-test: Used for comparing means of two independent groups when continuous variables follow a normal distribution.")
  cli::cat_line("- Mann-Whitney U-test: Used for comparing medians of two independent groups when continuous variables do not follow a normal distribution.")
  cli::cat_line("- Chi-squared test: Used for comparing proportions of categorical variables between two or more independent groups.")
  cli::cat_line("- Fisher's exact test: Used for comparing proportions of categorical variables when sample size is small (expected cell count < 5).")
  # step 3
  cli::cli_h3("\n3. `leo.table1.save`: Export")
  cli::cat_line("- Export the table as a Word document or CSV file.")
  # example code
  cli::cli_h1("\nExample Code: \n")
  eg_code <- '
    pacman::p_load(tableone, data.table)
    cohort <- fread("xxx.csv") # clinical cohort
    cohort_num_var <- leo.table1.step1(cohort, num_var = c("age", "tdi", "TyG")) # check normality for numeric variables
    table_1 <- leo.table1.step2(
      cohort,
      var_all = c("age", "gender", "ethnicity"),  # all variables to include in table 1
      var_cat = c("gender", "ethnicity"),         # categorical variables (can also use set_diff)
      var_non = c("age"),                         # non-normal continuous variables
      strata = "exposure_status",                 # stratify by exposure status (optional)
      compare_test = F, includeNA = F, showAllLevels = F)
    leo.table1.save(table_1, "table1.csv")        # save table 1 as a CSV file
    leo.table1.save(table_1, "table1.docx")       # save table 1 as a Word document
    '
  cli::cli_code(eg_code, language = "R")
}



#' Step 1: Check Normality for Variables
#'
#' This function checks the normality of each variable in the dataset and returns a summary.
#'
#' @param df A data frame or tibble containing the variables to be tested.
#' @param id A character string specifying the column name for patient ID. Default is "eid".
#' @param num_var user-provided numeric variables to test for normality. Default is NULL.
#'
#' @return A data frame with variables, p-values for normality test, and normality status.
#' @importFrom stats shapiro.test ks.test
#' @importFrom dplyr mutate select %>% n_distinct pull
#' @importFrom cli cli_alert_success cli_alert_warning
#' @importFrom nortest lillie.test
#' @export
#' @examples
#' # leo.table1.step1(t1d_cohort, id_to_exclude = c("eid", "iri_time"), num_var = c("age", "tdi", "score_diet", "TyG"))
leo.table1.step1 <- function(df, id_to_exclude = c("eid", "iri_time"), num_var = NULL) {
  # Filter numeric variables with >2 unique values
  if (is.null(num_var)) {
    num_var <- df %>%
      dplyr::select(-all_of(id_to_exclude)) %>%
      dplyr::select(where(~ is.numeric(.) && dplyr::n_distinct(.) > 2)) %>%
      names()
    cli::cli_alert_info("Automatically detected {length(num_var)} numeric variable{?s}: {paste(num_var, collapse = ', ')}")
  } else {
    cli::cli_alert_info("User-provided {length(num_var)} numeric variable{?s}: {paste(num_var, collapse = ', ')}")
  }

  # normality test
  res <- data.frame(var = num_var,
                    p_value = rep(NA, length(num_var)),
                    normality = rep(NA, length(num_var)),
                    method = rep(NA, length(num_var)))

  for (i in 1:nrow(res)) {
    var <- res[i, 1]
    df_tmp <- df %>%  dplyr::select(dplyr::all_of(var)) %>% na.omit()
    value_tmp <- df_tmp[[var]]
    n <- length(value_tmp)
    if (n < 3) stop("Too small sample size")
    if (n >= 3 & n < 5000) {
      test <- shapiro.test(value_tmp)
      res[i, 2] <- test$p.value
      res[i, 3] <- ifelse(test$p.value > 0.05, "Normal", "Non-normal")
      res[i, 4] <- "Shapiro-Wilk test"
    } else { # (i.e., n >= 5000)
      test <- nortest::lillie.test(value_tmp) # ref: https://zhuanlan.zhihu.com/p/332779277
      res[i, 2] <- test$p.value
      res[i, 3] <- ifelse(test$p.value > 0.05, "Normal", "Non-normal")
      res[i, 4] <- "Kolmogorov-Smirnov test"
    }
  }

  # Output messages for easy copying
  normal_var <- res %>% dplyr::filter(normality == "Normal") %>% dplyr::pull(var)
  non_normal_var <- res %>% dplyr::filter(normality == "Non-normal") %>% dplyr::pull(var)
  cli::cli_alert_success("We detected {length(normal_var)} normal var{?s} : {paste(normal_var, collapse = ', ')}")
  cli::cli_alert_success("We detected {length(non_normal_var)} non-normal var{?s} : {paste(non_normal_var, collapse = ', ')}")
  return(res)
}

#' Step 2: Generate Table 1 for Cohort Study
#'
#' This function generates a descriptive table for clinical cohort studies, with options to stratify by group and perform group comparisons.
#'
#' @param df A data frame or tibble containing the clinical variables.
#' @param strata Character, the column name to stratify by. Use "none" to disable stratification.
#' @param var_all Vector. All var to be included in table 1
#' @param var_non Vector. Non-normal continuous variables. (Process in `print`)
#' @param var_cat Vector. Categorical variables.
#' @param var_exact Vector. var_cat that needs fisher exact test. (Process in `print`)
#' @param showAllLevels Logical, if TRUE, all levels of factor variables are shown.
#' @param compare_test Logical, if TRUE, performs group comparison (only works when `strata` is not "none").
#' @param includeNA Logical, if TRUE, includes NAs as a level in categorical variables.
#' @param verbose Logical, if TRUE, prints the table.
#'
#' @return Printed table1
#' @importFrom tableone CreateTableOne
#' @importFrom cli cli_alert_success cli_alert_info
#' @export
leo.table1.step2 <- function(df, var_all, var_cat, var_non, strata = NULL, var_exact = NULL,
                             compare_test = F, includeNA = F, showAllLevels = F, verbose = T) {
  cli::cli_alert_info("Generating tableone object...")

  # Build arguments for CreateTableOne
  tbl_args <- list(
    vars       = var_all,
    data       = df,
    factorVars = var_cat,       # var that is num，but processed as cat
    includeNA  = includeNA,     # if TRUE, consider NA as a level
    test       = compare_test,  # if TRUE, perform group comparisons
    addOverall = TRUE
  )
  if (!is.null(strata)) tbl_args$strata <- strata # Conditionally add stratification

  # Invoke CreateTableOne with dynamic args
  table <- do.call(tableone::CreateTableOne, tbl_args)
  # table <- tableone::CreateTableOne(
  #   vars = var_all,
  #   strata = if (strata == "none") NULL else strata,
  #   data = df,
  #   factorVars = var_cat,
  #   includeNA = includeNA,
  #   test = compare_test,
  #   addOverall = TRUE
  # )

  # Print tableone object --- this is where to adjust the non-normal vars and exact vars
  cli::cli_alert_info("Printing table 1...")
  table_print <- print(
    table,
    nonnormal = var_non,
    exact = var_exact,
    catDigits = 1,  # 连续变量1位小数
    contDigits = 2, # 分类变量2位小数
    pDigits = 3,    # P value 3位小数
    # other settings
    test = compare_test,   # Show group comparisons
    missing = includeNA,   # Show missing values
    showAllLevels = showAllLevels,  # Show all levels of var_cat
    quote = FALSE,         # Whether to show everything in quotes.
    noSpaces = TRUE,       # 删除用于对齐的空格
    printToggle = verbose  # Whether to print the output
  )
  cli::cli_alert_success("Table 1 has been generated successfully!")

  return(table_print)
}

#' Save Table 1 to Word or CSV
#'
#' This function saves the Table 1 output as a Word or CSV file, depending on the file extension.
#'
#' @param table A table object (e.g., from CreateTableOne).
#' @param filename A character string indicating the output file name (e.g., "table1.csv" or "table1.docx").
#' @importFrom flextable save_as_docx align
#' @importFrom cli cli_alert_success
#' @importFrom xtable xtable
#' @importFrom tools file_ext
#' @export
leo.table1.save <- function(table, filename) {

  clean_names <- function(tbl) {
    # Clean row names: replace "X." with tab length (8 spaces)
    rownames(tbl) <- gsub("^X\\.", paste(rep("\u200B", 2), collapse = ""), rownames(tbl))   # Replace X. with 8 spaces
    rownames(tbl) <- gsub("\\.{2,}", ".", rownames(tbl))  # Clean excessive dots in row names
    return(tbl)
  }

  file_extension <- tools::file_ext(filename)
  if (file_extension == "csv") {
    write.csv(table, file = filename)
  } else if (file_extension == "docx") {
    table_xtable <- xtable::xtable(table)
    table_xtable <- clean_names(table_xtable)

    table_flextable <- flextable::align(flextable::as_flextable(table_xtable), j = 1, align = "left")
    flextable::save_as_docx(table_flextable, path = filename)
  } else {
    stop("Unsupported file extension. Please use .csv or .docx.")
  }

  return(cli::cli_alert_success("Table 1 has been saved to {filename}."))
}
