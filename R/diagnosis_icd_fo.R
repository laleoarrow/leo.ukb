#' Process first occurrence (FO) data in UKB ----
#'
#' Processes UK Biobank first occurrence data for survival analysis, handling abnormal dates,
#' censoring, and generating status/time variables.
#'
#' @param data Dataframe with UKB baseline / FO / death dates.
#' @param fo_date_col Column name (string) of first occurrence date to use.
#' @param censored Censoring date (character, "YYYY-MM-DD").
#' @param filter_NA Logical, whether to filter out rows with NA in the resulting status/time columns. Default TRUE.
#'
#' @return Dataframe with eid, {disease_prefix}_status, {disease_prefix}_time.
#' @export
#'
#' @examples
#' # Simulated UKB data
#' ukb_df <- tibble::tibble(
#'   eid = 1:6,
#'   p53_i0 = as.Date(c("2006-01-01", "2007-05-15", "2008-03-20", "2009-07-30", "2010-11-11", "2011-12-25")),
#'   p191 = as.Date(c("2020-01-01", "2021-06-15", NA, "2022-03-20", "2023-07-30", "2024-11-11")),
#'   p40000_i0 = as.Date(c(NA, "2022-05-05", "2023-08-08", NA, "2024-09-09", NA)),
#'   p40000_i1 = as.Date(c("2023-01-01", NA, NA, "2024-04-04", NA, "2025-12-12")),
#'   p20002_0_0 = as.Date(c("2015-01-01", "2016-02-02", "1900-01-01", NA, "2018-03-03", "2019-04-04"))
#' )
#' dignosis_process_fo(ukb_df, "p20002_0_0", censored = "2025-01-01") # This is the results.
dignosis_process_fo <- function(data, fo_date_col, censored = "2025-01-01", base_line="p53_i0", filter_NA = T) {
  # check
  required_cols <- c("eid", base_line, "p191", "p40000_i0", "p40000_i1", fo_date_col)
  missing_required <- setdiff(required_cols, names(data))
  if (length(missing_required) > 0) {
    leo.basic::leo_log("dignosis_process_fo: missing [{paste(missing_required, collapse = ', ')}] - re-check data", level = "danger");return(NULL)
  } else { leo.basic::leo_log("dignosis_process_fo: processing [{fo_date_col}] using [{censored}] as censored date.") }

  # subset and format date columns
  tmp <- data %>% dplyr::select(dplyr::all_of(required_cols))
  date_cols <- c(base_line, "p191", "p40000_i0", "p40000_i1", fo_date_col)
  tmp <- tmp %>% dplyr::mutate(dplyr::across(.cols = dplyr::all_of(date_cols),
                                             .fns  = ~ suppressWarnings(as.Date(.x, format = "%Y-%m-%d"))
                                             ))
  censored <- as.Date(censored, format = "%Y-%m-%d")  # censored as Date

  # combine death date p40000 - just in case
  leo.basic::leo_log("dignosis_process_fo: combining death dates p40000_i0 and p40000_i1.", level = "info")
  tmp <- tmp %>%
    dplyr::mutate(p40000 = dplyr::coalesce(p40000_i0, p40000_i1)) %>%
    dplyr::select(-p40000_i0, -p40000_i1)

  # build status and time column names
  disease_prefix    <- strsplit(fo_date_col, "_")[[1]][1]
  disease_status_col <- paste0(disease_prefix, "_status")
  disease_time_col   <- paste0(disease_prefix, "_time")

  # abnormal FO dates (UKB placeholders) as Date
  abnormal_dates <- as.Date(c("1900-01-01", "1901-01-01", "1902-02-02",
                              "1903-03-03", "1909-09-09", "2037-07-07"))
  tmp$is_abnormal_fo <- tmp[[fo_date_col]] %in% abnormal_dates
  n_abnormal <- sum(tmp$is_abnormal_fo, na.rm = TRUE)
  if (n_abnormal > 0) {
    leo.basic::leo_log("dignosis_process_fo: detected <{n_abnormal}> abnormal FO dates in [{fo_date_col}], recoding to NA.")
  } else {leo.basic::leo_log("dignosis_process_fo: no abnormal FO dates detected in [{fo_date_col}].")}

  tmp <- tmp %>%
    dplyr::mutate(
      # clean event date: abnormal codes -> NA
      event_date = .data[[fo_date_col]],
      event_date = dplyr::if_else(event_date %in% abnormal_dates, as.Date(NA), event_date),

      # base censoring date: min(p191, p40000, censored)
      censor_base = pmin(p191, p40000, censored, na.rm = TRUE),
      censor_base = dplyr::if_else(is.infinite(censor_base), censored, censor_base),

      # final follow-up date: event_date if occurs before censoring, else censor_base
      final_date = dplyr::if_else(!is.na(event_date) & event_date <= censor_base,
                                  event_date, censor_base),

      # time in years from baseline 
      # status: 1 if event_date exists and occurs before censoring
      tmp_time = as.numeric(difftime(final_date, .data[[base_line]], units = "days"))/365,
      tmp_status = dplyr::if_else(!is.na(event_date) & event_date <= censor_base, 1L, 0L),
      # set time/status to NA if abnormal FO
      !!disease_time_col := dplyr::if_else(is_abnormal_fo, as.numeric(NA), tmp_time),
      !!disease_status_col := dplyr::if_else(is_abnormal_fo, as.integer(NA), tmp_status)
    ) %>%
    dplyr::select(eid,
                  dplyr::all_of(disease_status_col),
                  dplyr::all_of(disease_time_col))
  # filter NA if specified
  if (filter_NA) {
    n_before <- nrow(tmp)
    tmp <- tmp %>%
      dplyr::filter(!is.na(.data[[disease_status_col]]) & !is.na(.data[[disease_time_col]]))
    n_after <- nrow(tmp)
    leo.basic::leo_log("dignosis_process_fo: filtered NA - rows before: {n_before}, after: {n_after}.", level = "info")
  }
  leo.basic::leo_log("dignosis_process_fo: created [{disease_status_col}] and [{disease_time_col}] - ALL DONE.", level = "success")
  return(tmp)
}

#' Process UKB ICD data for survival analysis
#'
#' Generate {icd_code}_status and {icd_code}_time from UKB ICD records.
#'
#' @param data A data.frame containing UKB baseline, censoring, ICD code and date columns.
#' @param icd_code A target ICD code prefix or a character vector of prefixes (e.g., "M10" or c("E10", "E11")).
#' @param icd ICD version, 10 or 9.
#' @param censored Censoring date ("YYYY-MM-DD").
#' @param icd_rank Rank of diagnosis occurrence to consider (1 = first occurrence, 2 = second occurrence, etc.).
#' @param baseline_date Baseline date column name (default "p53_i0").
#' @param filter_NA Logical, whether to filter out NA row. Default FALSE (Suitable for comorbidity check).
#'
#' @return A data.frame with eid, {icd_code}_status, {icd_code}_time.
#' @export
#'
#' @importFrom dplyr select all_of mutate across coalesce if_else rowwise ungroup
#' @importFrom rlang .data sym
#' @importFrom leo.basic leo_log
#' @examples
#' df_example <- data.frame(
#'   eid = 1:5,
#'   p53_i0 = as.Date(c("2006-01-01", "2007-06-01", "2005-03-15", "2004-10-30", "2008-05-10")),  # baseline
#'   p191 = as.Date(c("2020-01-01", "2019-12-31", "2021-06-30", "2022-07-01", "2023-09-10")),    # last follow-up
#'   p40000_i0 = as.Date(c(NA, NA, NA, "2021-06-01", NA)),                                       # death date (field 0)
#'   p40000_i1 = as.Date(c(NA, NA, NA, NA, NA)),                                                 # death date (field 1)
#'   p41270 = c("D261|E113|M102", "E113|M101", NA, "D261|M102|M103", "E113|M101|D261"),          # ICD10 codes
#'   p41280_a0 = as.Date(c("2010-05-01", "2011-07-20", NA, "2015-03-01", "2017-07-05")),         # 1st diagnosis date
#'   p41280_a1 = as.Date(c("2012-03-10", "2014-09-01", NA, "2016-06-01", "2018-08-10")),         # 2nd diagnosis date
#'   p41280_a2 = as.Date(c("2015-08-15", NA, NA, "2015-09-10", NA)),                             # 3rd diagnosis date
#'   p41280_a3 = as.Date(c(NA, NA, NA, NA, NA))                                                  # Additional diagnosis field for testing
#' )
#'
#' # Prefix matching: D261 matches only D261 (icd_rank = 1)
#' diagnosis_process_icd(df_example, icd_code = "D261", icd = 10, icd_rank = 1)
#'
#' # Prefix matching: M10 matches M101 / M102 (icd_rank = 1)
#' diagnosis_process_icd(df_example, icd_code = "M10", icd = 10, icd_rank = 1)
#'
#' # Prefix matching: E113 matches E113 only (icd_rank = 1)
#' diagnosis_process_icd(df_example, icd_code = "E113", icd = 10, icd_rank = 1)
#'
#' # Prefix matching: M10 matches M101 / M102 (icd_rank = 2)
#' diagnosis_process_icd(df_example, icd_code = "M10", icd = 10, icd_rank = 2)
#'
#' # Multiple ICD codes: c("M10", "E113") - takes global earliest date among all matches
#' diagnosis_process_icd(df_example, icd_code = c("M10", "E113"), icd = 10, icd_rank = 1)
#' # Output columns: M10_E113_status, M10_E113_time
diagnosis_process_icd <- function(data, icd_code, icd = 10, censored = "2025-01-01", icd_rank = 1, baseline_date = "p53_i0", filter_NA = FALSE) {
  icd <- as.integer(icd)
  if (!icd %in% c(9L, 10L)) stop("diagnosis_process_icd: `icd` must be 9 or 10.")

  # Resolve ICD columns by version
  icd_col <- if (icd == 10L) "p41270" else "p41271"
  date_prefix <- if (icd == 10L) "p41280_a" else "p41281_a"
  date_cols <- grep(paste0("^", date_prefix), names(data), value = T)

  # Check required columns
  required_cols <- c("eid", baseline_date, "p191", "p40000_i0", "p40000_i1", icd_col, date_cols)
  missing_required <- setdiff(required_cols, names(data))
  if (length(missing_required) > 0) {
    leo.basic::leo_log("diagnosis_process_icd: missing [{paste(missing_required, collapse = ', ')}] - re-check data", level = "danger");return(NULL)
  }
  leo.basic::leo_log("diagnosis_process_icd: processing ICD{icd} [{icd_code}] using [{censored}] as censored date.", level = "info")

  df <- data %>% dplyr::select(dplyr::all_of(required_cols))

  # Convert date columns
  leo.basic::leo_log("diagnosis_process_icd: converting date columns to <Date>.", level = "info")
  date_columns <- c(baseline_date, "p191", "p40000_i0", "p40000_i1", date_cols)
  df <- df %>% dplyr::mutate(dplyr::across(dplyr::all_of(date_columns),
                                           ~ suppressWarnings(as.Date(.x, "%Y-%m-%d"))))
  censored <- as.Date(censored, "%Y-%m-%d")

  # Merge death dates
  leo.basic::leo_log("diagnosis_process_icd: combining death dates.", level = "info")
  df <- df %>% dplyr::mutate(p40000 = dplyr::coalesce(p40000_i0, p40000_i1)) %>%
    dplyr::select(-p40000_i0, -p40000_i1)

  prefix <- paste(icd_code, collapse = "_")
  status_col <- paste0(prefix, "_status")
  time_col <- paste0(prefix, "_time")

  # Match ICD prefix(es) and extract diagnosis date by rank
  icd_pattern <- paste0("^(", paste(icd_code, collapse = "|"), ")")
  leo.basic::leo_log("diagnosis_process_icd: matching ICD{icd} pattern [{icd_pattern}].", level = "info")
  diagnosis_list <- df[[icd_col]]

  df$event_date <- sapply(seq_len(nrow(df)), function(i) {
    if (is.na(diagnosis_list[i])) return(as.Date(NA))
    diag_vec <- strsplit(diagnosis_list[i], "\\|")[[1]]
    match_pos <- which(grepl(icd_pattern, diag_vec))
    if (length(match_pos) == 0L) return(as.Date(NA))
    # filter valid positions within date_cols range
    match_pos <- match_pos[match_pos <= length(date_cols)]
    if (length(match_pos) == 0L) return(as.Date(NA))
    # collect all matched dates and sort globally
    matched_dates <- sapply(date_cols[match_pos], function(col) df[[col]][i])
    sorted_dates <- sort(matched_dates, na.last = TRUE)
    sorted_dates <- sorted_dates[!is.na(sorted_dates)]
    if (length(sorted_dates) == 0L || icd_rank > length(sorted_dates)) return(as.Date(NA))
    return(sorted_dates[icd_rank])
  }) %>% as.Date()

  # Compute survival time and status
  leo.basic::leo_log("diagnosis_process_icd: calculating survival time and status.", level = "info")
  df <- df %>%
    dplyr::mutate(censor_base = pmin(p191, p40000, censored, na.rm = T),
                  censor_base = dplyr::if_else(is.na(censor_base), censored, censor_base),
                  # censor_base = dplyr::if_else(is.infinite(censor_base), censored, censor_base),
                  final_date = dplyr::if_else(
                    !is.na(event_date) & event_date <= censor_base, event_date, censor_base
                  ),
                  !!time_col := as.numeric(difftime(final_date, !!rlang::sym(baseline_date), units = "days")) / 365,
                  !!status_col := dplyr::if_else(!is.na(event_date) & (event_date <= censor_base), 1L, 0L))

  result <- df %>% dplyr::select(dplyr::all_of(c("eid", status_col, time_col)))

  if (filter_NA) {
    n_before <- nrow(result)
    result <- result %>%
      dplyr::filter(!is.na(.data[[status_col]]) & !is.na(.data[[time_col]]))
    n_after <- nrow(result)
    leo.basic::leo_log("diagnosis_process_icd: filtered NA - rows before: {n_before}, after: {n_after}.", level = "info")
  }

  leo.basic::leo_log("diagnosis_process_icd: created [{status_col}], [{time_col}] - ALL DONE.", level = "success")

  return(result)
}

#' Process multiple ICD definitions in bulk
#'
#' Apply diagnosis_process_icd() to a list or table of ICD code sets and
#' return one status/time pair per diagnosis.
#'
#' @param data A data.frame containing UKB baseline, censoring, ICD code and date columns.
#' @param icd_map A named list (diagnosis -> vector of ICD codes) or a data.frame.
#'   If data.frame: first column is diagnosis name; remaining columns or a single
#'   "codes" column contain ICD codes (comma/space separated or one per cell).
#' @param icd ICD version, 10 or 9.
#' @param censored Censoring date ("YYYY-MM-DD").
#' @param icd_rank Rank of diagnosis occurrence to consider (1 = first occurrence).
#' @param baseline_date Baseline date column name (default "p53_i0").
#'
#' @return A data.frame with eid and {diagnosis}_status/{diagnosis}_time columns.
#' @export
diagnosis_process_icd_multi <- function(data, icd_map, icd = 10, censored = "2025-01-01", icd_rank = 1, baseline_date = "p53_i0") {
  if (is.data.frame(icd_map)) {
    dn <- as.character(icd_map[[1]])
    if (ncol(icd_map) == 2L) {
      codes_raw <- as.character(icd_map[[2]])
      code_list <- lapply(codes_raw, function(x) {
        if (is.na(x) || x == "") return(character(0))
        gsub("\\s+", "", unlist(strsplit(x, ",")))
      })
    } else {
      code_list <- apply(icd_map[, -1, drop = FALSE], 1, function(r) {
        r <- r[!is.na(r) & r != ""]
        gsub("\\s+", "", as.character(r))
      })
    }
    names(code_list) <- dn
    icd_map <- code_list
  }
  if (!is.list(icd_map) || is.null(names(icd_map))) {
    stop("diagnosis_process_icd_multi: `icd_map` must be a named list or a data.frame with diagnosis names.")
  }

  results <- list()
  for (dx in names(icd_map)) {
    codes <- unique(as.character(icd_map[[dx]]))
    codes <- codes[!is.na(codes) & codes != ""]
    if (length(codes) == 0L) next

    per_code <- lapply(codes, function(code) {
      diagnosis_process_icd(data, icd_code = code, icd = icd, censored = censored,
                           icd_rank = icd_rank, baseline_date = baseline_date)
    })
    df <- Reduce(function(x, y) dplyr::left_join(x, y, by = "eid"), per_code)
    status_cols <- paste0(codes, "_status")
    time_cols <- paste0(codes, "_time")

    status_mat <- as.matrix(df[, status_cols, drop = FALSE])
    time_mat <- as.matrix(df[, time_cols, drop = FALSE])

    dx_status <- apply(status_mat, 1, function(r) {
      if (all(is.na(r))) return(NA_integer_)
      if (any(r == 1, na.rm = TRUE)) return(1L)
      0L
    })

    dx_time <- vapply(seq_len(nrow(df)), function(i) {
      r_status <- status_mat[i, ]
      r_time <- time_mat[i, ]
      if (all(is.na(r_status))) return(NA_real_)
      if (any(r_status == 1, na.rm = TRUE)) {
        return(min(r_time[r_status == 1], na.rm = TRUE))
      }
      return(r_time[!is.na(r_time)][1])
    }, numeric(1))

    out <- df[, "eid", drop = FALSE]
    out[[paste0(dx, "_status")]] <- as.integer(dx_status)
    out[[paste0(dx, "_time")]] <- as.numeric(dx_time)
    results[[dx]] <- out
  }

  if (length(results) == 0L) return(data.frame())
  final <- Reduce(function(x, y) dplyr::left_join(x, y, by = "eid"), results)
  return(final)
}

#' Find the max date in given date columns
#' `r lifecycle::badge('stable')`
#'
#' Generate {icd_code}_status and {icd_code}_time from UKB ICD records.
#'
#' @param df A data.frame containing UKB statistics
#' @param date_columns A character vector of column names or patterns to select date columns.
#' @examples
#' ukb_sim <- tibble::tibble(
#'   eid = 1:8,
#'   p53_i0 = as.Date(c("2006-01-01", "2007-05-15", "2008-03-20", "2009-07-30",
#'                      "2010-11-11", "2011-12-25", "2005-02-14", "2004-09-09")),
#'   p191 = as.Date(c("2020-01-01", "2021-06-15", NA, "2022-03-20",
#'                    "2023-07-30", "2024-11-11", "2019-12-31", "2022-08-08")),
#'   p40000_i0 = as.Date(c(NA, "2022-05-05", "2023-08-08", NA, "2024-09-09", NA, NA, "2021-01-01")),
#'   p40000_i1 = as.Date(c("2023-01-01", NA, NA, "2024-04-04", NA, "2025-12-12", NA, NA)),
#'
#'   # self-reported first occurrence (example; includes an abnormal placeholder)
#'   p20002_0_0 = as.Date(c("2015-01-01", "2016-02-02", "1900-01-01", NA, "2018-03-03", "2019-04-04", "2020-06-06", "1901-01-01")),
#'
#'   # ICD10 codes (pipe-separated) and corresponding diagnosis dates
#'   p41270 = c(
#'     "D261|E113|M102",
#'     "E113|M101",
#'     NA,
#'     "D261|M102|M103",
#'     "E113|M101|D261",
#'     "M101|M100",
#'     "M101",
#'     "D261"
#'   ),
#'   p41280_a0 = as.Date(c("2012-05-01", "2021-07-20", NA, "2015-03-01", "2017-07-05", "2012-02-02", "2013-04-04", "2009-09-09")),
#'   p41280_a1 = as.Date(c("2012-03-10", "2014-09-01", NA, "2016-06-01", "2018-08-10", "2014-05-05", NA, NA)),
#'   p41280_a2 = as.Date(c("2015-08-15", NA, NA, "2015-09-10", NA, NA, NA, NA)),
#'   p41280_a3 = as.Date(c(NA, NA, NA, NA, NA, NA, NA, NA))
#' )
#'
#' get_max_date(ukb_sim, date_columns = c("p41280_a*"))
#' @export
get_max_date <- function(df, date_columns) {
  if (length(date_columns) == 0L) stop("`date_columns` is empty.")
  all_cols <- names(df)

  resolve_selector <- function(sel) {
    if (sel %in% all_cols) return(sel)
    if (grepl("[*?]", sel)) return(grep(utils::glob2rx(sel), all_cols, value = TRUE))
    grep(sel, all_cols, value = TRUE)
  }

  selected_cols <- unique(unlist(lapply(date_columns, resolve_selector), use.names = FALSE))
  if (length(selected_cols) == 0L) stop("No columns matched `date_columns`.")

  current_year <- as.integer(format(Sys.Date(), "%Y"))
  keep_years <- (current_year - 4L):current_year

  df_subset <- if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[, ..selected_cols]
  } else {
    df[selected_cols]
  }

  parsed <- lapply(df_subset, function(x) {
    d <- suppressWarnings(as.Date(as.character(x), format = "%Y-%m-%d"))
    y <- as.integer(format(d, "%Y"))
    d[is.na(d) | !(y %in% keep_years)] <- NA
    d
  })

  dates <- unlist(parsed, use.names = FALSE)
  if (length(dates) == 0L || all(is.na(dates))) return("None date in last 5 years")

  overall_max <- max(dates, na.rm = TRUE)
  return(format(as.Date(overall_max, origin = "1970-01-01"), "%F"))
}
