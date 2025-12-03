#' Process first occurrence (FO) data in UKB ----
#'
#' @param data Dataframe with UKB baseline / FO / death dates.
#' @param fo_date_col Column name (string) of first occurrence date to use.
#' @param censored Censoring date (character, "YYYY-MM-DD").
#'
#' @return Dataframe with eid, {disease_prefix}_status, {disease_prefix}_time.
#' @export
#'
#' @examples
data_process_fo <- function(data, fo_date_col, censored = "2025-01-01") {
  # check required columns
  required_cols <- c("eid", "p53_i0", "p191", "p40000_i0", "p40000_i1", fo_date_col)
  missing_required <- setdiff(required_cols, names(data))
  if (length(missing_required) > 0) {
    leo.basic::leo_log("data_process_fo: missing [{paste(missing_required, collapse = ', ')}] - re-check data", level = "danger")
    return(NULL)
  } else {
    leo.basic::leo_log("data_process_fo: processing [{fo_date_col}] using [{censored}] as censored date.")
  }

  # subset and format date columns
  leo.basic::leo_log("data_process_fo: converting date columns to <Date>.", level = "info")
  tmp <- data %>% dplyr::select(dplyr::all_of(required_cols))
  date_cols <- c("p53_i0", "p191", "p40000_i0", "p40000_i1", fo_date_col)
  tmp <- tmp %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(date_cols),
        .fns  = ~ suppressWarnings(as.Date(.x, format = "%Y-%m-%d"))
      )
    )
  censored <- as.Date(censored, format = "%Y-%m-%d")  # censored as Date

  # combine death date p40000 - just in case
  leo.basic::leo_log("data_process_fo: combining death dates p40000_i0 and p40000_i1.", level = "info")
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

  # check abnormal FO dates before mutate
  abnormal_mask <- tmp[[fo_date_col]] %in% abnormal_dates
  n_abnormal <- sum(abnormal_mask, na.rm = TRUE)
  if (n_abnormal > 0) {
    leo.basic::leo_log("data_process_fo: detected {n_abnormal} abnormal FO dates in [{fo_date_col}], recoding to NA.")
  } else {
    leo.basic::leo_log("data_process_fo: no abnormal FO dates detected in [{fo_date_col}].")
  }

  tmp <- tmp %>%
    dplyr::mutate(
      # clean event date: abnormal codes -> NA
      event_date = .data[[fo_date_col]],
      event_date = dplyr::if_else(event_date %in% abnormal_dates,
                                  as.Date(NA), event_date),

      # base censoring date: min(p191, p40000, censored)
      censor_base = pmin(p191, p40000, censored, na.rm = TRUE),
      censor_base = dplyr::if_else(is.infinite(censor_base), censored, censor_base),

      # final follow-up date: event_date if occurs before censoring, else censor_base
      final_date = dplyr::if_else(!is.na(event_date) & event_date <= censor_base,
                                  event_date, censor_base),

      # time in years from baseline (p53_i0)
      !!disease_time_col := as.numeric(difftime(final_date, p53_i0, units = "days")) / 365,

      # status: 1 if event_date exists and occurs before censoring
      !!disease_status_col := dplyr::if_else(!is.na(event_date) & event_date <= censor_base,
                                             1L, 0L)
    ) %>%
    dplyr::select(eid,
                  dplyr::all_of(disease_status_col),
                  dplyr::all_of(disease_time_col))

  leo.basic::leo_log("data_process_fo: created [{disease_status_col}] and [{disease_time_col}] - ALL DONE.", level = "success")
  return(tmp)
}
