#' Convert Visual Acuity Measurements to logMAR
#'
#' Handles numeric decimals and special terms: CF, HM, LP, NLP, and Chinese synonyms.
#' @param data A data frame containing visual acuity data.
#' @param acuity_cols Character vector of column names in \code{data} to be converted.
#' @return A data frame with the specified columns converted to numeric logMAR values.
#' @importFrom dplyr mutate across all_of case_when
#' @export
#' @examples
#' va_df <- tibble::tibble(
#'   id = 1:5,
#'   bcva_od = c("1.0","0.8","CF","LP","NLP"),
#'   bcva_os = c("0.6","HM","0.3","1.2","0.1")
#' )
#' convert_to_logMAR(va_df, c("bcva_od","bcva_os"))
convert_to_logMAR <- function(data, acuity_cols, term_map = NULL) {
  if (is.null(term_map))
  term_map <- c("CF" = 2.0, "HM" = 2.3, "LP" = 2.6, "NLP" = 2.9,
                "数指" = 2.0, "手动" = 2.3, "HN" = 2.3,
                "光感" = 2.6, "无光感" = 2.9, "义眼" = NA_real_)

  data %>% dplyr::mutate(dplyr::across(dplyr::all_of(acuity_cols), ~{
    s <- trimws(toupper(as.character(.x)))
    s <- gsub("[，]", ",", s); s <- gsub("[。]", ".", s); s <- gsub("’|'|\\s", "", s)
    s <- gsub("^CF.*", "CF", s, ignore.case = T)
    s <- gsub("^HM.*", "HM", s, ignore.case = T)
    nx <- suppressWarnings(as.numeric(gsub(",", ".", s)))

    dplyr::case_when(
      s %in% names(term_map) ~ as.numeric(term_map[s]),
      grepl("^[0-9]+([.,][0-9]+)?$", s) & is.finite(nx) & nx > 0 & nx <= 2 ~ -log10(nx),
      TRUE ~ NA_real_
    )
  }))
}

#' Convert IOP Measurements
#'
#' This function converts intraocular pressure (IOP) measurements by replacing "无" with NA and converting the rest to numeric.
#' @param data A data frame containing IOP data.
#' @param iop_cols A character vector of IOP column names in \code{data} to be converted.
#' @return A data frame with the specified IOP columns converted to numeric, with "无" as NA.
#' @export
#' @examples
#' # Simulated IOP data，含单位与缺失写法
#' iop_df <- tibble::tibble(
#'   id = 1:6,
#'   iop_od = c("15","无","18mmHg",">30","-","N/A"),
#'   iop_os = c("16","未测","14.5","12mmHg","NA","20")
#' )
#' convert_iop(iop_df, c("iop_od","iop_os"))
convert_iop <- function(data, iop_cols) {
  na_words <- c("", "-", "NA", "N/A", "无", "未测")
  data %>% mutate(across(all_of(iop_cols), ~{
    x <- trimws(as.character(.x))
    x[toupper(x) %in% toupper(na_words)] <- NA
    suppressWarnings(as.numeric(readr::parse_number(x)))
  }))
}

#' Combine Measurements from Two Eyes
#'
#' Combines two numeric vectors representing measurements from two eyes into a single vector using a specified method (mean, max, or min).
#' Handles NA values appropriately.
#' @param vector1 Numeric vector for the first eye.
#' @param vector2 Numeric vector for the second eye.
#' @param method Method to combine the two vectors: "mean", "max", or "min". Default is "mean".
#' @return A numeric vector with combined values.
#' @export
combine_two_eye <- function(vector1, vector2, method = c("mean", "max", "min")) {
  method <- match.arg(method)
  if (length(vector1) != length(vector2)) stop("Input vectors must have the same length.")

  result <- vector("numeric", length(vector1))
  for (i in seq_along(vector1)) {
    v1 <- vector1[i]
    v2 <- vector2[i]
    if (is.na(v1) && is.na(v2)) {
      result[i] <- NA
    } else if (is.na(v1)) {
      result[i] <- v2
    } else if (is.na(v2)) {
      result[i] <- v1
    } else {
      result[i] <- switch(method,
                          mean = mean(c(v1, v2), na.rm = TRUE),
                          max = max(v1, v2, na.rm = TRUE),
                          min = min(v1, v2, na.rm = TRUE)
      )}
  }
  return(result)
}
