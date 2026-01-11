#' Read UK Biobank HLA v2 header (ukb_hla_v2.txt)
#'
#' Downloads and parses the HLA v2 header file used to decode UK Biobank
#' Field 22182 (HLA imputation v2; best-guess probabilities).
#'
#' The returned character vector is the column order (length typically 362)
#' corresponding to the comma-separated values in field 22182.
#'
#' @details
#' The function tries `readLines()` first. If that fails, it falls back to
#' a system `wget` call (useful on clusters). If both fail, it errors.
#'
#' @return A character vector of allele column names (e.g., "A_101", "B_2705", ...).
#' @importFrom leo.basic leo_log
#' @examples
#' ukb_hla_header <- ukb_hla_header()
#' head(ukb_hla_header)
#' @export
ukb_hla_header <- function() {
  leo.basic::leo_log("Loading UKB HLA header; description source:",
                     "https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=2182")
  url <- "https://biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_hla_v2.txt"
  line <- tryCatch(
    readLines(url, n = 1, warn = FALSE)[1],
    error = function(e) {
      system(sprintf("wget -qO- %s", shQuote(url)), intern = TRUE)[1]
    }
  )
  fields <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(fields) != 362) {
    stop(sprintf("Unexpected header length: got %d fields; expected 362.", length(fields)))
  }
  return(fields)
}

#' HLA Genotyping in UKB
#'
#' Parse UK Biobank HLA v2 compound field (22182) and derive per-locus genotypes
#'
#' @param df A data.frame/tibble with `eid` and a compound column (default `p22182`).
#' @param header Allele column order, e.g. `ukb_hla_header()` (character vector).
#' @param col Compound column name. Default "p22182".
#' @param q_threshold Posterior threshold (UKB commonly uses 0.7). Default 0.7.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{calls}: data.table with per-allele calls (eid, locus, allele, allele_pretty, allele_copies, allele_q).
#'   \item \code{genotype}: data.table with per-locus genotype strings (eid, locus, genotype).
#' }
#'
#' @importFrom leo.basic leo_log
#' @importFrom dplyr mutate filter select group_by summarise %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split_fixed str_remove str_pad
#' @importFrom data.table as.data.table setorder fifelse
#' @examples
#' \dontrun{
#' library(tidyverse)
#' # x <- data.table::fread("~/Project/UKB/data/HLA.csv")
#' # x <- data.table::fread("HLA.csv")         # eid + p22182
#' header <- ukb_hla_header()
#' res <- ukb_hla_typing(x, header)
#' res$genotype %>% dplyr::filter(locus == "B") %>% dplyr::slice_head(n = 5)
#' }
#' @export
ukb_hla_typing <- function(df, header, col = "p22182", q_threshold = 0.7) {
  stopifnot(is.data.frame(df), "eid" %in% names(df), col %in% names(df), is.character(header))
  leo.basic::leo_log("Typing UKB HLA from", col, "with q_threshold =", q_threshold)

  # split compound field into columns (only need the p22182 column)
  df_hla <- df[[col]] %>%
    stringr::str_split_fixed(",", n = length(header)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>% setNames(header) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))

  # combine with eid
  df_hla <- df_hla %>% dplyr::bind_cols(df %>% dplyr::select(eid), .)

  # longer dataframe
  calls <- df_hla %>%
    tidyr::pivot_longer(-eid, names_to = "allele", values_to = "w") %>%
    dplyr::filter(!is.na(w), w > 0) %>%
    dplyr::mutate(locus = stringr::str_remove(allele, "_.*$"),
                  code = stringr::str_remove(allele, "^.*_"),
                  code = dplyr::if_else(nchar(code) == 3,
                                        stringr::str_pad(code, 4, pad = "0"),
                                        code))

  # Infer homo/hetero from the number of nonzero alleles per eid×locus
  DT <- data.table::as.data.table(calls)
  DT[, n_alleles_in_locus := .N, by = .(eid, locus)]
  DT[, is_homozygous := (n_alleles_in_locus == 1L)]
  DT[, allele_copies := data.table::fifelse(is_homozygous, 2L, 1L)]
  DT[, allele_q      := data.table::fifelse(is_homozygous, w / 2, w)]
  DT[, allele_pretty := data.table::fifelse(
    code == "9901" & locus %in% c("DRB3", "DRB4", "DRB5"),
    "ABSENT",
    paste0("HLA-", locus, "*", substr(code, 1, 2), ":", substr(code, 3, 4))
  )]
  calls_out <- DT[allele_q >= q_threshold,
                  .(eid, locus, allele, allele_pretty, allele_copies, allele_q)]

  # datatable of genotype per eid×locus
  data.table::setorder(calls_out, eid, locus, allele_pretty)
  genotype_out <- calls_out[
    , .(genotype = paste(rep(allele_pretty, allele_copies), collapse = "/")),
    by = .(eid, locus)
  ]

  return(list(calls = calls_out, genotype = genotype_out))
}
