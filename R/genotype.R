#' ukb_hla_header
#'
#' Read UKB HLA v2 header line.
#'
#' This function streams the UK Biobank HLA v2 header file, which is
#' described at: https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=2182
#'
#' @return A data.frame with one column \code{name}.
#' @importFrom leo.basic leo_log
#' @examples
#' ukb_hla_header <- ukb_hla_header()
#' head(ukb_hla_header)
#' @export
ukb_hla_header <- function() {
  leo.basic::leo_log("Loading UKB HLA header; description source:",
                     "https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=2182")
  url <- "https://biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_hla_v2.txt"
  line <- system(sprintf("wget -qO- %s", shQuote(url)), intern = TRUE)[1]
  fields <- strsplit(trimws(line), "\\s+")[[1]]
  data.frame(hla_order = fields, stringsAsFactors = FALSE)
}

#' ukb_hla_typing
#'
#' Parse UKB HLA v2 compound field (22182) and derive per-locus genotypes.
#'
#' @param df A data.frame/tibble with `eid` and a compound column (default `p22182`).
#' @param header Allele column order, e.g. `ukb_hla_header()$hla_order`.
#' @param col Compound column name. Default "p22182".
#' @param q_threshold Posterior threshold (UKB commonly uses 0.7). Default 0.7.
#'
#' @return A list with `calls` (long table) and `genotype` (eid x locus).
#' @importFrom leo.basic leo_log
#' @importFrom dplyr mutate filter select group_by summarise %>%
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split_fixed str_remove str_pad
#' @examples
#' \dontrun{
#' library(tidyverse)
#' # x <- data.table::fread("HLA.csv")         # eid + p22182
#' header <- ukb_hla_header()$hla_order
#' res <- ukb_hla_typing(x, header)
#' res$genotype %>% dplyr::filter(locus == "A") %>% dplyr::slice_head(n = 5)
#' }
#' @export
ukb_hla_typing <- function(df, header, col = "p22182", q_threshold = 0.7) {
  stopifnot(is.data.frame(df), "eid" %in% names(df), col %in% names(df), is.character(header))
  leo.basic::leo_log("Typing UKB HLA from", col, "with q_threshold =", q_threshold)

  df_hla <- df[[col]] %>%
    stringr::str_split_fixed(",", n = length(header)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    setNames(header) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
    dplyr::bind_cols(df %>% dplyr::select(eid), .)

  calls <- df_hla %>%
    tidyr::pivot_longer(-eid, names_to = "allele", values_to = "w") %>%
    dplyr::filter(!is.na(w), w > 0) %>%
    dplyr::mutate(
      locus = stringr::str_remove(allele, "_.*$"),
      code = stringr::str_remove(allele, "^.*_"),
      code = dplyr::if_else(nchar(code) == 3, stringr::str_pad(code, 4, pad = "0"), code),
      copies = dplyr::if_else(w > 1, 2L, 1L),
      q = dplyr::if_else(w > 1, w / 2, w),
      allele_pretty = dplyr::if_else(code == "9901" & locus %in% c("DRB3", "DRB4", "DRB5"),
                                     "ABSENT",
                                     paste0("HLA-", locus, "*", substr(code, 1, 2), ":", substr(code, 3, 4))
                                     )
    ) %>%
    dplyr::filter(q >= q_threshold) %>%
    dplyr::select(eid, locus, allele, allele_pretty, copies, q)

  genotype <- calls %>%
    dplyr::group_by(eid, locus) %>%
    dplyr::summarise(genotype = paste(rep(allele_pretty, copies), collapse = "/"), .groups = "drop")

  list(calls = calls, genotype = genotype)
}
# Example usage:
library(tidyverse)
# x <- data.table::fread("HLA.csv")         # eid + p22182
header <- ukb_hla_header()$hla_order
res <- ukb_hla_typing(x, header)
res$genotype %>% dplyr::filter(locus == "A") %>% dplyr::slice_head(n = 5)
