#' Read UK Biobank HLA v2 header (ukb_hla_v2.txt)
#' `r lifecycle::badge('stable')`
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
#' header <- ukb_hla_header()
#' head(header)
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
#' `r lifecycle::badge('stable')`
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

#' Specify carrier status for given HLA allele(s) in UK Biobank
#' `r lifecycle::badge('stable')`
#'
#' Extract per-individual carrier status (0/1/2 copies) for user-specified HLA allele queries
#' from the output of [ukb_hla_typing()].
#'
#' This function searches within \code{res$calls} (already filtered by the posterior threshold
#' in [ukb_hla_typing()]) and returns, for each query, the copy number carried by each individual.
#'
#' @param res Output of [ukb_hla_typing()], i.e. \code{list(calls=..., genotype=...)}.
#'
#' @param specify_alleles Character vector of allele queries. Multiple formats are supported:
#' \itemize{
#'   \item \strong{Family-level (2-digit) queries} (matches any \code{xx:yy} within the family):
#'     \itemize{
#'       \item \code{"B27"} or \code{"HLA-B27"} \cr
#'             Matches any \code{HLA-B*27:xx}.
#'       \item \code{"B*27"} or \code{"B*27:"} \cr
#'             Matches any \code{HLA-B*27:xx}.
#'       \item \code{"DQB1*06"} \cr
#'             Matches any \code{HLA-DQB1*06:xx}.
#'     }
#'
#'   \item \strong{Exact (4-digit) allele queries} (matches a specific \code{xx:yy} allele):
#'     \itemize{
#'       \item \code{"A*02:01"} or \code{"HLA-A*02:01"} \cr
#'             Matches exactly \code{HLA-A*02:01}.
#'       \item \code{"B*27:05"} or \code{"HLA-B*27:05"} \cr
#'             Matches exactly \code{HLA-B*27:05}.
#'     }
#'
#'   \item \strong{UK Biobank header-style queries} (column name in Field 22182 header):
#'     \itemize{
#'       \item \code{"DRB1_1501"}, \code{"B_2705"}, \code{"C_401"} \cr
#'             Converted internally to \code{HLA-<locus>*xx:yy} (e.g. \code{"B_2705"} -> \code{"HLA-B*27:05"}).
#'             Three-digit codes are left-padded to four digits (e.g. \code{"C_401"} -> \code{"HLA-C*04:01"}).
#'     }
#' }
#'
#' @details
#' The returned copy number (\code{copies}) is derived by summing \code{allele_copies} among the
#' matched calls in \code{res$calls}. Since \code{res$calls} has already been filtered by the
#' posterior threshold in [ukb_hla_typing()], this function reports carriers based on those
#' retained (high-confidence) allele calls.
#'
#' @return A table with one row per \code{eid × query}, containing:
#' \itemize{
#'   \item \code{eid} Individual identifier (if the eid has HLA data).
#'   \item \code{query} The original query string as provided in \code{specify_alleles}.
#'   \item \code{locus} Parsed locus (e.g., \code{"B"}, \code{"A"}, \code{"DRB1"}).
#'   \item \code{copies} Integer 0/1/2 indicating the total number of matched allele copies.
#'   \item \code{carrier} Logical; \code{TRUE} if \code{copies > 0}.
#'   \item \code{matched} Matched \code{allele_pretty} values (semicolon-separated); \code{NA} if non-carrier.
#'   \item \code{max_q} Maximum posterior (per-allele) among matched calls; \code{NA} if non-carrier.
#' }
#' , which is recommend to use code in the examples to convert to a wide format for downstream analysis.
#'
#' @importFrom leo.basic leo_log
#' @importFrom dplyr %>% distinct filter mutate summarise left_join bind_rows arrange group_by
#' @export
#' @examples
#' # `res` are from ukb_hla_typing()
#' # 1) Single query: family-level (2-digit)
#' res_s1 <- ukb_hla_specify(res, "B27")
#' # 2) Multiple queries: mixed formats
#' res_s2 <- ukb_hla_specify(res, c("B27", "A*02:01", "DRB1_1501"))
#' # 3) Queries with optional HLA- prefix and family-level with "*"
#' res_s3 <- ukb_hla_specify(res, c("HLA-B27", "DQB1*06"))
#'
#' ### Now normally we only need to know if one is a carrier for each allele
#' # Create a downstream-friendly wide table:
#' # eid, <query1>, <query1>_carrier, <query2>, <query2>_carrier, ...
#' carrier_wide <- res_s2 %>%
#'   transmute(eid, query,
#'             copies = as.integer(copies),
#'             carrier01 = as.integer(carrier)) %>%
#'   pivot_longer(c(copies, carrier01), names_to = "stat", values_to = "value") %>%
#'   mutate(key = if_else(stat == "copies", query, paste0(query, "_carrier"))) %>%
#'   select(eid, key, value) %>%
#'   pivot_wider(names_from = key, values_from = value,
#'               values_fill = list(value = 0L))
ukb_hla_specify <- function(res, specify_alleles = "B27") {
  leo.basic::leo_log("ukb_hla_specify: extracting {paste(specify_alleles, collapse = ', ')}")

  calls_df <- res$calls
  eid_df <- res$genotype %>% dplyr::distinct(eid)
  specify_alleles <- unique(as.character(specify_alleles))

  parse_query <- function(query) {
    # Return: query, locus, type, prefix, target
    q <- trimws(query); q <- sub("^HLA-", "", q, ignore.case = T)

    if (grepl("_", q, fixed = T)) {
      locus <- sub("_.*$", "", q)
      code <- sub("^.*_", "", q)
      code <- sprintf("%04d", as.integer(code))
      target <- paste0("HLA-", locus, "*", substr(code, 1, 2), ":", substr(code, 3, 4))
      return(list(query = query, locus = locus, type = "exact", prefix = NA_character_, target = target))
    }

    if (grepl("\\*", q)) {
      locus <- sub("\\*.*$", "", q)
      rest <- sub("^.*\\*", "", q)
      if (grepl("^\\d{2}$", rest) || grepl("^\\d{2}:$", rest)) {
        fam <- sub(":$", "", rest)
        return(list(query = query, locus = locus, type = "family",
                    prefix = paste0("HLA-", locus, "*", fam, ":"), target = NA_character_))
      }
      if (grepl("^\\d{2}:\\d{2}$", rest)) {
        return(list(query = query, locus = locus, type = "exact",
                    prefix = NA_character_, target = paste0("HLA-", locus, "*", rest)))
      }
    }

    if (grepl("^([A-Za-z0-9]+)(\\d{2})$", q)) {
      locus <- sub("\\d{2}$", "", q)
      fam <- sub("^.*?(\\d{2})$", "\\1", q)
      return(list(query = query, locus = locus, type = "family",
                  prefix = paste0("HLA-", locus, "*", fam, ":"), target = NA_character_))
    }

    stop("ukb_hla_specify: cannot parse query: ", query)
  }

  result_list <- vector("list", length(specify_alleles))

  for (i in seq_along(specify_alleles)) {
    spec <- parse_query(specify_alleles[i])

    matched_df <- if (spec$type == "exact") {
      calls_df %>% dplyr::filter(locus == spec$locus, allele_pretty == spec$target)
    } else {
      calls_df %>% dplyr::filter(locus == spec$locus, startsWith(allele_pretty, spec$prefix))
    }

    summary_df <- matched_df %>%
      dplyr::group_by(eid) %>%
      dplyr::summarise(copies = sum(allele_copies), max_q = max(allele_q),
                       matched = paste(sort(unique(allele_pretty)), collapse = ";"),
                       .groups = "drop")

    result_list[[i]] <- eid_df %>%
      dplyr::left_join(summary_df, by = "eid") %>%
      dplyr::mutate(query = spec$query, locus = spec$locus,
                    copies = dplyr::if_else(is.na(copies), 0L, as.integer(copies)),
                    carrier = copies > 0L,
                    matched = dplyr::if_else(carrier, matched, NA_character_))
  }

  out_df <- dplyr::bind_rows(result_list) %>% dplyr::arrange(query, eid)
  leo.basic::leo_log("ukb_hla_specify: done")
  return(out_df)
}
