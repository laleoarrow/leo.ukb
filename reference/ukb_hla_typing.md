# HLA Genotyping in UKB **\[stable\]**

Parse UK Biobank HLA v2 compound field (22182) and derive per-locus
genotypes

## Usage

``` r
ukb_hla_typing(df, header, col = "p22182", q_threshold = 0.7)
```

## Arguments

- df:

  A data.frame/tibble with `eid` and a compound column (default
  `p22182`).

- header:

  Allele column order, e.g.
  [`ukb_hla_header()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_header.md)
  (character vector).

- col:

  Compound column name. Default "p22182".

- q_threshold:

  Posterior threshold (UKB commonly uses 0.7). Default 0.7.

## Value

A list with:

- `calls`: data.table with per-allele calls (eid, locus, allele,
  allele_pretty, allele_copies, allele_q).

- `genotype`: data.table with per-locus genotype strings (eid, locus,
  genotype).

## Examples

``` r
if (FALSE) { # \dontrun{
library(tidyverse)
# x <- data.table::fread("~/Project/UKB/data/HLA.csv")
# x <- data.table::fread("HLA.csv")         # eid + p22182
header <- ukb_hla_header()
res <- ukb_hla_typing(x, header)
res$genotype %>% dplyr::filter(locus == "B") %>% dplyr::slice_head(n = 5)
} # }
```
