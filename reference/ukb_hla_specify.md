# Specify carrier status for given HLA allele(s) in UK Biobank **\[stable\]**

Extract per-individual carrier status (0/1/2 copies) for user-specified
HLA allele queries from the output of
[`ukb_hla_typing()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_typing.md).

## Usage

``` r
ukb_hla_specify(res, specify_alleles = "B27")
```

## Arguments

- res:

  Output of
  [`ukb_hla_typing()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_typing.md),
  i.e. `list(calls=..., genotype=...)`.

- specify_alleles:

  Character vector of allele queries. Multiple formats are supported:

  - **Family-level (2-digit) queries** (matches any `xx:yy` within the
    family):

    - `"B27"` or `"HLA-B27"`  
      Matches any `HLA-B*27:xx`.

    - `"B*27"` or `"B*27:"`  
      Matches any `HLA-B*27:xx`.

    - `"DQB1*06"`  
      Matches any `HLA-DQB1*06:xx`.

  - **Exact (4-digit) allele queries** (matches a specific `xx:yy`
    allele):

    - `"A*02:01"` or `"HLA-A*02:01"`  
      Matches exactly `HLA-A*02:01`.

    - `"B*27:05"` or `"HLA-B*27:05"`  
      Matches exactly `HLA-B*27:05`.

  - **UK Biobank header-style queries** (column name in Field 22182
    header):

    - `"DRB1_1501"`, `"B_2705"`, `"C_401"`  
      Converted internally to `HLA-<locus>*xx:yy` (e.g. `"B_2705"` -\>
      `"HLA-B*27:05"`). Three-digit codes are left-padded to four digits
      (e.g. `"C_401"` -\> `"HLA-C*04:01"`).

## Value

A table with one row per `eid × query`, containing:

- `eid` Individual identifier (if the eid has HLA data).

- `query` The original query string as provided in `specify_alleles`.

- `locus` Parsed locus (e.g., `"B"`, `"A"`, `"DRB1"`).

- `copies` Integer 0/1/2 indicating the total number of matched allele
  copies.

- `carrier` Logical; `TRUE` if `copies > 0`.

- `matched` Matched `allele_pretty` values (semicolon-separated); `NA`
  if non-carrier.

- `max_q` Maximum posterior (per-allele) among matched calls; `NA` if
  non-carrier.

, which is recommend to use code in the examples to convert to a wide
format for downstream analysis.

## Details

This function searches within `res$calls` (already filtered by the
posterior threshold in
[`ukb_hla_typing()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_typing.md))
and returns, for each query, the copy number carried by each individual.

The returned copy number (`copies`) is derived by summing
`allele_copies` among the matched calls in `res$calls`. Since
`res$calls` has already been filtered by the posterior threshold in
[`ukb_hla_typing()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_typing.md),
this function reports carriers based on those retained (high-confidence)
allele calls.

## Examples

``` r
# `res` are from ukb_hla_typing()
# 1) Single query: family-level (2-digit)
res_s1 <- ukb_hla_specify(res, "B27")
#> ℹ [15:28:08] ukb_hla_specify: extracting B27
#> Error: object 'res' not found
# 2) Multiple queries: mixed formats
res_s2 <- ukb_hla_specify(res, c("B27", "A*02:01", "DRB1_1501"))
#> ℹ [15:28:08] ukb_hla_specify: extracting B27, A*02:01, DRB1_1501
#> Error: object 'res' not found
# 3) Queries with optional HLA- prefix and family-level with "*"
res_s3 <- ukb_hla_specify(res, c("HLA-B27", "DQB1*06"))
#> ℹ [15:28:08] ukb_hla_specify: extracting HLA-B27, DQB1*06
#> Error: object 'res' not found

### Now normally we only need to know if one is a carrier for each allele
# Create a downstream-friendly wide table:
# eid, <query1>, <query1>_carrier, <query2>, <query2>_carrier, ...
carrier_wide <- res_s2 %>%
  transmute(eid, query,
            copies = as.integer(copies),
            carrier01 = as.integer(carrier)) %>%
  pivot_longer(c(copies, carrier01), names_to = "stat", values_to = "value") %>%
  mutate(key = if_else(stat == "copies", query, paste0(query, "_carrier"))) %>%
  select(eid, key, value) %>%
  pivot_wider(names_from = key, values_from = value,
              values_fill = list(value = 0L))
#> Error in res_s2 %>% transmute(eid, query, copies = as.integer(copies),     carrier01 = as.integer(carrier)) %>% pivot_longer(c(copies,     carrier01), names_to = "stat", values_to = "value") %>% mutate(key = if_else(stat ==     "copies", query, paste0(query, "_carrier"))) %>% select(eid,     key, value) %>% pivot_wider(names_from = key, values_from = value,     values_fill = list(value = 0L)): could not find function "%>%"
```
