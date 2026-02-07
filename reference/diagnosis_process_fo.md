# Process first occurrence (FO) data in UKB —-

Processes UK Biobank first occurrence data for survival analysis,
handling abnormal dates, censoring, and generating status/time
variables.

## Usage

``` r
diagnosis_process_fo(
  data,
  fo_date_col,
  censored = "2025-01-01",
  base_line = "p53_i0",
  filter_NA = T
)
```

## Arguments

- data:

  Dataframe with UKB baseline / FO / death dates.

- fo_date_col:

  Column name (string) of first occurrence date to use.

- censored:

  Censoring date (character, "YYYY-MM-DD").

- filter_NA:

  Logical, whether to filter out rows with NA in the resulting
  status/time columns. Default TRUE.

## Value

Dataframe with eid, disease_prefix_status, disease_prefix_time.

## Examples

``` r
# Simulated UKB data
ukb_df <- tibble::tibble(
  eid = 1:6,
  p53_i0 = as.Date(c("2006-01-01", "2007-05-15", "2008-03-20", "2009-07-30", "2010-11-11", "2011-12-25")),
  p191 = as.Date(c("2020-01-01", "2021-06-15", NA, "2022-03-20", "2023-07-30", "2024-11-11")),
  p40000_i0 = as.Date(c(NA, "2022-05-05", "2023-08-08", NA, "2024-09-09", NA)),
  p40000_i1 = as.Date(c("2023-01-01", NA, NA, "2024-04-04", NA, "2025-12-12")),
  p20002_0_0 = as.Date(c("2015-01-01", "2016-02-02", "1900-01-01", NA, "2018-03-03", "2019-04-04"))
)
diagnosis_process_fo(ukb_df, "p20002_0_0", censored = "2025-01-01") # This is the results.
#> ℹ [15:27:59] diagnosis_process_fo: processing [p20002_0_0] using [2025-01-01] as censored date.
#> ℹ [15:27:59] diagnosis_process_fo: combining death dates p40000_i0 and p40000_i1.
#> ℹ [15:27:59] diagnosis_process_fo: detected <1> abnormal FO dates in [p20002_0_0], recoding to NA.
#> ℹ [15:27:59] diagnosis_process_fo: filtered NA - rows before: 6, after: 5.
#> ✔ [15:27:59] diagnosis_process_fo: created [p20002_status] and [p20002_time] - ALL DONE.
#> # A tibble: 5 × 3
#>     eid p20002_status p20002_time
#>   <int>         <int>       <dbl>
#> 1     1             1        9.01
#> 2     2             1        8.73
#> 3     4             0       12.6 
#> 4     5             1        7.31
#> 5     6             1        7.28
```
