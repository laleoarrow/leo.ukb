# Find the max date in given date columns **\[stable\]**

Generate icd_code_status and icd_code_time from UKB ICD records.

## Usage

``` r
get_max_date(df, date_columns)
```

## Arguments

- df:

  A data.frame containing UKB statistics

- date_columns:

  A character vector of column names or patterns to select date columns.

## Examples

``` r
ukb_sim <- tibble::tibble(
  eid = 1:8,
  p53_i0 = as.Date(c("2006-01-01", "2007-05-15", "2008-03-20", "2009-07-30",
                     "2010-11-11", "2011-12-25", "2005-02-14", "2004-09-09")),
  p191 = as.Date(c("2020-01-01", "2021-06-15", NA, "2022-03-20",
                   "2023-07-30", "2024-11-11", "2019-12-31", "2022-08-08")),
  p40000_i0 = as.Date(c(NA, "2022-05-05", "2023-08-08", NA, "2024-09-09", NA, NA, "2021-01-01")),
  p40000_i1 = as.Date(c("2023-01-01", NA, NA, "2024-04-04", NA, "2025-12-12", NA, NA)),

  # self-reported first occurrence (example; includes an abnormal placeholder)
  p20002_0_0 = as.Date(c("2015-01-01", "2016-02-02", "1900-01-01", NA, "2018-03-03", "2019-04-04", "2020-06-06", "1901-01-01")),

  # ICD10 codes (pipe-separated) and corresponding diagnosis dates
  p41270 = c(
    "D261|E113|M102",
    "E113|M101",
    NA,
    "D261|M102|M103",
    "E113|M101|D261",
    "M101|M100",
    "M101",
    "D261"
  ),
  p41280_a0 = as.Date(c("2012-05-01", "2021-07-20", NA, "2015-03-01", "2017-07-05", "2012-02-02", "2013-04-04", "2009-09-09")),
  p41280_a1 = as.Date(c("2012-03-10", "2014-09-01", NA, "2016-06-01", "2018-08-10", "2014-05-05", NA, NA)),
  p41280_a2 = as.Date(c("2015-08-15", NA, NA, "2015-09-10", NA, NA, NA, NA)),
  p41280_a3 = as.Date(c(NA, NA, NA, NA, NA, NA, NA, NA))
)

get_max_date(ukb_sim, date_columns = c("p41280_a*"))
#> [1] "None date in last 5 years"
```
