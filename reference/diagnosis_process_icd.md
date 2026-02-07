# Process UKB ICD data for survival analysis

Generate icd_code_status and icd_code_time from UKB ICD records.

## Usage

``` r
diagnosis_process_icd(
  data,
  icd_code,
  icd = 10,
  censored = "2025-01-01",
  icd_rank = 1,
  baseline_date = "p53_i0",
  filter_NA = FALSE
)
```

## Arguments

- data:

  A data.frame containing UKB baseline, censoring, ICD code and date
  columns.

- icd_code:

  A target ICD code prefix or a character vector of prefixes (e.g.,
  "M10" or c("E10", "E11")).

- icd:

  ICD version, 10 or 9.

- censored:

  Censoring date ("YYYY-MM-DD").

- icd_rank:

  Rank of diagnosis occurrence to consider (1 = first occurrence, 2 =
  second occurrence, etc.).

- baseline_date:

  Baseline date column name (default "p53_i0").

- filter_NA:

  Logical, whether to filter out NA row. Default FALSE (Suitable for
  comorbidity check).

## Value

A data.frame with eid, icd_code_status, icd_code_time.

## Examples

``` r
df_example <- data.frame(
  eid = 1:5,
  p53_i0 = as.Date(c("2006-01-01", "2007-06-01", "2005-03-15", "2004-10-30", "2008-05-10")),  # baseline
  p191 = as.Date(c("2020-01-01", "2019-12-31", "2021-06-30", "2022-07-01", "2023-09-10")),    # last follow-up
  p40000_i0 = as.Date(c(NA, NA, NA, "2021-06-01", NA)),                                       # death date (field 0)
  p40000_i1 = as.Date(c(NA, NA, NA, NA, NA)),                                                 # death date (field 1)
  p41270 = c("D261|E113|M102", "E113|M101", NA, "D261|M102|M103", "E113|M101|D261"),          # ICD10 codes
  p41280_a0 = as.Date(c("2010-05-01", "2011-07-20", NA, "2015-03-01", "2017-07-05")),         # 1st diagnosis date
  p41280_a1 = as.Date(c("2012-03-10", "2014-09-01", NA, "2016-06-01", "2018-08-10")),         # 2nd diagnosis date
  p41280_a2 = as.Date(c("2015-08-15", NA, NA, "2015-09-10", NA)),                             # 3rd diagnosis date
  p41280_a3 = as.Date(c(NA, NA, NA, NA, NA))                                                  # Additional diagnosis field for testing
)

# Prefix matching: D261 matches only D261 (icd_rank = 1)
diagnosis_process_icd(df_example, icd_code = "D261", icd = 10, icd_rank = 1)
#> ℹ [15:53:42] diagnosis_process_icd: processing ICD10 [D261] using [2025-01-01] as censored date.
#> ℹ [15:53:42] diagnosis_process_icd: converting date columns to <Date>.
#> ℹ [15:53:42] diagnosis_process_icd: combining death dates.
#> ℹ [15:53:42] diagnosis_process_icd: matching ICD10 pattern [^(D261)].
#> ℹ [15:53:42] diagnosis_process_icd: calculating survival time and status.
#> ✔ [15:53:42] diagnosis_process_icd: created [D261_status], [D261_time] - ALL DONE.
#>   eid D261_status D261_time
#> 1   1           1  4.331507
#> 2   2           0 12.591781
#> 3   3           0 16.304110
#> 4   4           1 10.339726
#> 5   5           0 15.345205

# Prefix matching: M10 matches M101 / M102 (icd_rank = 1)
diagnosis_process_icd(df_example, icd_code = "M10", icd = 10, icd_rank = 1)
#> ℹ [15:53:42] diagnosis_process_icd: processing ICD10 [M10] using [2025-01-01] as censored date.
#> ℹ [15:53:42] diagnosis_process_icd: converting date columns to <Date>.
#> ℹ [15:53:42] diagnosis_process_icd: combining death dates.
#> ℹ [15:53:42] diagnosis_process_icd: matching ICD10 pattern [^(M10)].
#> ℹ [15:53:42] diagnosis_process_icd: calculating survival time and status.
#> ✔ [15:53:42] diagnosis_process_icd: created [M10_status], [M10_time] - ALL DONE.
#>   eid M10_status  M10_time
#> 1   1          1  9.624658
#> 2   2          1  7.257534
#> 3   3          0 16.304110
#> 4   4          1 10.868493
#> 5   5          1 10.257534

# Prefix matching: E113 matches E113 only (icd_rank = 1)
diagnosis_process_icd(df_example, icd_code = "E113", icd = 10, icd_rank = 1)
#> ℹ [15:53:42] diagnosis_process_icd: processing ICD10 [E113] using [2025-01-01] as censored date.
#> ℹ [15:53:42] diagnosis_process_icd: converting date columns to <Date>.
#> ℹ [15:53:42] diagnosis_process_icd: combining death dates.
#> ℹ [15:53:42] diagnosis_process_icd: matching ICD10 pattern [^(E113)].
#> ℹ [15:53:42] diagnosis_process_icd: calculating survival time and status.
#> ✔ [15:53:42] diagnosis_process_icd: created [E113_status], [E113_time] - ALL DONE.
#>   eid E113_status E113_time
#> 1   1           1  6.191781
#> 2   2           1  4.136986
#> 3   3           0 16.304110
#> 4   4           0 16.597260
#> 5   5           1  9.158904

# Prefix matching: M10 matches M101 / M102 (icd_rank = 2)
diagnosis_process_icd(df_example, icd_code = "M10", icd = 10, icd_rank = 2)
#> ℹ [15:53:42] diagnosis_process_icd: processing ICD10 [M10] using [2025-01-01] as censored date.
#> ℹ [15:53:42] diagnosis_process_icd: converting date columns to <Date>.
#> ℹ [15:53:42] diagnosis_process_icd: combining death dates.
#> ℹ [15:53:42] diagnosis_process_icd: matching ICD10 pattern [^(M10)].
#> ℹ [15:53:42] diagnosis_process_icd: calculating survival time and status.
#> ✔ [15:53:42] diagnosis_process_icd: created [M10_status], [M10_time] - ALL DONE.
#>   eid M10_status M10_time
#> 1   1          0 14.00822
#> 2   2          0 12.59178
#> 3   3          0 16.30411
#> 4   4          1 11.59452
#> 5   5          0 15.34521

# Multiple ICD codes: c("M10", "E113") - takes global earliest date among all matches
diagnosis_process_icd(df_example, icd_code = c("M10", "E113"), icd = 10, icd_rank = 1)
#> ℹ [15:53:42] diagnosis_process_icd: processing ICD10 [M10 and E113] using [2025-01-01] as censored date.
#> ℹ [15:53:42] diagnosis_process_icd: converting date columns to <Date>.
#> ℹ [15:53:42] diagnosis_process_icd: combining death dates.
#> ℹ [15:53:42] diagnosis_process_icd: matching ICD10 pattern [^(M10|E113)].
#> ℹ [15:53:42] diagnosis_process_icd: calculating survival time and status.
#> ✔ [15:53:42] diagnosis_process_icd: created [M10_E113_status], [M10_E113_time] - ALL DONE.
#>   eid M10_E113_status M10_E113_time
#> 1   1               1      6.191781
#> 2   2               1      4.136986
#> 3   3               0     16.304110
#> 4   4               1     10.868493
#> 5   5               1      9.158904
# Output columns: M10_E113_status, M10_E113_time
```
