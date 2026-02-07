# Process multiple ICD definitions in bulk

Apply diagnosis_process_icd() to a list or table of ICD code sets and
return one status/time pair per diagnosis.

## Usage

``` r
diagnosis_process_icd_multi(
  data,
  icd_map,
  icd = 10,
  censored = "2025-01-01",
  icd_rank = 1,
  baseline_date = "p53_i0"
)
```

## Arguments

- data:

  A data.frame containing UKB baseline, censoring, ICD code and date
  columns.

- icd_map:

  A named list (diagnosis -\> vector of ICD codes) or a data.frame. If
  data.frame: first column is diagnosis name; remaining columns or a
  single "codes" column contain ICD codes (comma/space separated or one
  per cell).

- icd:

  ICD version, 10 or 9.

- censored:

  Censoring date ("YYYY-MM-DD").

- icd_rank:

  Rank of diagnosis occurrence to consider (1 = first occurrence).

- baseline_date:

  Baseline date column name (default "p53_i0").

## Value

A data.frame with eid and diagnosis_status/diagnosis_time columns.
