# Generate New Clinical Indicators

This function applies specific functions to generate new clinical
indicators. It is designed to work with vectors (e.g., columns of a
dataframe).

## Usage

``` r
leo_new_clinical_indicators(
  df,
  types = c("eGDR", "TyG"),
  type_id = TRUE,
  remove_assist = TRUE,
  remove_other = TRUE,
  ...
)
```

## Arguments

- df:

  A dataframe containing the required columns for calculation.

- types:

  A vector of indicator types, options are "eGDR", "TyG", "sleep_score",
  "LGI"...

- type_id:

  Logical; if `FALSE`, the function expects semantic column names (e.g.
  `"age"`, `"BMI"`). If `TRUE` (default), UKB field-id style names (e.g.
  `"p21003_i0"`) are used.

- remove_assist:

  Logical, whether to remove the assisting columns (default `TRUE`).

- remove_other:

  Logical, whether to keep only `eid` and the newly generated indicator
  columns in the returned dataframe (default `TRUE`). If `FALSE`, all
  original columns in `df` are kept.

- ...:

  Additional parameters passed to the specific indicator calculation
  functions.

## Value

A dataframe containing the new clinical indicators.

## Details

Supported indicators (use values in `types`):

- Basics: age, BMI, gender, ethnicity, ethnicity_finer.

- Socioeconomic: tdi, education, household_income, household_income_2,
  career.

- Inflammation: LGI.

- Insulin resistance: hba1c, hba1c_percent, eGDR, TyG.

- Lifestyle: smoking_status, drinking_status, diet_us,
  physical_activity, sleep_score.

## See also

[`leo_eGDR`](https://laleoarrow.github.io/leo.ukb/reference/leo_eGDR.md),
[`leo_TyG`](https://laleoarrow.github.io/leo.ukb/reference/leo_TyG.md)
for Insulin Resistance indicators.

## Examples

``` r
df <- data.frame(
  hba1c = c(48, 55, 60),
  waist = c(80, 85, 90),
  hypertension = c(1, 0, 1),
  triglycerides = c(1.5, 1.8, 2.0),
  glucose = c(5.0, 5.5, 6.0)
)
leo.ukb::leo_new_clinical_indicators(df, types = c("eGDR", "TyG"))
#> ── Generating new clinical indicator: eGDR ─────────────────────────────────────
#> ℹ Checking if all required column is valid for eGDR calculation.
#> ✖ Missing: p30750_i0, p48_i0. Skip eGDR
#> ── Generating new clinical indicator: TyG ──────────────────────────────────────
#> ℹ Checking if all required column is valid for TyG calculation.
#> ✖ Missing: p30870_i0, p30740_i0. Skip TyG
#> Error in dplyr::select(., eid, all_of(generated_types)): Can't select columns that don't exist.
#> ✖ Column `eid` doesn't exist.
```
