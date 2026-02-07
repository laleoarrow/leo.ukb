# Step 1: Check Normality for Variables

This function checks the normality of each variable in the dataset and
returns a summary.

## Usage

``` r
leo.table1.step1(df, id_to_exclude = c("eid", "iri_time"), num_var = NULL)
```

## Arguments

- df:

  A data frame or tibble containing the variables to be tested.

- num_var:

  user-provided numeric variables to test for normality. Default is
  NULL.

- id:

  A character string specifying the column name for patient ID. Default
  is "eid".

## Value

A data frame with variables, p-values for normality test, and normality
status.

## Examples

``` r
# leo.table1.step1(t1d_cohort, id_to_exclude = c("eid", "iri_time"), num_var = c("age", "tdi", "score_diet", "TyG"))
```
