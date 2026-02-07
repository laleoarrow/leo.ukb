# Extract and format linear regression results for publication **\[experimental\]**

This function extracts key results from linear regression models and
formats them in a publication-ready table with standardized effect sizes
and confidence intervals.

## Usage

``` r
leo_regression_result(
  regression_result,
  round_digits = 3,
  include_ci = TRUE,
  include_stats = TRUE
)
```

## Arguments

- regression_result:

  List object returned by leo_linear_regression function

- round_digits:

  Number of digits to round numeric values (default: 3)

- include_ci:

  Whether to include confidence intervals (default: TRUE)

- include_stats:

  Whether to include model statistics (R², sample size) (default: TRUE)

## Value

A data frame with formatted regression results suitable for publication

## Examples

``` r
# Run regression
result <- leo_linear_regression(prs_zh8, "va1y", "PRS", covariates = c("age_onset", "gender"))
#> Error: object 'prs_zh8' not found

# Format results for publication
pub_table <- leo_regression_result(result)
#> Error: object 'result' not found
print(pub_table)
#> Error: object 'pub_table' not found
```
