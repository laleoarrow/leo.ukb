# Linear regression analysis **\[experimental\]**

This function performs linear regression analysis between PRS and
clinical outcomes with options for handling missing values through
complete-case analysis or imputation.

## Usage

``` r
leo_linear_regression(
  data,
  y_var,
  x_var,
  covariates = NULL,
  impute_na = FALSE,
  impute_method = "mean",
  ...
)
```

## Arguments

- data:

  A data frame containing the analysis data

- y_var:

  The name of the outcome variable (response)

- x_var:

  The name of the predictor variable (PRS)

- covariates:

  Optional vector of covariate names to adjust for

- impute_na:

  Logical indicating whether to impute missing values (default: FALSE)

- impute_method:

  Imputation method to use if impute_na = TRUE

- ...:

  Additional arguments passed to the imputation function

## Value

A list containing regression results including model, coefficients, and
fit statistics

## Examples

``` r
# Complete case analysis
result_cc <- leo_linear_regression(prs_zh8, "va1y", "PRS",
                                     covariates = c("age_onset", "gender"))
#> Error: object 'prs_zh8' not found

# With mean imputation
result_imputed <- leo_linear_regression(prs_zh8, "va1y", "PRS",
                                         covariates = c("age_onset", "gender"),
                                         impute_na = TRUE, impute_method = "mean")
#> Error: object 'prs_zh8' not found
```
