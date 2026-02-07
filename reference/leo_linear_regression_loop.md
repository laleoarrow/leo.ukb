# Batch linear regression analysis for multiple predictors

This function performs linear regression analysis for multiple predictor
variables against a single outcome, with consistent formatting of
results.

## Usage

``` r
leo_linear_regression_loop(data, y_var, x_vars, covariates = NULL, ...)
```

## Arguments

- data:

  Data frame containing all variables

- y_var:

  Outcome variable name

- x_vars:

  Vector of predictor variable names to test

- covariates:

  Covariates to adjust for in all models

- ...:

  Additional arguments passed to leo_linear_regression

## Value

A data frame with consolidated results from all regression models

## Examples

``` r
# Analyze multiple PRS methods against the same outcome
x_vars <- c("PRS_plink", "PRS_catboost", "PRS_lasso")
results <- leo_linear_regression_loop(prs_zh8, "va1y", x_vars,
                                     covariates = c("age_onset", "gender"))
#> ℹ [15:28:04] leo_linear_regression_loop function placeholder - to be implemented
#> ℹ [15:28:04] Would analyze 3 predictors against outcome va1y
```
