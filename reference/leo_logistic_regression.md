# General looped logistic regression: y (binary) ~ x (+ covariates)

Loop over multiple exposures (x) to fit logistic regressions against a
single binary outcome (y). Input data frames must have the first column
as individual ID. Covariates are optional.

## Usage

``` r
leo_logistic_regression(
  x,
  y,
  cov = NULL,
  x_col = NULL,
  y_col = NULL,
  id_col = "id",
  y_level = c("0", "1"),
  y_label = c("0", "1"),
  scale_x = F,
  na_action = "na.omit",
  verbose = TRUE
)
```

## Arguments

- x:

  Data frame (first col = ID) of exposures; remaining columns are
  candidate predictors (numeric or factor).

- y:

  Data frame (first col = ID) of binary outcome; one column is used (see
  y_col). If numeric, it can be 0/1 or any two distinct numeric codes;
  if factor/character, must have exactly two levels.

- cov:

  Optional data frame (first col = ID) of covariates; remaining columns
  are covariates.

- x_col:

  Character vector of exposure column names in x; default: all columns
  from the 2nd.

- y_col:

  Outcome column name in y; default: the 2nd column of y.

- id_col:

  Shared ID column name in x/y/cov; default "id".

- y_level:

  Two labels found in y that define reference/contrast order (first =
  reference); default c("0","1").

- y_label:

  Descriptive names to print in logs parallel to y_level; default
  c("0","1").

- scale_x:

  Logical; if TRUE, scale() numeric exposures; default FALSE.

- na_action:

  One of "na.omit" or "na.exclude"; default "na.omit".

- verbose:

  Logical; if FALSE, suppress logs; default TRUE.

## Value

A tibble; each exposure may produce multiple rows when categorical:

- `x_name`: exposure variable name

- `x_type`: "numeric" or "factor"

- `x_level`: "(continuous)" for numeric; for factor: level
  name（参考水平标注为 "(ref)"）

- `ref_level`: reference level（factor 才有；numeric 为 NA）

- `y_name`: outcome variable name

- `n`: number of non-missing samples used

- `r2_mcfadden`: McFadden's pseudo-R² = 1 - (model deviance / null
  deviance)

- `beta`, `se`, `z`, `p`: coefficient on logit scale (per 1SD if
  `scale_x=TRUE`)

- `ci_l`, `ci_u`: 95\\

- `OR`, `OR_CI_l`, `OR_CI_u`: odds ratio and 95\\

- `p_adj`: BH-adjusted p-value across all rows in this run

## Examples

``` r
# --- Example 1: numeric outcome 0/1, numeric exposure (scaled), no covariates ---
set.seed(10)
x <- tibble::tibble(id = 1:300, PRS = rnorm(300))
linpred <- 0.8 * scale(x$PRS)[,1]
y <- tibble::tibble(id = 1:300, outcome = rbinom(300, 1, prob = stats::plogis(linpred)))
res1 <- leo_logistic_regression(x = x, y = y, x_col = "PRS", y_col = "outcome",
                                id_col = "id", scale_x = T,
                                y_level = c("0","1"), y_label = c("control","case"))
#> ℹ [15:53:46] Performing logistic regression with settings:
#> y = outcome (0=control, 1=case)
#> x = PRS
#> cov = (none)
#> ✔ [15:53:46] Completed. Processed 1 exposure(s). N range: 300–300.
head(res1)
#> # A tibble: 1 × 17
#>   x_name x_type  x_level    ref_level y_name     n r2_mcfadden  beta    se     z
#>   <chr>  <chr>   <chr>      <chr>     <chr>  <int>       <dbl> <dbl> <dbl> <dbl>
#> 1 PRS    numeric (continuo… NA        outco…   300      0.0980 0.803 0.137  5.87
#> # ℹ 7 more variables: p <dbl>, ci_l <dbl>, ci_u <dbl>, OR <dbl>, OR_CI_l <dbl>,
#> #   OR_CI_u <dbl>, p_adj <dbl>

# --- Example 2: factor exposure with explicit reference level, plus covariates ---
set.seed(11)
x2 <- tibble::tibble(id = 1:400,
                     group = factor(sample(c("A","B","C"), 400, T), levels = c("A","B","C")))  # A is ref
cov2 <- tibble::tibble(id = 1:400, age = rnorm(400, 50, 10))
y2  <- tibble::tibble(id = 1:400, outcome = rbinom(400, 1, 0.4))
res2 <- leo_logistic_regression(x = x2, y = y2, cov = cov2,
                                x_col = "group", y_col = "outcome",
                                id_col = "id", scale_x = F,
                                y_level = c("0","1"), y_label = c("control","case"))
#> ℹ [15:53:46] Performing logistic regression with settings:
#> y = outcome (0=control, 1=case)
#> x = group
#> cov = age
#> ! [15:53:46] Exposure 'group' is categorical; please ensure factor() with first level as reference. Levels: A, B, C (ref = A).
#> ℹ [15:53:46] For categorical exposure 'group', interpretation is vs reference 'A'.
#> ✔ [15:53:46] Completed. Processed 1 exposure(s). N range: 400–400.
res2
#> # A tibble: 3 × 17
#>   x_name x_type x_level ref_level y_name     n r2_mcfadden    beta     se      z
#>   <chr>  <chr>  <chr>   <chr>     <chr>  <int>       <dbl>   <dbl>  <dbl>  <dbl>
#> 1 group  factor (ref)   A         outco…   400     0.00404 NA      NA     NA    
#> 2 group  factor B       A         outco…   400     0.00404 -0.269   0.250 -1.07 
#> 3 group  factor C       A         outco…   400     0.00404 -0.0341  0.248 -0.138
#> # ℹ 7 more variables: p <dbl>, ci_l <dbl>, ci_u <dbl>, OR <dbl>, OR_CI_l <dbl>,
#> #   OR_CI_u <dbl>, p_adj <dbl>
```
