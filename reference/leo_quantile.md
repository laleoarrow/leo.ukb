# Bin a continuous variable by quantiles

A lightweight wrapper of [`cut()`](https://rdrr.io/r/base/cut.html)
using quantile-based breaks (equal-frequency).

## Usage

``` r
leo_quantile(
  x,
  bins = 4,
  labels = NULL,
  right = TRUE,
  include_lowest = TRUE,
  ...
)
```

## Arguments

- x:

  Numeric/character vector.

- bins:

  Integer; number of quantile bins (default 4).

- labels:

  Character vector; default uses "Q1..Qk".

- right:

  Logical; passed to [`cut()`](https://rdrr.io/r/base/cut.html).

- include_lowest:

  Logical; passed to [`cut()`](https://rdrr.io/r/base/cut.html).

- ...:

  Reserved for future use.

## Value

Character vector of quantile bins; NA preserved.

## Examples

``` r
x <- c(1:10, NA_integer_)
q4 <- leo_quantile(x, bins = 4); q4
#>  [1] "Q1" "Q1" "Q1" "Q2" "Q2" "Q3" "Q3" "Q4" "Q4" "Q4" NA  
table(q4, useNA = "ifany")
#> q4
#>   Q1   Q2   Q3   Q4 <NA> 
#>    3    2    2    3    1 

x_tie <- c(rep(1L, 6), rep(2L, 6))
q4_tie <- leo_quantile(x_tie, bins = 4); q4_tie
#> ! leo_quantile(): duplicated break points detected (ties). Collapsing duplicates.
#>  [1] "Q1" "Q1" "Q1" "Q1" "Q1" "Q1" "Q2" "Q2" "Q2" "Q2" "Q2" "Q2"
table(q4_tie, useNA = "ifany")
#> q4_tie
#> Q1 Q2 
#>  6  6 
```
