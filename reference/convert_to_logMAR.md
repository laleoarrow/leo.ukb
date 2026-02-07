# Convert Visual Acuity Measurements to logMAR

Handles numeric decimals and special terms: CF, HM, LP, NLP, and Chinese
synonyms.

## Usage

``` r
convert_to_logMAR(data, acuity_cols, term_map = NULL)
```

## Arguments

- data:

  A data frame containing visual acuity data.

- acuity_cols:

  Character vector of column names in `data` to be converted.

## Value

A data frame with the specified columns converted to numeric logMAR
values.

## Examples

``` r
va_df <- tibble::tibble(
  id = 1:5,
  bcva_od = c("1.0","0.8","CF","LP","NLP"),
  bcva_os = c("0.6","HM","0.3","1.2","0.1")
)
convert_to_logMAR(va_df, c("bcva_od","bcva_os"))
#> # A tibble: 5 × 3
#>      id bcva_od bcva_os
#>   <int>   <dbl>   <dbl>
#> 1     1  0       0.222 
#> 2     2  0.0969  2.3   
#> 3     3  2       0.523 
#> 4     4  2.6    -0.0792
#> 5     5  2.9     1     
```
