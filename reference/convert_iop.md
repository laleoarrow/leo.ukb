# Convert IOP Measurements

This function converts intraocular pressure (IOP) measurements by
replacing "无" with NA and converting the rest to numeric.

## Usage

``` r
convert_iop(data, iop_cols)
```

## Arguments

- data:

  A data frame containing IOP data.

- iop_cols:

  A character vector of IOP column names in `data` to be converted.

## Value

A data frame with the specified IOP columns converted to numeric, with
"无" as NA.

## Examples

``` r
# Simulated IOP data，含单位与缺失写法
iop_df <- tibble::tibble(
  id = 1:6,
  iop_od = c("15","无","18mmHg",">30","-","N/A"),
  iop_os = c("16","未测","14.5","12mmHg","NA","20")
)
convert_iop(iop_df, c("iop_od","iop_os"))
#> # A tibble: 6 × 3
#>      id iop_od iop_os
#>   <int>  <dbl>  <dbl>
#> 1     1     15   16  
#> 2     2     NA   NA  
#> 3     3     18   14.5
#> 4     4     30   12  
#> 5     5     NA   NA  
#> 6     6     NA   20  
```
