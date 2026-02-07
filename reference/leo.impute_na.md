# Basic imputation function for missing values

This function provides multiple methods for imputing missing values in a
dataset, including mean, median, random forest, KNN, and multiple
imputation.

## Usage

``` r
leo.impute_na(data, method = "mean", ...)
```

## Arguments

- data:

  A data frame containing the data with missing values

- method:

  Imputation method: "mean", "median", "rf" (random forest), "knn"
  (K-nearest neighbors), or "mice" (multiple imputation)

- ...:

  Additional arguments passed to the specific imputation function

## Value

A data frame with missing values imputed

## Examples

``` r
# Create sample data with missing values
set.seed(123)
sample_data <- data.frame(
  age = c(25, 30, NA, 40, 45),
  score = c(80, NA, 90, 85, NA),
  group = factor(c("A", "B", "A", NA, "B"))
)

# Mean imputation
leo.impute_na(sample_data, method = "mean")
#> ℹ [15:28:01] Missing values per column (count/percentage):
#> ℹ [15:28:01]   age: 1/5 (20%)
#> ℹ [15:28:01]   score: 2/5 (40%)
#> ℹ [15:28:01]   group: 1/5 (20%)
#> ! [15:28:01] Remaining missing values after imputation:
#> ! [15:28:01]   group: 1
#>   age score group
#> 1  25    80     A
#> 2  30    85     B
#> 3  35    90     A
#> 4  40    85  <NA>
#> 5  45    85     B

# Random forest imputation
leo.impute_na(sample_data, method = "rf")
#> ℹ [15:28:01] Missing values per column (count/percentage):
#> ℹ [15:28:01]   age: 1/5 (20%)
#> ℹ [15:28:01]   score: 2/5 (40%)
#> ℹ [15:28:01]   group: 1/5 (20%)
#> Missing value imputation by random forests
#> 
#> Variables to impute:     age, group, score
#> Variables used to impute:    age, score, group
#> 
#> iter 1 
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> ℹ [15:28:01] Imputation completed using rf method
#>   age score group
#> 1  25 80.00     A
#> 2  30 84.91     B
#> 3  40 90.00     A
#> 4  40 85.00     A
#> 5  45 84.91     B
```
