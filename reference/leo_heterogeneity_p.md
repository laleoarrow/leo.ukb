# Calculate P-value for Heterogeneity from Subgroup Summary Statistics **\[stable\]**

This function performs a test for heterogeneity (interaction) based on
summary statistics (Hazard Ratios and P-values) from a subgroup analysis
of a Cox model. It uses Cochran's Q test.

## Usage

``` r
leo_heterogeneity_p(hrs, p_values, subgroup_names = NULL)
```

## Arguments

- hrs:

  A numeric vector of Hazard Ratios (HR) for each subgroup.

- p_values:

  A numeric vector of the corresponding P-values for each HR.

- subgroup_names:

  An optional character vector of subgroup names for display.

## Value

A list containing the Q statistic, degrees of freedom (df), the P-value
for heterogeneity, and the I-squared statistic. It also prints a summary
to the console.

## Examples

``` r
# Example: Subgroup analysis by smoking status (Smoker vs. Non-smoker)
hrs_example <- c(3.576652, 2.657723366)
p_values_example <- c(0.012254053, 0.002129672)
names_example <- c("Age>65", "Age<=65")
leo_heterogeneity_p(
  hrs = hrs_example,
  p_values = p_values_example,
  subgroup_names = names_example
)
#> --- Test for Heterogeneity in Subgroup Analysis ---
#> 
#> Input Data:
#>   Subgroup       HR     P_value log_HR SE_log_HR
#> 1   Age>65 3.576652 0.012254053 1.2744    0.5088
#> 2  Age<=65 2.657723 0.002129672 0.9775    0.3182
#> 
#> Results of Cochran's Q Test:
#> Q-statistic: 0.2448 
#> Degrees of Freedom (df): 1 
#> P-value for Heterogeneity: 0.621 
#> 
#> Heterogeneity Metrics:
#> I^2 statistic: 0% 
#> Interpretation of I^2: ~0-40% (might not be important); 30-60% (may represent moderate heterogeneity); 50-90% (may represent substantial heterogeneity); 75-100% (considerable heterogeneity).
#> 
#> ---------------------------------------------------
```
