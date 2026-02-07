# Calculate Physical Activity Indicator **\[stable\]**

Here, using IPAQ minutes/week as input, classify regular activity
according to AHA (150/75/equivalent combination) rules:

- Moderate: \\\ge 150\\ min/week; or

- Vigorous: \\\ge 75\\ min/week; or

- Equivalent: \\moderate + 2 \times vigorous \ge 150.\\

## Usage

``` r
leo_physical_activity(p884_i0, p894_i0, p904_i0, p914_i0, ...)
```

## Arguments

- p884_i0:

  Moderate activity days/week (10+ min). UKB 884; -1/-3 -\> NA.

- p894_i0:

  Moderate activity minutes/day (typical day). UKB 894; -1/-3 -\> NA.

- p904_i0:

  Vigorous activity days/week (10+ min). UKB 904; -1/-3 -\> NA.

- p914_i0:

  Vigorous activity minutes/day (typical day). UKB 914; -1/-3 -\> NA.

- ...:

  Reserved for future use.

## Value

Factor with levels c("Non-regular", "Regular"); NA if insufficient data.

## Note

*Adults should pursue at least 150 minutes per week of
moderate-intensity physical activity, or 75 minutes per week of
vigorous-intensity aerobic physical activity, or an equivalent
combination of moderate- and vigorous-intensity aerobic activities.*
(Circulation, 2010,
[source](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.109.192703))

- **≥150 min/wk moderate intensity**

- **≥75 min/wk vigorous intensity**

- **Combination**
