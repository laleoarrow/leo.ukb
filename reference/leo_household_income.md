# Recode household income (UKB 738) into readable categories

Recode UKB field 738 "Average total household income before tax" into 5
brackets plus "Unknown" (-1/-3/NA/blank).

## Usage

``` r
leo_household_income(p738_i0, ...)
```

## Arguments

- p738_i0:

  Household income before tax (UKB 738, instance 0). Values: 1–5; -1/-3
  -\> Unknown.

- ...:

  Reserved for future use.

## Value

Character vector: "0-18,000", "18,000-30,999", "31,000-51,999",
"52,000-100,000", "Above 100,000", or "Unknown".
