# Calculate Triglyceride-Glucose Index (TyG)

This function calculates the Triglyceride-Glucose Index (TyG) using
`triglyceride` and `glucose` values.

## Usage

``` r
leo_TyG(triglycerides, glucose, ...)
```

## Arguments

- triglycerides:

  A numeric vector of triglyceride values (in mmol/L).

- glucose:

  A numeric vector of glucose values (in mmol/L).

## Value

A numeric vector of TyG index values.

## Note

The formula for TyG is: \$\$TyG = \log \left( \frac{triglycerides \times
glucose}{2} \right)\$\$ Triglycerides and glucose are converted from
mmol/L to mg/dL using the factors:

- Triglycerides: \$\$mg/dL = mmol/L \times 88.5704\$\$

- Glucose: \$\$mg/dL = mmol/L \times 18.0168\$\$ – In UKB, triglycerides
  (Field ID: 30870) and glucose (Field ID: 30740) are measured in
  mmol/L.

## References

<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3662289/> (Cardiovasc
Diabetol, 2013)
