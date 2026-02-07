# Recode ethnicity into finer categories (top-level 1–6) https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000

Recode ethnicity into finer categories (top-level 1–6)
https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000

## Usage

``` r
leo_ethnicity_finer(ethnicity, ...)
```

## Arguments

- ethnicity:

  Numeric/character vector, usually from p21000.

## Value

Factor with levels: "White", "Mixed", "Asian or Asian British", "Black
or Black British", "Chinese", "Other ethnic group"; NA for
-1/-3/missing/other unexpected codes.
