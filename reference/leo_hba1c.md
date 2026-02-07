# Recode HbA1c (IFCC, mmol/mol) https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30750

UK Biobank HbA1c is stored in mmol/mol (IFCC). This function performs:

- numeric coercion,

- recoding UKB special codes (-3/-1) to NA,

- dropping non-positive values to NA (conservative guard).

## Usage

``` r
leo_hba1c(hba1c, ...)
```

## Arguments

- hba1c:

  Numeric/character vector of HbA1c values (mmol/mol, IFCC), usually
  from p30750.

- ...:

  Reserved for future use.

## Value

Numeric vector of HbA1c in mmol/mol (IFCC), with special codes as NA.
