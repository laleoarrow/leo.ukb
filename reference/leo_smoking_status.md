# Recode smoking status (Never/Previous/Current) https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116

Recode smoking status (Never/Previous/Current)
https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116

## Usage

``` r
leo_smoking_status(smoking_status, ...)
```

## Arguments

- smoking_status:

  Numeric/character vector, usually from p20116.

## Value

Factor with levels "Never", "Previous", "Current"; NA for
negative/special codes or others.
