# Recode drinking status (Never/Previous/Current) https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20117

Recode drinking status (Never/Previous/Current)
https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20117

## Usage

``` r
leo_drinking_status(drinking_status, ...)
```

## Arguments

- drinking_status:

  Numeric/character vector, usually from p20117.

## Value

Factor with levels "Never", "Previous", "Current"; NA for
negative/special codes or others.
