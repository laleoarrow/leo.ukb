# Recode UKB qualifications into 3 categories

Recode UKB field 6138 (Qualifications; multiple selections separated by
"\|") into Degree / No degree / Unknown.

## Usage

``` r
leo_education(p6138_i0, ...)
```

## Arguments

- p6138_i0:

  UKB 6138 qualifications string
  (https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=6138).

- ...:

  Reserved for future use.

## Value

Character vector with values: "Degree", "No degree", or "Unknown".
