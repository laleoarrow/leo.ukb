# Recode employment status

Recode employment status into paid vs not-paid employment
https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6142

## Usage

``` r
leo_career(p6142_i0, ...)
```

## Arguments

- p6142_i0:

  Current employment status (UKB 6142, instance 0). Values: 1 = In paid
  employment or self-employed 2 = Retired 3 = Looking after home and/or
  family 4 = Unable to work because of sickness or disability 5 =
  Unemployed 6 = Doing unpaid or voluntary work 7 = Full or part-time
  student -7 = None of the above -3 = Prefer not to answer

- ...:

  Reserved for future use.

## Value

Factor with categories "Paid employment" and "Not paid employment"; NA
for special codes.
