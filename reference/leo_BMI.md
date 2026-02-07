# Recode BMI into categories (Underweight/Normal weight/Overweight/Obese) https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001

Recode BMI into categories (Underweight/Normal weight/Overweight/Obese)
https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001

## Usage

``` r
leo_BMI(BMI, ...)
```

## Arguments

- BMI:

  Numeric/character vector of BMI (kg/m^2) from p21001.

## Value

Factor with levels "Underweight", "Normal weight", "Overweight",
"Obese"; NA if invalid or missing.
