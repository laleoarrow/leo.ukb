# Combine Measurements from Two Eyes

Combines two numeric vectors representing measurements from two eyes
into a single vector using a specified method (mean, max, or min).
Handles NA values appropriately.

## Usage

``` r
combine_two_eye(vector1, vector2, method = c("mean", "max", "min"))
```

## Arguments

- vector1:

  Numeric vector for the first eye.

- vector2:

  Numeric vector for the second eye.

- method:

  Method to combine the two vectors: "mean", "max", or "min". Default is
  "mean".

## Value

A numeric vector with combined values.
