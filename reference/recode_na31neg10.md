# Helper function to recode UKB special values (-3/-1 to NA and -10 to 0)

Helper function to recode UKB special values (-3/-1 to NA and -10 to 0)

## Usage

``` r
recode_na31neg10(x)
```

## Details

For UK Biobank fields where -10 encodes "Less than once" and should be
treated as 0; -1 = "Do not know", -3 = "Prefer not to answer".
