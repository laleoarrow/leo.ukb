# Read UK Biobank HLA v2 header (ukb_hla_v2.txt) **\[stable\]**

Downloads and parses the HLA v2 header file used to decode UK Biobank
Field 22182 (HLA imputation v2; best-guess probabilities).

## Usage

``` r
ukb_hla_header()
```

## Value

A character vector of allele column names (e.g., "A_101", "B_2705",
...).

## Details

The returned character vector is the column order (length typically 362)
corresponding to the comma-separated values in field 22182.

The function tries
[`readLines()`](https://rdrr.io/r/base/readLines.html) first. If that
fails, it falls back to a system `wget` call (useful on clusters). If
both fail, it errors.

## Examples

``` r
header <- ukb_hla_header()
#> ℹ [16:10:45] Loading UKB HLA header; description source: https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=2182
head(header)
#> [1] "A_101" "A_102" "A_103" "A_201" "A_202" "A_203"
```
