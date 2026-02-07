# Clean HLA carrier output to a wide format

Convert the long table returned by
[`ukb_hla_specify()`](https://laleoarrow.github.io/leo.ukb/reference/ukb_hla_specify.md)
into a wide table suitable for downstream analysis:
`eid, <allele>, <allele>_carrier, ...`

## Usage

``` r
ukb_hla_clean(res)
```

## Arguments

- res:

  A data.frame with columns `eid`, `query`, `copies`, `carrier`.

## Value

A wide-format data.frame with one row per `eid`.

## Examples

``` r
hla <- fread("/Users/leoarrow/Project/UKB/data/HLA.csv") %>%
  ukb_hla_typing(header) %>%
  ukb_hla_specify("B27") %>%
  ukb_hla_clean()
#> Error in fread("/Users/leoarrow/Project/UKB/data/HLA.csv") %>% ukb_hla_typing(header) %>%     ukb_hla_specify("B27") %>% ukb_hla_clean(): could not find function "%>%"
```
