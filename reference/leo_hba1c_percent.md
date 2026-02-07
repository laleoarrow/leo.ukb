# Convert HbA1c from IFCC (mmol/mol) to NGSP (%) https://ngsp.org/ifcc.asp

Conversion formula: \$\$NGSP(\\) = 0.09148 \times IFCC(mmol/mol) +
2.152\$\$

## Usage

``` r
leo_hba1c_percent(hba1c, ...)
```

## Arguments

- hba1c:

  Numeric/character vector of HbA1c values (mmol/mol, IFCC), usually
  from p30750.

- ...:

  Reserved for future use.

## Value

Numeric vector of HbA1c in % (NGSP), NA preserved.
