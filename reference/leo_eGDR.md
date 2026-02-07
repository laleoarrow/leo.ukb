# Calculate Estimated Glucose Disposal Rate (eGDR)

This function calculates the eGDR using waist circumference,
hypertension status, and HbA1c values. It works element-wise for vectors
of HbA1c, waist, and hypertension.

## Usage

``` r
leo_eGDR(hba1c, waist, hypertension, ...)
```

## Arguments

- hba1c:

  A numeric vector of HbA1c values (in mmol/mol).

- waist:

  A numeric vector of waist circumferences (in cm).

- hypertension:

  A numeric vector of hypertension status (1 for yes, 0 for no).

## Value

A numeric vector of eGDR values.

## Note

### eGDR: Estimated Glucose Disposal Rate:

    # eGDR is calculated using the following formula:
    eGDR = 21.158 + (-0.09 * waist circumference [cm]) + (−3.407 * hypertension [yes=1/no=0]) + (−0.551 * HbA1c (%))

In UKB, HbA1c is measured in **mmol/mol (IFCC)**, and needs to be
converted to **% (NGSP)**.

- **HbA1c (mmol/mol)** (Field ID: 30750)

  - Conversion formula: NGSP (%) = 0.09148 \* IFCC + 2.152
    [Ref](https://ngsp.org/ifcc.asp)

- **Waist circumference** (Field ID: 48)

- **hypertension** need to extract from icd10 or first occurrence
  records first.

## References

1.  <https://diabetesjournals.org/care/article/36/8/2280/32950/Use-of-the-Estimated-Glucose-Disposal-Rate-as-a>
    (Diabetes Care, 2013 July)

2.  <https://pmc.ncbi.nlm.nih.gov/articles/PMC11439291/#Sec2>
    (Cardiovasc Diabetol, 2024 Sep)
