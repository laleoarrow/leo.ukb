# Calculate Low-Grade Inflammation Score (LGI / INFLA-score)

This function calculates the low-grade inflammation score (INFLA-score)
using four biomarkers:

- high-sensitivity C-reactive protein (CRP),

- white blood cell (WBC) count,

- platelet (PLT) count,

- neutrophil-to-lymphocyte ratio (NLR = neutrophil / lymphocyte).

## Usage

``` r
leo_LGI(crp, wbc, plt, neutrophil, lymphocyte, ...)
```

## Arguments

- crp:

  Numeric vector of hs-CRP values (mg/L).

- wbc:

  Numeric vector of white blood cell counts (×10^9/L).

- plt:

  Numeric vector of platelet counts (×10^9/L).

- neutrophil:

  Numeric vector of neutrophil counts (×10^9/L).

- lymphocyte:

  Numeric vector of lymphocyte counts (×10^9/L).

- ...:

  Reserved for future use.

- x:

  Numeric vector of biomarker values.

## Value

Numeric vector of LGI / INFLA scores (-16 to +16).

Integer vector in -4, 4; NA if input is NA.

## Details

For each biomarker, within-sample deciles (1–10) are mapped to scores:

- deciles 1–4 -\> -4, -3, -2, -1

- deciles 5–6 -\> 0, 0

- deciles 7–10 -\> 1, 2, 3, 4

The LGI / INFLA-score is the sum of the four component scores and ranges
from -16 to +16. Higher values indicate higher chronic low-grade
inflammation.

## Note

Suggested UK Biobank fields (instance 0):

- CRP: Field 30710 (e.g. p30710_i0)

- WBC: Field 30000 (e.g. p30000_i0)

- Platelets: Field 30080 (e.g. p30080_i0)

- Neutrophil: Field 30140 (e.g. p30140_i0)

- Lymphocyte: Field 30120 (e.g. p30120_i0)

The fourth component uses NLR instead of GrL, taking advantage of
differential white cell counts in UK Biobank. Internal helper: score
biomarker by deciles for LGI
