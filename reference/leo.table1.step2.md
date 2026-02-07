# Step 2: Generate Table 1 for Cohort Study

This function generates a descriptive table for clinical cohort studies,
with options to stratify by group and perform group comparisons.

## Usage

``` r
leo.table1.step2(
  df,
  var_all,
  var_cat,
  var_non,
  strata = NULL,
  var_exact = NULL,
  compare_test = F,
  includeNA = F,
  showAllLevels = F,
  verbose = T
)
```

## Arguments

- df:

  A data frame or tibble containing the clinical variables.

- var_all:

  Vector. All var to be included in table 1

- var_cat:

  Vector. Categorical variables.

- var_non:

  Vector. Non-normal continuous variables. (Process in `print`)

- strata:

  Character, the column name to stratify by. Use "none" to disable
  stratification.

- var_exact:

  Vector. var_cat that needs fisher exact test. (Process in `print`)

- compare_test:

  Logical, if TRUE, performs group comparison (only works when `strata`
  is not "none").

- includeNA:

  Logical, if TRUE, includes NAs as a level in categorical variables.

- showAllLevels:

  Logical, if TRUE, all levels of factor variables are shown.

- verbose:

  Logical, if TRUE, prints the table.

## Value

Printed table1
