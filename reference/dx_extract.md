# Extract UKB data using DNAnexus Table Exporter

This function takes a file containing UKB field IDs and runs the Table
Exporter app on DNAnexus RAP to extract the data. The output will be
saved to `data_files/` on the RAP project (or the specified project).
Reference:
<https://github.com/UK-Biobank/UKB-RAP-Notebooks-Access/blob/main/RStudio/A110_Export_participant_data.Rmd>

## Usage

``` r
dx_extract(
  field_id = NULL,
  category_id = NULL,
  output_prefix = NULL,
  title_pattern_to_keep = NULL,
  title_pattern_to_exclude = NULL,
  expand = TRUE,
  project = NULL,
  dataset = NULL,
  entity = "participant",
  coding_option = "RAW",
  output_format = "CSV",
  header_style = "FIELD-NAME",
  instance_type = "mem1_ssd1_v2_x4",
  dry_run = FALSE
)
```

## Arguments

- field_id:

  Path to a file or a vector of Field IDs (e.g., `c("p31", "41270")`).

- category_id:

  Path to a file or a vector of numeric Category IDs (e.g., `1712`,
  `c(1712, 100)`). Category IDs must be integers and correspond to UKB
  data categories in the hierarchy schema.

- output_prefix:

  Output file prefix. If NULL, auto-generated based on
  fields/categories.

- title_pattern_to_keep:

  Optional regex to keep category-expanded fields by `title`
  (case-insensitive). Only applies when `category_id` is provided.

- title_pattern_to_exclude:

  Optional regex to drop category-expanded fields by `title`
  (case-insensitive). Applied after `title_pattern_to_keep`.

- expand:

  Logical. If TRUE (default), validates and expands field IDs (e.g.
  `21003`) into all instances/arrays present in the dataset (e.g.
  `p21003_i0`, `p21003_i1`) using the official dataset dictionary
  verification strategy.

- project:

  Optional target project ID (e.g., `"project-XXX"`). If `dataset`
  includes a project prefix, it must match `project` (otherwise error).
  If `dataset` is a name only and `project` is provided, the dataset is
  resolved within that project. To list accessible projects:
  `projects <- dx_run(c("find", "projects", "--brief"), intern = TRUE)`
  `projects <- dx_run(c("find", "projects"), intern = TRUE)` \# more
  details

- dataset:

  The dataset ID/name to extract from (e.g.,
  `"app12345_20240101.dataset"` or `"project-Gk2...:record-Fp3..."`).
  Defaults to the latest dataset in the project if NULL. To list
  datasets across all projects:
  `datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief", "--all-projects"), intern = TRUE)`

- entity:

  Entity type. Options:

  - `"participant"` (default): Main phenotype/covariate data.

  - `"hesin"`: Hospital episode statistics.

  - `"death"`: Death registry records.

  - `"image"`: Imaging metadata.

- coding_option:

  How to encode categorical values. Options:

  - `"RAW"` (default): Return raw numeric codes (e.g., 1, 2, 3).

  - `"REPLACE"`: Replace codes with human-readable labels (e.g., "Male",
    "Female").

- output_format:

  Output file format. Options:

  - `"CSV"` (default): Comma-separated values.

  - `"TSV"`: Tab-separated values.

- header_style:

  Column header naming style. Options:

  - `"FIELD-NAME"` (default): Use field names like `p21003_i0`.

  - `"FIELD-TITLE"`: Use descriptive titles like
    `"Age when attended assessment centre"`.

  - `"UKB-FORMAT"`: Use raw UKB field IDs like `21003-0.0`.

- instance_type:

  The DNAnexus instance type for the job. Default: `"mem1_ssd1_v2_x4"`.
  For large extracts, `"mem1_hdd1_v2_x8"` is recommended.

- dry_run:

  Logical. If TRUE, only skips the Table Exporter job submission. All
  other steps (dataset detection, dictionary validation, file upload)
  still run. Useful for testing command generation without launching
  jobs. Default is FALSE.

## Value

A list containing the job ID and expected output path (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
dx_login() # Make sure you have logged in to DNAnexus
dx_status() # Check the status of your login

# 1. Basic extraction (extracts covariates defined in fields.txt)
dx_extract("tmp/fields.txt", output_prefix = "ukb_data", dry_run = TRUE)

# 2. Extraction using vectors (Field & Category)
dx_extract(field_id = c("p53", "p191", "p40000", # basic demographics
                        "p41270", "p41271", "p41280", "p41281"), # ICD codes
           category_id = 1712, # First occurrences
           title_pattern_to_exclude = "^Source",
           output_prefix = "diagnosis20260205")

# 3. Just Category (numeric ID or vector of IDs)
dx_extract(category_id = 1712)  # Extract all fields in category 1712
dx_extract(category_id = c(1712, 100))  # Multiple categories

# 3b. First occurrences (FO) title filters
# Keep only titles starting with "Date"
dx_extract(category_id = 1712, title_pattern_to_keep = "^Date")
# Exclude titles starting with "Source"
dx_extract(category_id = 1712, title_pattern_to_exclude = "^Source")

# 4. Extract with human-readable labels and descriptive headers
dx_extract("tmp/fields.txt",
           output_prefix = "ukb_data_readable",
           coding_option = "REPLACE",
           header_style = "FIELD-TITLE")

# 5a. Explicitly specifying the dataset (Recommended for reproducibility)
datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE)
dx_extract("tmp/fields.txt", dataset = datasets[1])

# 5b. Extract in a specified project
projects_all <- dx_run(c("find", "projects"), intern = TRUE); projects_all
projects_all_brief <- dx_run(c("find", "projects", "--brief"), intern = TRUE); projects_all_brief
datasets_all <- dx_run(c("find", "data", "--name", "*.dataset", "--brief", "--all-projects"), intern = TRUE)
dx_extract(field_id = c("p53", "p191", "p40000", # basic demographics
                        "p41270", "p41271", "p41280", "p41281"), # ICD codes
           category_id = 1712, # First occurrences
           title_pattern_to_exclude = "^Source",
           output_prefix = "diagnosis20260205",
           project = projects_all_brief[2])

# 6. Dry run mode - Test without actually submitting the job
# Useful for debugging and verifying the command before execution
result <- dx_extract(
  field_id = c("p31", "p21003"),
  category_id = 1712,
  output_prefix = "test_run",
  dry_run = TRUE
)
# Shows all commands that would be executed
# Returns: list(job_id = "job-DRYRUN123456", output_path = "...", dry_run = TRUE)
} # }
```
