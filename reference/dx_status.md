# Check DNAnexus Job Status & Dashboard

Checks the status of a specific job or displays a dashboard of recent
jobs.

## Usage

``` r
dx_status(all_projects = TRUE, job_id = NULL, limit = 5)
```

## Arguments

- all_projects:

  Logical. If TRUE, shows jobs from all projects in the account
  (account-wide). Default is FALSE (only current project).

- job_id:

  Optional. A string containing the job ID (e.g., "job-Gk2...") or a
  list returned by
  [`dx_extract`](https://laleoarrow.github.io/leo.ukb/reference/dx_extract.md).
  If provided, `all_projects` is ignored.

- limit:

  Numeric. Number of recent jobs to show in the dashboard (default: 3).

## Value

The state of the job as a string (if `job_id` provided) or a data frame
of recent jobs (if `job_id` is NULL), invisibly.

## Details

When called without a `job_id`, it displays a "htop"-style dashboard of
the last several jobs.

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. View dashboard of recent jobs in current project
dx_status()

# 2. View dashboard for ONLY the current project
dx_status(all_projects = FALSE)

# 3. Check status of a specific job
dx_status(job_id = "job-Gk2PzX0Jj1X4Y5Z6")
} # }
```
