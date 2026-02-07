# Robust dx runner

Executes dx CLI commands with automatic path detection and logging.
Supports macOS Homebrew installations and dry-run mode for testing.

## Usage

``` r
dx_run(
  args,
  intern = TRUE,
  ignore.stderr = FALSE,
  dry_run = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- args:

  Character vector of dx command arguments (e.g.,
  `c("describe", "file-123")`).

- intern:

  Logical. If TRUE, capture stdout as character vector. If FALSE, return
  exit code.

- ignore.stderr:

  Logical. If TRUE, suppress stderr output.

- dry_run:

  Logical. If TRUE, only logs the command without executing it. Returns
  mock success.

- verbose:

  Logical. If TRUE (default), logs the dx command before execution.

- ...:

  Additional arguments passed to
  [`system2()`](https://rdrr.io/r/base/system2.html).

## Value

If `intern=TRUE`, returns character vector of stdout. If `intern=FALSE`,
returns exit code. In dry_run mode, returns `character(0)` if
intern=TRUE, or `0L` if intern=FALSE.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Get current project info
  env_info <- dx_run(c("env"), intern = TRUE)

  # Extract project ID from environment
  # Output format: "Current workspace\tproject-XXX"
  workspace_line <- grep("Current workspace\\t", env_info, value = TRUE)
  project_id <- sub(".*\\t(project-[^\\s]+).*", "\\1", workspace_line)

  # Now use the extracted project ID
  result <- dx_run(c("describe", project_id), intern = TRUE)

  # Find datasets in current project
  datasets <- dx_run(c("find", "data", "--name", "*.dataset", "--brief"), intern = TRUE)

  # Upload file without capturing output (dry run)
  exit_code <- dx_run(c("upload", "file.txt", "--path", "/data/"), intern = FALSE, dry_run = TRUE)

  # Test command without execution (dry run)
  dx_run(c("run", "app-table-exporter", "-y"), dry_run = TRUE)

  # Ignore stderr for commands that may produce warnings
  dx_run(c("find", "data", "--brief"), ignore.stderr = TRUE)
} # }
```
