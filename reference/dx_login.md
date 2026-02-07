# Log in to DNAnexus

This is a wrapper for `dx login`. It supports both interactive login (if
no token provided) and token-based login.

## Usage

``` r
dx_login(token = NULL)
```

## Arguments

- token:

  Optional. A DNAnexus authentication token.

## Value

None (runs system command).
