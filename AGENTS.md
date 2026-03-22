# Repo Rules

## Tutorials

- Every public tutorial/article must have two separate versions: one English and one Chinese.
- English and Chinese tutorials should live in separate source files and separate generated pages.
- When a tutorial is added or updated, both language versions should be updated together.
- The pkgdown navbar should expose the two tutorial entries separately as `Tutorial` and `中文教程`, rather than hiding them inside one shared dropdown.

## Workspace Cleanliness

- After any local vignette, pkgdown, or package-check build, remove generated root-level artifacts that are not source files.
- Do not leave `doc/`, `Meta/`, `docs-site/`, local preview outputs, or vignette-generated `.html`, `.knit.md`, `.utf8.md`, or other compiled tutorial byproducts in the repository root after validation.
- Do not leave temporary local web servers or other long-running preview processes alive after checking pkgdown or tutorial pages.
- Keep only source directories and necessary configuration under the package root; generated verification artifacts should be rebuilt when needed, not stored as persistent workspace clutter.
- Before finishing a task, explicitly check that the package root is clean and that no build-only directories or preview processes remain.
