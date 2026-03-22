# LEO·UKB <img src="man/figures/logo.png?v=5" align="right" height="120" alt="leo.ukb logo" />

<!-- badges: start -->
[![build](https://img.shields.io/github/actions/workflow/status/laleoarrow/leo.ukb/pkgdown.yaml?branch=main&style=flat&label=build&labelColor=111827&logo=githubactions&logoColor=white)](https://github.com/laleoarrow/leo.ukb/actions/workflows/pkgdown.yaml)
[![docs](https://img.shields.io/badge/docs-github.io-14b8a6?style=flat&labelColor=111827&logo=githubpages&logoColor=white)](https://laleoarrow.github.io/leo.ukb/)
[![version 0.5.1](https://img.shields.io/badge/version-0.5.1-0ea5e9?style=flat&labelColor=111827&logo=git&logoColor=white)](https://github.com/laleoarrow/leo.ukb/releases)
[![dev experimental](https://img.shields.io/badge/dev-experimental-f59e0b?style=flat&labelColor=111827&logo=flask&logoColor=white)](https://github.com/laleoarrow/leo.ukb)
[![license Proprietary](https://img.shields.io/badge/license-Proprietary-64748b?style=flat&labelColor=111827&logo=bookstack&logoColor=white)](LICENSE)
[![R >= 3.5](https://img.shields.io/badge/R-%3E%3D%203.5-6366f1?style=flat&labelColor=111827&logo=r&logoColor=white)](https://cran.r-project.org/)
<!-- badges: end -->

`leo.ukb` is an internal toolkit for working with UK Biobank (UKB) data. It provides a small set of helper functions that streamline common workflows—field parsing/decoding, format conversion, batch processing, and preparing analysis-ready datasets—for large-scale biomedical studies.

> **Note**: This package is intended for academic research purposes. The code functionalities are experimental and subject to change.

<div align="center">
  <img src="man/figures/banner.png" alt="leo.ukb banner" width="700"/>
</div>

- [UK Biobank](https://www.ukbiobank.ac.uk/)
- [UKB Showcase](https://biobank.ndph.ox.ac.uk/showcase/)
- [UKB Schema](https://biobank.ndph.ox.ac.uk/showcase/schema.cgi)

## Installation

```r
# Install from GitHub (requires access permission)
remotes::install_github("laleoarrow/leo.ukb")
```

## Documentation

- Documentation home: [https://laleoarrow.github.io/leo.ukb/](https://laleoarrow.github.io/leo.ukb/)

## License

This software is proprietary. All rights reserved. See [LICENSE](LICENSE) for details.

## Authors

Ao Lu. Author, maintainer, copyright holder.

## Citation

Lu A (2026). leo.ukb: LEO package for UK Biobank database analysis. R package version 0.5.1, https://laleoarrow.github.io/leo.ukb/.

```bibtex
@Manual{leo.ukb,
  title = {leo.ukb: LEO package for UK Biobank database analysis},
  author = {Ao Lu},
  year = {2026},
  note = {R package version 0.5.1},
  url = {https://laleoarrow.github.io/leo.ukb/}
}
```
