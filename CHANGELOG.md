# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-11-05

### Changed

- `solve!` now measures convergence with the relative L2 residual
  `‖r‖₂ / ‖b‖₂` instead of its square. `Param.tol` and the `res` value returned
  by `solve!` follow this convention. This is a breaking change for callers
  that relied on the previous squared-residual semantics: to reproduce the old
  stopping criterion, use `tol = old_tol^2` and compare `res^2`.
- The test suite now uses the L2 convention and covers the resulting iteration
  counts.

### Fixed

- Corrected the documented Helmholtz sign to `(∇² - coef) x = b`.
- Documented that the `msk` keyword cannot yet be used through the multilevel
  `setup` (see issue #17).

## [0.1.1] - 2026-11-05

### Added

- Completed the README with installation, quick-start, API, operator, smoother
  and boundary-condition documentation.

## [0.1.0] - 2024-07-22

### Added

- Initial release.

[0.2.0]: https://github.com/ClimFlows/CFMultigrid.jl/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/ClimFlows/CFMultigrid.jl/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/ClimFlows/CFMultigrid.jl/releases/tag/v0.1.0
