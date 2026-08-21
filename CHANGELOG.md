# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - 2026-08-22

### Added

- Configuration-selectable adaptive ODE and deterministic explicit-Euler solvers that
  share reporting grids and derivative logic.
- Validated distance- and gravity-based spatial contact generation with explicit units,
  normalization settings, and a physical-validation report.
- Generic within-patch group stratification using modeller-supplied interaction matrices.
- Seeded first-order and total-order Sobol sensitivity studies with verified artifact
  reuse and compact provenance.
- Bounded multi-start least-squares calibration for prepared observations, including
  optional population-conserving fitted initial states and local diagnostics.
- A side-effect-free `patchsim.simulate()` API for external analysis and fitting code.
- SIS model template, complete workflow documentation, and the JOSS manuscript source.

### Changed

- SALib is available through the optional `analysis` dependency group rather than the
  base simulation installation.
- Run artifacts identify the selected solver, and the CLI exposes validation, contact
  generation, sensitivity, and calibration workflows.
- Documentation now states the implemented network orientation, transition semantics,
  solver behavior, validation limits, and artifact contracts.
- **Breaking:** `NetworkModel.simulate_discrete()` uses each interval's width and rejects
  invalid time grids or divergent steps instead of assuming unit spacing.
- **Breaking:** transition expressions reject floor division and modulo because they
  introduce discontinuities into rate arithmetic.

### Security

- Transition rate expressions are now parsed and evaluated against an arithmetic-only
  allowlist instead of `eval`. Expressions come from user-supplied YAML, and configs are
  shared as self-contained project directories, so a config could previously execute
  arbitrary code when a simulation was run.

### Fixed

- Rate expressions that overflow, divide by zero, or produce a complex value are now
  rejected instead of returning `inf` or a complex number into the simulation.
- Rate expressions accept numpy scalars of any dtype; previously only `float64` worked,
  because it is the one numpy type that subclasses Python's `float`.
- Patch-specific parameter overrides now affect CLI ODE runs instead of being parsed
  and ignored.
- Model validation rejects global or patch-specific parameter names that duplicate a
  compartment name.
- Explicit-Euler roundoff checks scale with each compartment update rather than unrelated
  population magnitudes.

## 0.1.0b1 - 2026-04-08

### Added
- Config-relative path resolution for `PatchFile`, `SeedFile`, `NetworkFile`, and `OutputDir`.
- Built-in YAML model templates for `sir`, `seir`, and `sirs` under package templates.
- `list-models` now reports both Python reference models and YAML templates.

### Changed
- Validation now checks transition source/target compartments explicitly.
- Validation now checks transition expression identifiers against `compartments ∪ Parameters`.
- Compartment source-of-truth reconciled: when `compartments` is provided in config, it must match `SeedFile` columns.
- Package metadata updated to `Development Status :: 4 - Beta`.
- Runtime dependencies now include conservative minimum version floors.

### Removed
- Removed unused `patchsim.utils.loader` module and its documentation reference.
