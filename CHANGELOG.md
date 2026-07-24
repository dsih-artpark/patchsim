# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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

### Changed
- **Breaking:** `NetworkModel.simulate_discrete()` now takes one Euler step per interval
  in `t_range`, using that interval's own width, and raises on divergence rather than
  returning negative populations. It previously assumed unit spacing regardless of the
  grid supplied, so any caller passing a non-unit `t_range` received silently incorrect
  results and will now see different values. Time grids must be finite and strictly
  increasing; unevenly spaced grids are supported.

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
