# Changelog

All notable changes to this project will be documented in this file.

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
