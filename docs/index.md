# PatchSim documentation

PatchSim is a configuration-driven framework for compartmental epidemiology on one
or more connected patches. The CLI loads YAML and CSV inputs, runs either
LSODA-based ODE integration or deterministic explicit Euler, and writes a
time-series CSV, a plot, and a run log.

PatchSim is under active development. Use the documentation for the installed
version, and validate every configuration before running it.

## Start here

- [Getting started](getting-started.md): create, validate, and run a real project.
- [Configuration](configuration.md): field names, file contracts, and expressions.
- [CLI reference](cli-reference.md): commands, output streams, and exit behavior.
- [Mathematical model](mathematical-model.md): transition and network equations.
- [Group stratification](group-stratification.md): custom groups and interaction matrices.
- [Contact generation](contact-generation.md): spatial kernels, units, and validation.
- [Simulation workflow](simulation-workflow.md): solver steps, Sobol analysis, and calibration examples.
- [Results](results.md): output paths, column naming, and rerun behavior.

## Documentation conventions

- Commands use `patchsim` as installed from a package. From a source checkout, use
  `uv run patchsim` instead.
- Relative paths in a configuration file are resolved from that file's directory,
  not from the shell's current directory.
- YAML field names and compartment names are case-sensitive unless a page says
  otherwise.
- Output columns such as `I_0` use the zero-based row order from `PatchFile`; the
  suffix is not the patch identifier.
- Examples describe behavior in the current code. Planned plugin and calibration
  features are not documented as available.

```{toctree}
:maxdepth: 2

getting-started.md
configuration.md
cli-reference.md
mathematical-model.md
rate-multiplication.md
network-design.md
group-stratification.md
contact-generation.md
simulation-workflow.md
results.md
architecture.md
```
