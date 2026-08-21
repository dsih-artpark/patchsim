# CLI reference

## Invocation and streams

After installation:

```bash
patchsim COMMAND [OPTIONS]
```

From a source checkout:

```bash
uv run patchsim COMMAND [OPTIONS]
```

The CLI uses these conventions:

- Human-readable command results go to standard output.
- Logs and error diagnostics go to standard error and, for a run, to a file under
  `OutputDir/logs/`.
- `--json` writes the command result as JSON to standard output.
- A failed command exits non-zero and reports the underlying validation or runtime
  error.
- Paths inside a config file are resolved relative to the config file.

## Commands

| Command | Purpose |
| --- | --- |
| `init` | Create a project from a built-in template |
| `validate` | Load and validate a configuration and its input files |
| `run` | Run the configured solver and write artifacts |
| `sensitivity` | Run or reuse a configured Sobol study |
| `calibrate` | Run or reuse a bounded least-squares calibration |
| `list-models` | List built-in project templates |
| `generate-contacts` | Generate a validated spatial network CSV |

## `init`

```bash
patchsim init NAME [--template {seir,sir,sirs,sis}] [--force]
```

`NAME` is both the target directory and the generated `ModelName`. The default
template is `sir`.

```bash
patchsim init baseline
patchsim init latency-study --template seir
```

`--force` replaces a non-empty target directory after safety checks. It is
destructive: use it only when the existing directory can be discarded. PatchSim
refuses to replace the filesystem root, the current directory, or an ancestor of
the current directory.

The scaffold contains `config.yaml`, patch/network/seed CSV files, and an `output`
directory. It does not create a `.patchsim` directory.

## `validate`

```bash
patchsim validate -c CONFIG [--json]
patchsim validate --schema
```

`-c/--config` is required unless `--schema` is used.

Validation checks:

- required configuration fields;
- patch and population columns;
- positive patch populations;
- seed compartments, values, patch identifiers, and population totals;
- optional group populations, seed coverage, interaction labels, weights, and units;
- transition arrow syntax, compartment names, and expression identifiers;
- patch-parameter identifiers; and
- day-zero network patch identifiers and non-negative weights.

Plain output:

```text
Configuration is valid: config.yaml
```

JSON output includes `ok`, `config`, `model_name`, `num_patches`, and `patches`.
Grouped validation additionally includes `num_groups`, ordered `groups`, and
interaction diagnostics suitable for saving with a study.
When `Calibration` is present, validation reports `n`, `p`, starts, the maximum
forward simulations, and time-alignment or single-start warnings. A time warning
does not make the simulation configuration invalid, but `calibrate` will refuse
to run until the observations are aligned.
`--schema` prints the JSON Schema used for editor and tooling integration.

## `run`

```bash
patchsim run -c CONFIG [--json]
```

The command validates the same inputs, runs `Solver: ode` or
`Solver: discrete`, and writes:

| Path under `OutputDir` | Contents |
| --- | --- |
| `runs/all_patches_MODEL_SOLVER.csv` | Reporting time and every compartment |
| `plots/patch_timeseries_MODEL_SOLVER.png` | One subplot per patch |
| `logs/MODEL_run_YYYYMMDD_HHMMSS.log` | Resolved inputs and run details |

The JSON summary contains:

```json
{
  "ok": true,
  "config": "/path/to/config.yaml",
  "model_name": "baseline",
  "output_dir": "/path/to/output/baseline",
  "csv_path": "/path/to/output/baseline/runs/all_patches_baseline_ode.csv",
  "plot_path": "/path/to/output/baseline/plots/patch_timeseries_baseline_ode.png",
  "num_patches": 2,
  "patches": ["A", "B"],
  "t_max": 60,
  "solver": "ode",
  "time_step": 1.0
}
```

Resolved paths may be absolute. The JSON summary is not the simulation time
series.

:::{warning}
The CSV and PNG names are deterministic per solver. A rerun with the same
`ModelName`, solver, and `OutputDir` overwrites them, including a rerun with a
different `TimeStep`. Use one output directory per retained scenario or run.
Logs remain separate because their names include a timestamp.
:::

## `sensitivity`

```bash
patchsim sensitivity -c CONFIG [--json]
```

The command validates the simulation and `Sensitivity` block, reports the
planned evaluation count to standard error, and computes first-order and
total-order Sobol indices. It requires the `analysis` optional dependency.

The result summary contains `reused`, planned and completed evaluation counts,
elapsed seconds, and absolute paths to `samples.csv`, `responses.csv`,
`indices.csv`, and `manifest.json`. JSON mode keeps standard output valid JSON.

An existing named study is reused only when its request fingerprint and table
hashes match. A reused invocation reports zero completed evaluations. Changed
or incomplete studies fail rather than being overwritten. See the
{ref}`worked workflow <sensitivity-study>`.

## `calibrate`

```bash
patchsim calibrate -c CONFIG [--json]
```

The command validates the simulation and `Calibration` block, reports `n`, `p`,
the start count, warnings, and the maximum actual forward simulations to standard
error, then runs bounded deterministic TRF least squares from each declared
start. It uses existing SciPy dependencies.

The summary contains `reused`, `n`, `p`, the start count, the zero-based
selected-start index, actual forward simulations, elapsed seconds, warnings,
and absolute paths to `estimates.csv`, `fitted-seeds.csv`, `attempts.csv`,
`residuals.csv`, and `manifest.json`. JSON mode keeps standard output valid JSON.

One failed start does not discard successful starts. If every start fails, no
study directory is published. A matching named study is reused only after its
request and artifact hashes pass verification. See the {ref}`worked calibration
workflow <calibration-study>`.

## `list-models`

```bash
patchsim list-models [--json]
```

Plain output:

```text
Built-in models and templates:
- seir (yaml-template)
- sir (yaml-template)
- sirs (yaml-template)
- sis (yaml-template)
```

These are configuration templates, not separately implemented solver classes.

## `generate-contacts`

```bash
patchsim generate-contacts SOURCE OUTPUT \
  --id-column NAME \
  --kernel {distance,gravity} \
  --decay FLOAT \
  --min-distance-km FLOAT \
  --normalize {none,row}
```

Row normalization additionally requires `--self-share`. Unnormalized output
requires `--scale` and `--self-weight`. Gravity additionally requires
`--population-column`.

The command writes `OUTPUT` and `OUTPUT.validation.json`. It refuses to replace
either artifact unless `--force` is supplied. Vector input requires the `geo`
optional dependency, and polygon input requires `--centroid-crs`.

See [Contact generation](contact-generation.md) for formulas, units, complete
options, and validation guidance.

## Version

```bash
patchsim --version
```

The version is read from the installed package metadata.
