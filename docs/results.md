# Results

Each successful `patchsim run` writes a CSV, a PNG, and a timestamped log under
the configured `OutputDir`.

## Directory layout

For `ModelName: baseline` and `Solver: ode`:

```text
OutputDir/
  logs/
    baseline_run_YYYYMMDD_HHMMSS.log
  plots/
    patch_timeseries_baseline_ode.png
  runs/
    all_patches_baseline_ode.csv
```

The `OutputDir` path is resolved relative to the configuration file.

PatchSim does not create or require a `.patchsim` directory. All current run
artifacts are user-facing files.

## Retention and overwrite behavior

The CSV and PNG filenames are deterministic for each solver. Running the same
`ModelName` and solver into the same `OutputDir` overwrites those two files,
including runs with another `TimeStep`. The discrete solver uses `_discrete`
instead of `_ode`. Log filenames contain a timestamp, so logs from separate runs
are retained unless two runs start in the same second.

Use an output convention that identifies the scenario or retained run:

```yaml
OutputDir: output/baseline
```

```yaml
OutputDir: output/high-transmission
```

For an external orchestration system, generate a unique directory before
invoking PatchSim. Do not depend on PatchSim to version run artifacts.

## CSV columns

The CSV starts with `time`, followed by every compartment for every patch:

```csv
time,S_0,I_0,R_0,S_1,I_1,R_1
0.0,999.0,1.0,0.0,500.0,0.0,0.0
```

- `time` contains `k * TimeStep` for `k = 0` through `TMax - 1`.
- `COMPARTMENT_INDEX` identifies a compartment and zero-based patch index.
- Patch indices follow `PatchFile` row order.

For example, with:

```csv
patch,population
A,1000
B,500
```

`I_0` is infectious population in patch `A`, and `I_1` is infectious population
in patch `B`.

Grouped runs add the zero-based group index:

```csv
time,S_0_0,I_0_0,R_0_0,S_0_1,I_0_1,R_0_1
```

`COMPARTMENT_PATCH_GROUP` follows `PatchFile` and `GroupFile` ordering. The JSON
validation and run summaries include the ordered `groups` list.

Compartment values are floating point for both built-in solvers. The discrete
solver is deterministic explicit Euler, not an integer-state method.

## Plot

The PNG contains one subplot per patch. Ungrouped plots have one line per
compartment; grouped plots have one line per compartment/group combination.
Subplot titles use patch identifiers from `PatchFile`. Axes are numeric time and
compartment count.

The plot is a convenience view. Use the CSV for analysis and reproducibility.

## Log

The timestamped log records:

- model name and global parameters;
- resolved patch, seed, network, and output paths;
- Python and platform information;
- patch order;
- configured transitions;
- final CSV and PNG paths; and
- selected solver and `TimeStep`.

The log does not currently contain internal LSODA diagnostics, package-version
provenance, or a copy of the configuration. Archive the config and input CSVs
alongside retained outputs.

`Logging: false` does not disable this file in the current runner.

## JSON run summary

`patchsim run -c config.yaml --json` emits a summary with the resolved output
directory, CSV path, plot path, patch list, `TMax`, solver, and `TimeStep`. It
does not embed the time series. Logs continue to use standard error and the log
file.

## Basic analysis

```python
import pandas as pd

df = pd.read_csv("output/baseline/runs/all_patches_baseline_ode.csv")

peak_i0 = df["I_0"].max()
peak_time_i0 = df.loc[df["I_0"].idxmax(), "time"]
initial_total_0 = df.loc[0, ["S_0", "I_0", "R_0"]].sum()
final_total_0 = df.loc[df.index[-1], ["S_0", "I_0", "R_0"]].sum()

print(f"Patch 0 peak I: {peak_i0:.3f} at t={peak_time_i0:g}")
print(f"Patch 0 total change: {final_total_0 - initial_total_0:.3e}")
```

Interpret metrics in the context of the configured model. For example, final
`R` is not a general attack-rate measure for SIRS, SIS, or models with additional
flows.

## Sanity checks

Before using a result:

1. Confirm patch-index mapping against `PatchFile`.
2. Confirm the first CSV row matches `SeedFile`.
3. Check finite values and unexpected negative values.
4. Check conserved totals only if the transition structure should conserve them.
5. Retain the exact config and all referenced CSV files with the result.
