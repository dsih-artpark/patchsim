# Simulation workflow

## 1. Resolve configuration

`load_config` reads a YAML mapping and resolves `PatchFile`, `SeedFile`,
`NetworkFile`, and `OutputDir` relative to the config file.

## 2. Validate inputs

`setup_simulation` loads the CSV files and checks:

- positive patch populations;
- matching seed compartments and patch identifiers;
- non-negative seed values whose sums match patch populations;
- transition arrow syntax and known identifiers;
- patch-parameter identifiers; and
- day-zero network identifiers and non-negative weights.

`patchsim validate` stops after this setup stage.

## 3. Construct state

Patch order comes from `PatchFile`. Initial-state keys combine each compartment
with its zero-based patch index, such as `S_0` and `I_1`.

The network is a dense matrix in the same patch order. Only `day == 0` rows are
loaded.

## 4. Integrate

`patchsim run` creates reporting times with:

```python
np.arange(TMax, dtype=float)
```

The high-level model constructs the coupled ODE and solves it with
`scipy.integrate.odeint`. Internal LSODA steps are adaptive; output is sampled at
the reporting times.

## 5. Write artifacts

PatchSim writes:

- `runs/all_patches_MODEL_ode.csv`;
- `plots/patch_timeseries_MODEL_ode.png`; and
- `logs/MODEL_run_TIMESTAMP.log`.

The CSV and PNG are replaced by a rerun with the same `ModelName` and
`OutputDir`. See [Results](results.md) for retention conventions.

## Reproducible run checklist

Retain:

1. the exact YAML config;
2. every referenced CSV file;
3. the installed PatchSim version;
4. the generated CSV, PNG, and log; and
5. the command and environment used to run the model.

The current run log records Python and platform details but does not capture the
PatchSim package version or immutable copies of the inputs.
