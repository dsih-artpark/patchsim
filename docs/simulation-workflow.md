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
np.arange(TMax, dtype=float) * TimeStep
```

Both built-in solvers use the same coupled derivative function. `Solver: ode`
uses `scipy.integrate.odeint`; internal LSODA steps are adaptive and output is
sampled at the reporting times. `Solver: discrete` takes one explicit-Euler step
per reporting interval.

## 5. Write artifacts

For `Solver: ode`, PatchSim writes:

- `runs/all_patches_MODEL_ode.csv`;
- `plots/patch_timeseries_MODEL_ode.png`; and
- `logs/MODEL_run_TIMESTAMP.log`.

The discrete solver replaces `_ode` with `_discrete` in the CSV and PNG names.
The log and JSON run summary record the selected solver and `TimeStep`.

The CSV and PNG are replaced by a rerun with the same `ModelName`, `OutputDir`,
and solver. See [Results](results.md) for retention conventions.

(sensitivity-study)=

## 6. Optional sensitivity study

Install the analysis dependency:

```bash
python -m pip install "patchsim[analysis]"
```

Add one `Sensitivity` block to the same simulation config:

```yaml
Sensitivity:
  Name: beta-gamma
  Method: sobol
  BaseSamples: 256
  Seed: 20260728
  Parameters:
    beta: [0.04, 0.12]
    gamma: [0.05, 0.20]
  Metrics:
    peak_infectious:
      Columns: [I_0, I_1]
      Reduce: max
    final_removed:
      Columns: [R_0, R_1]
      Reduce: final
```

The bounds define independent uniform input distributions. Both names must be
global `Parameters`; a name also present in `PatchParameters` is rejected.
Metric columns must exactly match the time-series columns for the configured
patch and group structure. PatchSim sums each metric's columns at every
reporting point, then takes `max` or `final`. A maximum is therefore a maximum
on the reporting grid, not a continuous-time optimum.

Validate before committing the compute budget:

```bash
patchsim validate -c config.yaml
```

With two parameters and `BaseSamples: 256`, the first-order/total-order Sobol
design requires:

```text
256 * (2 + 2) = 1024 model evaluations
```

Run the study:

```bash
patchsim sensitivity -c config.yaml
```

PatchSim writes:

```text
OutputDir/
  sensitivity/
    beta-gamma/
      samples.csv
      responses.csv
      indices.csv
      manifest.json
```

`indices.csv` contains `S1`, `ST`, and their seeded bootstrap confidence
interval half-widths for every metric/parameter pair. The intervals are the
estimate plus or minus the corresponding `_conf` value. `S1` measures the
parameter's first-order contribution; `ST` includes all interactions involving
that parameter. A material `ST - S1` can indicate aggregate interactions, but
this version does not attribute interactions to parameter pairs.

Repeat the study with `BaseSamples` doubled and a different `Name`. Compare
changes in the indices and confidence intervals; a power-of-two sample count
does not certify convergence. Sobol analysis also assumes the configured
independent bounds are scientifically meaningful. Correlated or jointly
constrained parameters need another method.

Reusing the same `Name` with the exact same config, inputs, method, and relevant
versions verifies the saved hashes and returns the existing artifacts without
new solves. PatchSim refuses an incomplete, modified, or different study at
that path. It does not create a hidden cache or overwrite a study.

(calibration-study)=

## 7. Optional calibration study

Prepare observations before giving them to PatchSim. The file must contain
numeric model-relative times, named observables, and values in the same units as
the corresponding model states:

```csv
time,observable,value
1,removed,4.0
2,removed,7.2
3,removed,9.76
4,removed,11.808
5,removed,13.4464
```

This synthetic example uses a one-patch discrete model with population `100`,
transition `"I -> R": "theta * I"`, and a deliberately low starting seed of
`I: 8`. The observations were generated with `theta = 0.2` and `I(0) = 20`.
Add this block to that simulation config:

```yaml
Calibration:
  Name: theta-i0
  Method: least_squares
  Observations: data/observations.csv
  MaxEvaluations: 200
  Observables:
    removed:
      Columns: [R_0]
      Scale: 1.0
  Parameters:
    theta: [0.01, 0.5]
  InitialConditions:
    - Patch: A
      Remainder: S
      Fit:
        I: [1, 50]
  Starts:
    - Parameters:
        theta: 0.4
      InitialConditions:
        - Patch: A
          Values:
            I: 40
```

The configured `theta` and seed `I` form start `0`; `Starts` adds a second
complete point. On every trial PatchSim sets `S(0) = 100 - I(0)`, preserving the
patch population. Initial states omitted from `Fit` remain fixed. If an initial
state was used to construct `SeedFile`, do not also include the same value as an
observation unless the state is fitted and the observation has an appropriate
residual scale.

Validate the data alignment and compute ceiling:

```bash
patchsim validate -c config.yaml
```

With five observations, two fitted variables, two starts, and
`MaxEvaluations: 200`, validation reports `n=5`, `p=2`, and at most `400`
forward simulations. This count includes simulations used for numerical
Jacobians. `n > p` does not by itself establish identifiability.

Run the fit:

```bash
patchsim calibrate -c config.yaml
```

The selected result should recover values close to the generating
`theta = 0.2` and `I(0) = 20`. Inspect every row of `attempts.csv`, the aligned
raw and standardized residuals, active bounds, and the local Jacobian
diagnostics before interpreting the estimate. Multiple starts reveal
sensitivity to initialization; they do not provide uncertainty intervals or a
global-optimum guarantee.

PatchSim writes:

```text
OutputDir/
  calibration/
    theta-i0/
      estimates.csv
      fitted-seeds.csv
      attempts.csv
      residuals.csv
      manifest.json
```

Reusing the same `Name` with the same config, observations, simulation inputs,
starts, method, and recorded versions verifies the artifacts and performs zero
new simulations. Use the fitted seed and global estimates for a separate normal
run; calibration does not overwrite the source config or seed file.

## Reproducible run checklist

Retain:

1. the exact YAML config;
2. every referenced CSV file;
3. the installed PatchSim version;
4. the generated CSV, PNG, and log; and
5. the command and environment used to run the model.

The current run log records Python and platform details but does not capture the
PatchSim package version or immutable copies of the inputs.

For sensitivity studies, `manifest.json` records the normalized configuration,
source and input hashes, method settings, versions, and artifact hashes. It
detects changed files but does not copy them, so retain the original YAML and
CSV inputs.

Calibration manifests add observation hashes, explicit starts, actual forward
simulation counts, selected-attempt termination, and local Jacobian diagnostics.
Retain the source observation file as well as the model inputs.
